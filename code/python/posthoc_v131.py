"""Post hoc analyses added after pre-submission review (v1.3.1).
Usage: python posthoc_v131.py <ML_RESULTS_v1.2> <ML_RESULTS_v1.3> <METABRIC clinical dir> <out dir>
Not part of the pre-specified plan."""
import sys, warnings, numpy as np, pandas as pd
from sklearn.metrics import f1_score
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
warnings.filterwarnings('ignore')
r12, r13, clin, out = sys.argv[1:5]
L = {'LumA', 'LumB'}

# 1. Pooled patient-level bootstrap of the principal within-cohort comparisons
p = pd.read_csv(f'{r12}/per_seed/predictions_all.tsv.gz', sep='\t')
rng = np.random.default_rng(2026); rows = []
for scope in ['METABRIC_within', 'SCANB_within']:
    q = p[(p.scope == scope) & (p.seed != 0)]
    w = q.pivot_table(index=['seed', 'SAMPLE_ID', 'y_true'], columns='arm', values='y_pred').reset_index()
    ids = w.SAMPLE_ID.unique(); pos = w.groupby('SAMPLE_ID').indices
    for ref, role in [('compact_TOP80', 'primary'), ('lr_full3000', 'secondary'), ('pam50_genefu', 'secondary')]:
        d = lambda x: f1_score(x.y_true, x.lr_pam50genes, average='macro') - f1_score(x.y_true, x[ref], average='macro')
        bs = np.array([d(w.iloc[np.concatenate([pos[i] for i in rng.choice(ids, len(ids))])]) for _ in range(2000)])
        rows.append(dict(cohort=scope.split('_')[0], comparison=f'lr_pam50genes minus {ref}', role=role, patients=len(ids),
                         test_predictions=len(w), diff=d(w), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                         p=max(2 * min((bs <= 0).mean(), (bs >= 0).mean()), 1 / 2000)))
b = pd.DataFrame(rows)
for c, g in b[b.role == 'secondary'].groupby('cohort'):  # Holm across secondary comparisons per cohort
    o = g.p.sort_values(); adj = np.minimum(1, np.maximum.accumulate(o.values * (len(o) - np.arange(len(o)))))
    b.loc[o.index, 'p_holm'] = adj
b.to_csv(f'{out}/pooled_bootstrap_within.tsv', sep='\t', index=False)

# 2. Survival data with the four-call groups
def load(cohort):
    d = pd.read_csv(f'{r13}/tumour_level_scores_{cohort}.tsv.gz', sep='\t')
    d = d[d.published.isin(L)].copy()
    if cohort == 'METABRIC':
        d['T'] = d.os_time.clip(upper=120); d['E'] = ((d.dss_event == 1) & (d.os_time <= 120)).astype(int)
    else:
        d['T'] = d.os_time; d['E'] = d.os_event
    return d

def grouped(d, calls):
    d = d[d[calls].isin(L).all(axis=1)].copy()
    a = (d[calls] == 'LumA').all(axis=1); bb = (d[calls] == 'LumB').all(axis=1)
    d['U'] = (~a & ~bb).astype(int); d['B'] = bb.astype(int)
    return d

def cox(d, cov, strata=None):
    cols = ['T', 'E', 'U', 'B'] + cov + ([strata] if strata else [])
    y = d[cols].dropna()
    m = CoxPHFitter().fit(y, 'T', 'E', strata=[strata] if strata else None); s = m.summary
    ph = proportional_hazard_test(m, y, time_transform='rank').summary
    r = dict(n=len(y), n_A=int(((y.U == 0) & (y.B == 0)).sum()), n_U=int(y.U.sum()), n_B=int(y.B.sum()),
             events=int(y.E.sum()), events_U=int(y.loc[y.U == 1, 'E'].sum()))
    for g in ['U', 'B']:
        r.update({f'HR_{g}': s.loc[g, 'exp(coef)'], f'lo_{g}': s.loc[g, 'exp(coef) lower 95%'],
                  f'hi_{g}': s.loc[g, 'exp(coef) upper 95%'], f'p_{g}': s.loc[g, 'p']})
    r.update(PH_p_U=ph.loc['U', 'p'], PH_min_p_any_term=ph.p.min())
    return r

FOUR = ['published', 'genefu', 'oof', 'cross_naive']
SOURCES = {'overall (four calls)': FOUR, 'implementation (published, genefu)': ['published', 'genefu'],
           'model (published, within-cohort model)': ['published', 'oof'],
           'transfer (within-cohort, other-cohort model)': ['oof', 'cross_naive']}
rows = []
for cohort, ep in [('METABRIC', 'DSS'), ('SCANB', 'OS')]:
    d = load(cohort)
    for name, calls in SOURCES.items():
        rows.append(dict(cohort=cohort, endpoint=ep, source=name, excluded_nonluminal=int((~d[calls].isin(L).all(axis=1)).sum()),
                         **cox(grouped(d, calls), ['age'])))
pd.DataFrame(rows).to_csv(f'{out}/H2_instability_sources.tsv', sep='\t', index=False)

# 3. Cox adjustment ladder
pt = pd.read_csv(f'{clin}/data_clinical_patient.txt', sep='\t', comment='#')
sp = pd.read_csv(f'{clin}/data_clinical_sample.txt', sep='\t', comment='#')
c = sp.merge(pt, on='PATIENT_ID'); yn = {'YES': 1, 'NO': 0}; pn = {'Positive': 1, 'Negative': 0}
c = c.assign(grade3=(c.GRADE == 3).astype(float).where(c.GRADE.notna()), log_size=np.log(c.TUMOR_SIZE),
             her2=c.HER2_STATUS.map(pn), pr=c.PR_STATUS.map(pn), postmeno=c.INFERRED_MENOPAUSAL_STATE.map({'Post': 1, 'Pre': 0}),
             chemo=c.CHEMOTHERAPY.map(yn), endocrine=c.HORMONE_THERAPY.map(yn), radio=c.RADIO_THERAPY.map(yn), cohort=c.COHORT)
extra = ['grade3', 'log_size', 'her2', 'pr', 'postmeno', 'chemo', 'endocrine', 'radio', 'cohort']
m = grouped(load('METABRIC'), FOUR).merge(c[['SAMPLE_ID'] + extra], on='SAMPLE_ID', how='left')
s = grouped(load('SCANB'), FOUR); s['grade3'] = (s.nhg == 3).astype(float).where(s.nhg.notna())
LAD = {'M1': ['age'], 'M2': ['age', 'node_pos'], 'M3': ['age', 'node_pos', 'grade3'],
       'M4': ['age', 'node_pos', 'grade3', 'log_size', 'her2', 'pr', 'postmeno', 'chemo', 'endocrine', 'radio']}
cc = m.dropna(subset=LAD['M4'] + ['cohort']); rows = []
for k, cov in LAD.items():
    st = 'cohort' if k == 'M4' else None
    rows.append(dict(cohort='METABRIC', endpoint='DSS', model=k, sample='all available', covariates='+'.join(cov) + (', stratified by METABRIC sub-cohort' if st else ''), **cox(m, cov, st)))
    if k != 'M4':
        rows.append(dict(cohort='METABRIC', endpoint='DSS', model=k, sample='M4 complete cases', covariates='+'.join(cov), **cox(cc, cov)))
        rows.append(dict(cohort='SCANB', endpoint='OS', model=k, sample='all available', covariates='+'.join(cov), **cox(s, cov)))
pd.DataFrame(rows).to_csv(f'{out}/H2_cox_ladder.tsv', sep='\t', index=False)
m[LAD['M4'] + ['cohort']].isna().mean().rename('fraction_missing').to_csv(f'{out}/H2_covariate_missingness_METABRIC.tsv', sep='\t')
print('done')
