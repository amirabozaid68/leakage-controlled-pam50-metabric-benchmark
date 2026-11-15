"""Builds Additional_file_1.tex (Figures S1-S3, Tables S1-S5).
Usage: python make_additional_file.py <ML_RESULTS_v1.3> <output folder containing FigS*.pdf> <posthoc results folder>"""
import sys, pandas as pd
R, OUTDIR, PH = sys.argv[1], sys.argv[2], sys.argv[3]
NM = {"lr_pam50genes": "LR, 50 PAM50 genes", "lr_full3000": "LR, 3,000 genes", "compact_TOP80": "LR, TOP80 panel"}
MD = {"naive": "Own-cohort", "balanced": "Balanced", "oracle": "Oracle (true labels)"}
CO = {"METABRIC": "METABRIC", "SCANB": "SCAN-B"}
EN = {"os_event": "Overall", "dss_event": "Disease-specific"}
def f(x, d=3): return "--" if pd.isna(x) else f"{x:.{d}f}"
def p(x): return "$<0.0001$" if x < 1e-4 else f"{x:.4f}" if x < 0.01 else f"{x:.3f}"

s = pd.read_csv(f"{R}/H1/H1_per_seed.tsv", sep="\t"); s = s[s.scope == "transfer"]
t = pd.read_csv(f"{R}/H1/H1_paired_balanced_minus_naive.tsv", sep="\t")
rows = []
for src in ["SCANB", "METABRIC"]:
    tgt = "METABRIC" if src == "SCANB" else "SCANB"
    for arm in NM:
        for mode in MD:
            q = s[(s.source == src) & (s.arm == arm) & (s["mode"] == mode)]
            extra = ""
            if mode == "balanced":
                tt = t[(t.source == src) & (t.arm == arm)]
                extra = f"{tt.delta.mean():+.3f} ({int((tt.ci_low > 0).sum())}/{int((tt.ci_high < 0).sum())})"
            rows.append(f"{CO[src]} $\\rightarrow$ {CO[tgt]} & {NM[arm]} & {MD[mode]} & {f(q.macro_f1.mean())} $\\pm$ "
                        f"{f(q.macro_f1.std())} & {f(q.LumA_to_LumB.mean(), 1)} & {f(q.LumB_to_LumA.mean(), 1)} & {extra} \\\\")
        rows.append("\\addlinespace")
S1 = "\n".join(rows)

h = pd.read_csv(f"{R}/H2/H2_cox.tsv", sep="\t")
GN = {"stable_LumA": "Stable LumA", "stable_LumB": "Stable LumB", "unstable": "Unstable"}
S2 = "\n".join(f"{CO[r.cohort]} & {EN[r.endpoint]} & {r.adjusted_for.replace('+node_pos', ' + node').replace('age', 'Age')} & "
               f"{GN[r.reference]} & {GN[r.group]} & {r.n} & {r.HR:.2f} ({r.HR_low:.2f}--{r.HR_high:.2f}) & {p(r.p)} \\\\"
               for r in h.itertuples())

lr = pd.read_csv(f"{R}/H3/H3_likelihood_ratio_tests.tsv", sep="\t")
AN = {"score_oof_z": "Within-cohort score", "score_cross_naive_z": "Other-cohort score (own centring)",
      "score_cross_balanced_z": "Other-cohort score (balanced)"}
BN = {"age+pub_LumB": "Age + published label", "age+genefu_LumB": "Age + genefu call", "age+prolif_z": "Age + proliferation"}
S3a = "\n".join(f"{CO[r.cohort]} & {EN[r.endpoint]} & {BN[r.base_model]} & {AN[r.added]} & {r.LR_stat:.2f} & {p(r.LR_p)} & "
                f"{r.HR_per_SD:.2f} ({r.HR_low:.2f}--{r.HR_high:.2f}) \\\\" for r in lr.itertuples())

c = pd.read_csv(f"{R}/H3/H3_concordance.tsv", sep="\t")
c = c[c.endpoint.isin(["os_event", "dss_event"])].drop(columns=[x for x in ["p", "n"] if x in c.columns])
u = pd.read_csv(f"{R}/H3/H3_age_adjusted_HR.tsv", sep="\t")
MN = {"pub_LumB": "Published LumB label", "genefu_LumB": "genefu LumB call", "score_oof_z": "Within-cohort score",
      "score_cross_naive_z": "Other-cohort score (own centring)", "score_cross_balanced_z": "Other-cohort score (balanced)",
      "prolif_z": "Proliferation score"}
cc = c.merge(u, on=["cohort", "endpoint", "marker"])
S3b = "\n".join(f"{CO[r.cohort]} & {EN[r.endpoint]} & {MN[r.marker]} & {r.c_index:.3f} & "
                f"{r.HR:.2f} ({r.HR_low:.2f}--{r.HR_high:.2f}) & {p(r.p)} \\\\" for r in cc.itertuples())

def table(cols, head, caption, body):
    return (f"\\begin{{longtable}}{{{cols}}}\n\\caption*{{{caption}}}\\\\\n\\toprule\n{head} \\\\\n\\midrule\\endfirsthead\n"
            f"\\toprule {head} \\\\ \\midrule\\endhead\n{body}\n\\bottomrule\n\\end{{longtable}}\n\n")

doc = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=2cm]{geometry}
\usepackage{graphicx,booktabs,longtable,caption,array}
\usepackage[T1]{fontenc}\usepackage{helvet}\renewcommand{\familydefault}{\sfdefault}
\captionsetup{font=small,labelfont=bf}
\begin{document}
\noindent{\Large\bfseries Additional file 1}\\[4pt]
\textbf{Data-driven gene panels do not improve PAM50-based breast cancer subtype classification: a leakage-controlled cross-cohort benchmark and analysis of Luminal A/Luminal B instability}\\[2pt]
Amira Abozaid, Mohamed Farouq Assal, Yasuhiro Kojima, Ahmed Osman, Ahmed B. Zaky

\vspace{1em}
\begin{figure}[h]\centering\includegraphics[width=\textwidth]{FigS1_all_metrics_within.pdf}
\caption*{\textbf{Figure S1. All performance measures within each cohort.} Mean $\pm$ SD over five stratified 70/15/15 splits for macro-F1, balanced accuracy, macro ROC-AUC and macro PR-AUC. Top row METABRIC, bottom row SCAN-B. The y-axis range differs between panels.}\end{figure}

\begin{figure}[h]\centering\includegraphics[width=\textwidth]{FigS2_per_class_recall.pdf}
\caption*{\textbf{Figure S2. Per-class recall within each cohort.} Mean recall over five test splits for each subtype and approach in (a) METABRIC and (b) SCAN-B.}\end{figure}

\begin{figure}[h]\centering\includegraphics[width=\textwidth]{FigS3_genefu_vs_published_labels.pdf}
\caption*{\textbf{Figure S3. genefu calls against the published PAM50 labels.} genefu PAM50 (pam50.robust) applied to every labelled tumour of (a) METABRIC and (b) SCAN-B. Numbers are tumour counts and shading is the row proportion.}\end{figure}
\clearpage
\scriptsize
\setlength{\tabcolsep}{3.5pt}
"""
doc += table("lllcccc", "Direction & Model & Centring & Macro-F1 & LumA$\\rightarrow$B (\\%) & LumB$\\rightarrow$A (\\%) & Balanced $-$ own",
             "\\textbf{Table S1. Cross-cohort transfer under three centring schemes (H1).} Mean $\\pm$ SD over five training splits. "
             "Luminal errors are the percentages of the target cohort's Luminal A tumours called Luminal B and of its Luminal B tumours "
             "called Luminal A. The last column gives the mean paired macro-F1 difference between balanced and own-cohort centring, with "
             "the number of splits whose bootstrap 95\\% interval lay above / below zero.", S1)
doc += table("llllllll", "Cohort & Endpoint & Adjusted for & Reference & Group & $n$ & HR (95\\% CI) & $p$",
             "\\textbf{Table S2. Cox models for boundary-unstable luminal tumours (H2).} Hazard ratios (95\\% CI) for each group against "
             "the reference group among published luminal tumours. Tumours with a non-luminal call are excluded. METABRIC follow-up is "
             "censored at 10 years.", S2)
doc += table("lllllll", "Cohort & Endpoint & Base model & Added & LR statistic & $p$ & HR per SD (95\\% CI)",
             "\\textbf{Table S3a. Does the continuous luminal score add prognostic information (H3)?} Likelihood-ratio tests of the base "
             "Cox model against the base model plus the score (1 degree of freedom), among published luminal tumours. The pre-specified "
             "primary test is age + published label + within-cohort score for overall survival.", S3a)
doc += table("llllll", "Cohort & Endpoint & Marker & C-index & HR (95\\% CI) & $p$",
             "\\textbf{Table S3b. Discrimination and age-adjusted hazard ratios of each marker.} Harrell's C-index (higher marker value "
             "predicting shorter survival) and age-adjusted Cox hazard ratio among published luminal tumours. Scores and the "
             "proliferation score are per SD.", S3b)
ld = pd.read_csv(f"{PH}/H2_cox_ladder.tsv", sep="\t")
CV = {"age": "age", "node_pos": "node", "grade3": "grade 3", "log_size": "size", "her2": "HER2", "pr": "PR", "postmeno": "menopause",
      "chemo": "chemotherapy", "endocrine": "endocrine", "radio": "radiotherapy"}
def cov(s): return " + ".join(CV[x] for x in s.split(",")[0].split("+")) + (" (strata)" if "stratified" in s else "")
S4 = "\n".join(f"{CO[r.cohort]} & {r.model} & {cov(r.covariates)} & {r.sample} & {r.n} ({r.events}) & {r.n_U} ({r.events_U}) & "
                f"{r.HR_U:.2f} ({r.lo_U:.2f}--{r.hi_U:.2f}) & {p(r.p_U)} & {r.HR_B:.2f} ({r.lo_B:.2f}--{r.hi_B:.2f}) & {r.PH_p_U:.2f} \\\\"
                for r in ld.sort_values(["cohort", "model", "sample"]).itertuples())
src = pd.read_csv(f"{PH}/H2_instability_sources.tsv", sep="\t")
S5 = "\n".join(f"{CO[r.cohort]} & {r.source} & {r.n_A} & {r.n_U} ({r.events_U}) & {r.n_B} & {r.excluded_nonluminal} & "
                f"{r.HR_U:.2f} ({r.lo_U:.2f}--{r.hi_U:.2f}) & {p(r.p_U)} & {r.HR_B:.2f} ({r.lo_B:.2f}--{r.hi_B:.2f}) \\\\"
                for r in src.itertuples())
doc += table("ll>{\\raggedright\\arraybackslash}p{2.7cm}>{\\raggedright\\arraybackslash}p{1.4cm}llllll", "Cohort & Model & Covariates & Sample & $n$ (events) & Unstable (events) & HR unstable & $p$ & HR stable B & PH $p$",
             "\\textbf{Table S4. Clinically adjusted Cox models for boundary-unstable tumours (post hoc).} Hazard ratios against stable "
             "Luminal A. METABRIC: breast cancer-specific survival censored at 10 years. SCAN-B: overall survival. M4 is stratified by "
             "METABRIC sub-cohort. Missing covariates were handled by complete-case analysis, and models M1 to M3 are also shown on the "
             "M4 complete cases. PH $p$: scaled Schoenfeld residual test for the unstable group. In SCAN-B the age term departed from "
             "proportional hazards ($p = 0.006$ to 0.016).", S4)
doc += table("l>{\\raggedright\\arraybackslash}p{4.0cm}lllllll", "Cohort & Source of instability & Stable A & Unstable (events) & Stable B & Excluded & HR unstable & $p$ & HR stable B",
             "\\textbf{Table S5. Survival by source of Luminal A/Luminal B instability (post hoc).} Unstable tumours disagree on the listed "
             "calls, stable groups agree on them. Excluded: tumours with a non-luminal call among the listed calls. Age-adjusted hazard "
             "ratios against stable Luminal A. METABRIC: breast cancer-specific survival. SCAN-B: overall survival.", S5)
doc += "\\end{document}\n"
open(f"{OUTDIR}/Additional_file_1.tex", "w").write(doc)
print("written")
