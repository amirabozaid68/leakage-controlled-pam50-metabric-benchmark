# Pre-specified analysis plan, v1.3

Written on 14 November 2026, before any v1.3 analysis was run. The v1.2 benchmark results (Zenodo 10.5281/zenodo.18983314) motivate the hypotheses below. No survival or clinical-marker data have been analysed for this project before this plan.

## Data

- **METABRIC.** 1,608 four-class tumours. Overall survival (OS) is censored at 120 months. Disease-specific survival (DSS, `VITAL_STATUS = Died of Disease`) is a secondary endpoint.
- **SCAN-B (GSE96058).** 3,052 four-class tumours. OS from the GEO characteristics. Ki67 status and Nottingham histological grade (NHG) are used where recorded.
- **Gene space.** 18,515 genes measured in both cohorts.
- **Model for H2 and H3.** Logistic regression on the 50 PAM50 genes, the best model in v1.2.

## H1. Subtype-balanced centring reduces the cross-cohort luminal shift

**Method.** Training data are centred on the mean of class means (the training labels give a balanced reference). The target cohort is centred without labels: start from its own mean, predict subtype probabilities, recompute the centre as the probability-weighted mean of class means, and repeat until it converges (at most 20 iterations). Each gene is scaled by the SD of its own cohort.

**Comparisons.** Per-cohort centring as in v1.2 ("naive"), balanced centring, and an oracle that balances with the true target labels. The oracle is reported only as a ceiling.

**Scope.** Models: LR on the 50 PAM50 genes (primary), LR on 3,000 genes, TOP80. Both directions, seeds 123 to 127, same source splits as v1.2.

**Success.** For the PAM50-gene model, balanced centring raises transfer macro-F1 in both directions, with a paired bootstrap interval above zero in at least 3 of 5 seeds, and it at least halves the dominant luminal error rate (LumB to LumA on METABRIC, LumA to LumB on SCAN-B).

## H2. Boundary-unstable luminal tumours form a prognostically intermediate group

**Population.** Tumours with a published LumA or LumB label.

**Four calls per tumour:**
- the published label;
- genefu PAM50 on the whole cohort (v1.2 calls);
- the out-of-fold call from 5-fold within-cohort CV of the PAM50-gene model;
- the call of the PAM50-gene model trained on the other cohort with naive centring.

**Groups:**
- stable LumA: all four calls are LumA;
- stable LumB: all four calls are LumB;
- unstable: all four calls are luminal but they disagree;
- non-luminal discordant: at least one call is Basal, Her2 or Normal. Reported, but not part of the test.

**Analysis.** Cox regression of OS on group, adjusted for age, with stable LumA as reference and then re-referenced to stable LumB. Kaplan-Meier curves. Mean proliferation score (mean z of the 11 PAM50 proliferation genes: BIRC5, CCNB1, CDC20, NUF2, CEP55, NDC80, MKI67, PTTG1, RRM2, TYMS, UBE2C). Ki67-high and NHG 3 fractions in SCAN-B.

**Success.** In both cohorts, the unstable group's HR lies between stable LumA and stable LumB, and its HR against stable LumA exceeds 1 with p < 0.05.

**Sensitivity.** Adding lymph node status. DSS in METABRIC.

## H3. A continuous luminal score carries prognostic information beyond the binary label

**Score.** P(LumB) / (P(LumA) + P(LumB)) among published luminal tumours. It comes from the out-of-fold within-cohort model, and from the other-cohort model with naive and with balanced centring. Each score is scaled to unit SD.

**Primary test.** In each cohort, a likelihood-ratio test of age + published LumB label + out-of-fold score against age + published LumB label (Cox, OS).

**Secondary tests:**
- the same comparison against the genefu binary call;
- Harrell's C-index of each score and each binary call;
- whether the score adds to the proliferation score;
- whether the other-cohort score is prognostic in the target cohort.

**Success.** The likelihood-ratio test gives p < 0.05 in both cohorts, and the other-cohort score is prognostic (p < 0.05, age-adjusted) in both target cohorts.

## Reporting

All pre-specified tests are reported whatever their result. Any analysis not listed here is labelled exploratory.
