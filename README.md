# Leakage-controlled PAM50 benchmark: METABRIC and SCAN-B (release v1.3)

Code and results for the manuscript *Gene panels do not improve on PAM50, and the Luminal A/Luminal B boundary is not reproducible across cohorts: a leakage-controlled benchmark in METABRIC and SCAN-B* (Abozaid A, Assal MF, Kojima Y, Osman A, Zaky AB).

Release v1.3 adds the pre-specified luminal-boundary analyses (H1 to H3) to the v1.2 benchmark. The plan for those analyses, `ANALYSIS_PLAN_v1.3.md`, was committed to this repository before they were run. See `CHANGELOG.md` for what changed in each release. The earlier analyses are kept unchanged in `archive/`.

## What the analysis does

Four-class intrinsic subtype classification (Basal-like, HER2-enriched, Luminal A, Luminal B) in two cohorts:

| Cohort | Platform | Four-class tumours | Survival | Source |
|---|---|---|---|---|
| METABRIC | Illumina HT-12 microarray | 1,608 | overall and disease-specific | cBioPortal |
| SCAN-B | RNA-seq | 3,052 | overall | GEO GSE96058 (downloaded by the notebook) |

Both cohorts are restricted to the 18,515 genes measured in both.

**Benchmark (v1.2).** These approaches are compared on identical splits:

1. genefu PAM50 (`intrinsic.cluster.predict`, `pam50.robust`), run on each test set alone
2. logistic regression on the 50 PAM50 genes
3. logistic regression on 3,000 training-set variance-filtered genes
4. L1-selected compact panels of 20, 40, 50, 80 (primary) and 120 genes, followed by logistic regression
5. the v1.0 centroid re-implementation, kept as a historical reference

Within each cohort there are five stratified 70/15/15 splits (seeds 123 to 127), with all preprocessing fitted on training data only. Across cohorts, every model refitted in one cohort predicts every four-class tumour of the other, in both directions.

**Pre-specified luminal analyses (v1.3):**

- **H1.** Does subtype-balanced centring repair the cross-cohort Luminal A/Luminal B shift?
- **H2.** Do boundary-unstable luminal tumours (conflicting calls from the published label, genefu, a within-cohort model and an other-cohort model) have intermediate survival?
- **H3.** Does a continuous Luminal B probability add prognostic information beyond the binary label?

`results_v1.3/PRESPECIFIED_VERDICTS.tsv` lists the pass or fail outcome of each test.

## How to reproduce

1. Run `code/R/Transition_to_python_BeatPAM50.Rmd` in RStudio (Knit) on the cBioPortal METABRIC download. It writes `EXPORT_FOR_COLAB/` with `X_samples_x_genes.tsv.gz`, `y_labels.tsv` and `meta.tsv`. The expression matrix is not in this release because of its size.
2. Put `EXPORT_FOR_COLAB/` in Google Drive (`MyDrive/EXPORT_FOR_COLAB`), together with the cBioPortal file `data_clinical_patient.txt`.
3. In Google Colab, run `code/python/PAM50_benchmark_v1.2.ipynb` (Runtime > Run all). It installs R and genefu, downloads GSE96058 and writes `ML_RESULTS_v1.2/`.
4. Then run `code/python/PAM50_benchmark_v1.3.ipynb`, which writes `ML_RESULTS_v1.3/`. Both notebooks checkpoint heavy steps.
5. Draw the figures: `python code/python/make_figures.py ML_RESULTS_v1.2 figures ML_RESULTS_v1.3`
6. Run the post hoc analyses (v1.3.1): `python code/python/posthoc_v131.py ML_RESULTS_v1.2 ML_RESULTS_v1.3 <folder with the cBioPortal files data_clinical_patient.txt and data_clinical_sample.txt> results_v1.3/posthoc`
7. Build Additional file 1: `python code/python/make_additional_file.py ML_RESULTS_v1.3 figures results_v1.3/posthoc`

## Folder contents

| Folder | Contents |
|---|---|
| `code/python/` | v1.2 and v1.3 Colab notebooks, post hoc script, figure script, Additional file script |
| `code/R/` | METABRIC cohort reconstruction and export |
| `code/environment/` | run configuration, Python package versions, R session information |
| `data_processed/` | METABRIC labels, metadata and cohort flow table |
| `results_v1.2/` | benchmark outputs |
| `results_v1.3/` | H1 to H3 outputs, tumour-level scores and the pre-specified verdicts |
| `results_v1.3/posthoc/` | post hoc analyses added in v1.3.1 (not pre-specified) |
| `figures_v1.3/` | all manuscript and supplementary figures, PDF and 600-dpi PNG |
| `figures_v1.2/` | figures as released in v1.2 |
| `archive/v1.0/`, `archive/v1.1/` | earlier analyses, unchanged |

## Figure mapping

| Manuscript | File in `figures_v1.3/` | Source results |
|---|---|---|
| Figure 1 | `Fig1_study_design.pdf` | design schematic |
| Figure 2 | `Fig2_within_cohort_macroF1.pdf` | `results_v1.2/per_seed/` |
| Figure 3 | `Fig3_cross_cohort_transfer.pdf` | `results_v1.2/per_seed/` |
| Figure 4 | `Fig4_luminal_boundary.pdf` | `results_v1.2/per_seed/predictions_all.tsv.gz` |
| Figure 5 | `Fig5_panel_stability_and_size.pdf` | `results_v1.2/stability/` |
| Figure 6 | `Fig6_centring_and_luminal_shift.pdf` | `results_v1.3/H1/` |
| Figure 7 | `Fig7_boundary_unstable_luminal.pdf` | `results_v1.3/H2/` |
| Additional file 1, Figures S1 to S3 | `FigS1_*`, `FigS2_*`, `FigS3_*` | `results_v1.2/` |
| Additional file 1, Tables S1 to S3 | built by `make_additional_file.py` | `results_v1.3/H1/`, `H2/`, `H3/` |
| Additional file 1, Tables S4 and S5 | built by `make_additional_file.py` | `results_v1.3/posthoc/` |

## Licence and contact

Code is released under the MIT License. Contact: Amira Abozaid, Egypt-Japan University of Science and Technology (E-JUST).
