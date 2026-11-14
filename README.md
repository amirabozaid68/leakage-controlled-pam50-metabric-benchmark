# Leakage-controlled PAM50 benchmark: METABRIC and SCAN-B (release v1.2)

Code and results for the manuscript *Data-driven gene panels do not outperform the PAM50 genes for breast cancer intrinsic subtyping: a leakage-controlled benchmark across microarray and RNA-seq cohorts* (Abozaid A, Assal MF, Kojima Y, Osman A, Zaky AB).

Release v1.2 replaces the main analysis of v1.0. See `CHANGELOG.md` for what changed and why. The earlier releases are kept unchanged in `archive/`.

## What the analysis does

Four-class intrinsic subtype classification (Basal-like, HER2-enriched, Luminal A, Luminal B) in two cohorts:

| Cohort | Platform | Four-class tumours | Source |
|---|---|---|---|
| METABRIC | Illumina HT-12 microarray | 1,608 | cBioPortal |
| SCAN-B | RNA-seq | 3,052 | GEO GSE96058 (downloaded by the notebook) |

Both cohorts are restricted to genes measured in both. Approaches compared on identical splits:

1. genefu PAM50 (`intrinsic.cluster.predict`, `pam50.robust`), run on each test set alone
2. logistic regression on the 50 PAM50 genes
3. logistic regression on 3,000 training-set variance-filtered genes
4. L1-selected compact panels of 20, 40, 50, 80 (primary) and 120 genes, followed by logistic regression
5. the v1.0 centroid re-implementation, kept as a historical reference

Within each cohort: five stratified 70/15/15 splits (seeds 123 to 127), all preprocessing fitted on training data only, C tuned on validation macro-F1, one test pass after refitting on training plus validation. Across cohorts: every model refitted in one cohort predicts every four-class tumour of the other, in both directions, with genefu on the whole target cohort as reference.

## How to reproduce

1. Run `code/R/Transition_to_python_BeatPAM50.Rmd` in RStudio (Knit) on the cBioPortal METABRIC download. It writes `EXPORT_FOR_COLAB/` with `X_samples_x_genes.tsv.gz`, `y_labels.tsv` and `meta.tsv`. The expression matrix is not in this release because of its size.
2. Put `EXPORT_FOR_COLAB/` in Google Drive (`MyDrive/EXPORT_FOR_COLAB`).
3. Open `code/python/PAM50_benchmark_v1.2.ipynb` in Google Colab and choose Runtime > Run all. The notebook installs R and genefu, downloads GSE96058 into Drive, and writes `ML_RESULTS_v1.2/`. Heavy steps are checkpointed, so Run all can be repeated after a disconnect.
4. Run `code/python/make_figures_v1.2.py ML_RESULTS_v1.2 figures_v1.2` to draw the figures.

## Folder contents

| Folder | Contents |
|---|---|
| `code/python/` | v1.2 Colab notebook and figure script |
| `code/R/` | METABRIC cohort reconstruction and export (unchanged from v1.0) |
| `code/environment/` | run configuration, Python package versions, R session information |
| `data_processed/` | METABRIC labels, metadata and cohort flow table |
| `results_v1.2/` | all v1.2 outputs (see below) |
| `figures_v1.2/` | manuscript figures, PDF and 600-dpi PNG |
| `archive/v1.0/`, `archive/v1.1/` | earlier analyses, unchanged |

Key files in `results_v1.2/`:

| File | Used for |
|---|---|
| `summary_mean_sd.tsv` | Tables 2 and 3 (mean and SD by approach, cohort and scope) |
| `per_seed/per_seed_metrics_all.tsv` | per-split metrics, Figures 2 and 3 |
| `per_seed/predictions_all.tsv.gz` | per-tumour predictions and probabilities for every approach |
| `per_seed/tuned_hyperparameters.tsv` | chosen C per approach and split |
| `paired_tests.tsv` | paired bootstrap intervals and exact McNemar tests |
| `luminal_confusions.tsv` | Luminal A/Luminal B reassignment counts, Figure 4 |
| `genefu/` | genefu calls, scripts and agreement with published labels, Supplementary Figure S3 |
| `stability/` | panel genes per split, inclusion frequencies, stable panels, Figure 5 |
| `biology/` | stable-panel coefficients, enrichment results and the gene-set libraries used |
| `splits_by_sample_id.tsv` | every train/validation/test assignment by sample ID |

## Figure mapping

| Manuscript | File |
|---|---|
| Figure 1 | `figures_v1.2/Fig1_study_design.pdf` |
| Figure 2 | `figures_v1.2/Fig2_within_cohort_macroF1.pdf` |
| Figure 3 | `figures_v1.2/Fig3_cross_cohort_transfer.pdf` |
| Figure 4 | `figures_v1.2/Fig4_luminal_boundary.pdf` |
| Figure 5 | `figures_v1.2/Fig5_panel_stability_and_size.pdf` |
| Supplementary Figures S1 to S3 | `figures_v1.2/FigS1_*`, `FigS2_*`, `FigS3_*` |

## Licence and contact

Code is released under the MIT License. Contact: Amira Abozaid, Egypt-Japan University of Science and Technology (E-JUST).
