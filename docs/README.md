# METABRIC PAM50 Benchmark Release v1

## Project title
Leakage-Controlled PAM50 Benchmarking in METABRIC Reveals a Stable Compact Gene Panel for Intrinsic Subtype Classification

## What the release contains
This release contains the code, processed inputs, main analysis outputs, secondary analysis outputs, manuscript figures, and reproducibility documentation for a leakage-controlled benchmarking workflow for breast cancer intrinsic subtype classification in METABRIC.

The project started in **R** with cohort reconstruction, clinical harmonization, PAM50 label definition, and export of aligned processed inputs, then moved to **Google Colab / Python** for leakage-controlled benchmarking, robustness analyses, compact panel extraction, PAM50 comparison, interpretability analyses, and figure generation.

Release structure:

- `code/R/`
  R code used for cohort reconstruction and export of processed inputs.
- `code/python/`
  Colab notebook or Python script used for benchmarking and figure generation.
- `code/environment/`
  Software environment records for R and Python.
- `data_processed/`
  Processed aligned inputs used by the benchmarking workflow.
- `results_main/`
  Main outputs used in the manuscript primary claim.
- `results_secondary/`
  Secondary and exploratory outputs.
- `figures_main/`
  Main manuscript figures.
- `figures_supplementary/`
  Supplementary figures.
- `docs/`
  Documentation for run order, checksums, and release notes.

## How the workflow runs
The workflow is split into two stages.

First, an **R Markdown** workflow reconstructs the METABRIC cohort from cBioPortal raw files, performs strict `SAMPLE_ID` matching, harmonizes PAM50 labels, and exports aligned processed files into `EXPORT_FOR_COLAB/`.

Second, a **Colab Python** workflow takes these processed inputs, performs leakage-controlled preprocessing and benchmarking, exports PAM50 centroids from `genefu`, evaluates PAM50 and compact logistic-regression models on identical splits, quantifies robustness across seeds, derives a stable compact panel, performs overlap and enrichment analyses, and produces manuscript-ready figures.

## Input files
The main processed input files are:

- `X_samples_x_genes.tsv.gz`
  Processed expression matrix with rows as samples and columns as genes.
- `y_labels.tsv`
  Sample-level labels used for downstream classification tasks.
- `meta.tsv`
  Sample-level metadata including `SAMPLE_ID`, patient identifiers, and label annotations.

These files are created by the R stage and used as the input to the Colab benchmarking stage.

## Run order
Step 1: Run `Transition_to_python.Rmd`
Output: `EXPORT_FOR_COLAB/`

Step 2: Run `METABRIC_PAM50_benchmark.ipynb` or `.py` in Colab
Input: `EXPORT_FOR_COLAB/`
Output: `ML_RESULTS_UPGRADED/`

Step 3: Use the files in `results_main/` and `figures_main/` for the manuscript

## Main outputs used in manuscript
Examples of key manuscript outputs include:

- `cohort_flow_table.tsv`
- `pam50_centroids.tsv`
- `pam50_genes.tsv`
- `pam50_eval_primary4_val_test.csv`
- `compact_eval_primary4_val_test_TOP80.csv`
- `compact_lr_primary4_aucpr_val_test_TOP80.csv`
- `robustness_pam50_vs_compactLR_4class_per_seed.csv`
- `robustness_pam50_vs_compactLR_4class_summary.csv`
- `stable_compact_panel_TOP80_5seeds_5.csv`
- `TOP80_overlap_with_PAM50.csv`
- `TOP80_nonoverlap_vs_PAM50.csv`
- `enrichment_TOP80_nonoverlap_enrichr.csv`
- `STABLE_TOP80_multinomialLR_coefficients_4class_seed123split.csv`
- `STABLE_TOP80_top15_genes_per_class_4class_seed123split.csv`

## How figures map to manuscript figure numbers
Current figure-file mapping is as follows:

- **Figure 1** → `fig1_workflow.pdf`
- **Figure 2** → `fig2_robustness.pdf`
- **Figure 3** → `fig3_confusion.pdf`
- **Figure 4** → `fig4_coefficients_Heatmap.pdf`

If PNG versions are also included in the release, they should use the same base names.

## Notes on manuscript assembly
The current LaTeX manuscript uses the following figure files directly:

- `fig1_workflow.pdf`
- `fig2_robustness.pdf`
- `fig3_confusion.pdf`
- `fig4_coefficients_Heatmap.pdf`

Make sure these filenames are preserved exactly during release packaging and manuscript submission.

## Contact
**Amira Abozaid**
Biotechnology Program, Institute of Basic and Applied Sciences, Egypt-Japan University of Science and Technology (E-JUST)
Molecular Biology Division, Zoology Department, Faculty of Science, Alexandria University
Email: Amira.abozaid@ejust.edu.eg

For questions about the release, code, or manuscript assembly, contact the corresponding author.
