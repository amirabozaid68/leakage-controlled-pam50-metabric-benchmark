# Changelog

## v1.2 (November 2026): current analysis

- **SCAN-B added.** GSE96058 (3,052 four-class tumours, RNA-seq) is downloaded by the notebook. The full benchmark is repeated inside SCAN-B, and every model is tested across cohorts in both directions.
- **Shared gene space.** Both cohorts are restricted to genes measured in both, so METABRIC numbers differ slightly from v1.1.
- **Tuning grid extended.** C is now chosen from 10^-4 to 3. In v1.1 the compact model chose the smallest grid value (0.01) in every split.
- **genefu probabilities.** For two SCAN-B tumours genefu returned no valid probability, and its crisp call is used for them (listed in `results_v1.2/genefu/genefu_undefined_probabilities.tsv`).

## v1.1 (November 2026): control arms

v1.1 lives in `archive/v1.1/`.

- **Canonical genefu PAM50 added.** It scored 0.836 macro-F1 in METABRIC against 0.712 for the v1.0 centroid re-implementation. Most of the v1.0 gain over PAM50 came from that re-implementation.
- **Two control arms added:** logistic regression on the 50 PAM50 genes, and on 3,000 genes. The PAM50-gene model was the best supervised model.
- **Panel genes recorded per split** in the primary four-class task.
- **Enrichment rerun** against the genes the selector could choose. The v1.0 enrichment used a whole-genome background.

## v1.0 (Zenodo 10.5281/zenodo.18983306): superseded

v1.0 lives in `archive/v1.0/`.

- **Stability and biology came from the wrong task.** The stability, overlap, enrichment and coefficient files in `results_main/` were produced by the secondary five-class run (Claudin-low included), not by the primary four-class task. The four-class robustness tables named in its manifest were not included.
- **The PAM50 comparator was a re-implementation**, not canonical genefu.
- **Its main conclusion does not hold** under the v1.1 and v1.2 controls.
