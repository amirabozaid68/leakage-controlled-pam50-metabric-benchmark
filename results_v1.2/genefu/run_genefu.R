
suppressPackageStartupMessages(library(genefu))
suppressWarnings(try(data(pam50.robust, package = "genefu"), silent = TRUE))
if (!exists("pam50.robust")) data(pam50, package = "genefu")
d <- Sys.getenv("GDIR")
tests <- read.delim(file.path(d, "test_sets.tsv"))
run <- function(M, tag, cohort) {
  annot <- data.frame(probe = colnames(M), EntrezGene.ID = NA, row.names = colnames(M))
  r <- intrinsic.cluster.predict(sbt.model = pam50.robust, data = M, annot = annot,
                                 do.mapping = FALSE, do.prediction.strength = FALSE, verbose = FALSE)
  data.frame(cohort = cohort, run = tag, SAMPLE_ID = rownames(M), call = as.character(r$subtype),
             r$subtype.proba, check.names = FALSE)
}
res <- NULL
for (co in c("METABRIC", "SCANB")) {
  X <- read.delim(file.path(d, paste0("pam50_expr_", co, ".tsv")), check.names = FALSE)
  rownames(X) <- X$SAMPLE_ID; X$SAMPLE_ID <- NULL; X <- as.matrix(X)
  stopifnot(all(rownames(pam50.robust$centroids) %in% colnames(X)))
  res <- rbind(res, run(X, "full_cohort", co))
  for (s in unique(tests$seed[tests$cohort == co])) {
    ids <- as.character(tests$SAMPLE_ID[tests$cohort == co & tests$seed == s])
    res <- rbind(res, run(X[ids, , drop = FALSE], paste0("seed_", s), co))
  }
}
write.table(res, file.path(d, "genefu_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("genefu calls written:", nrow(res), "\n")
