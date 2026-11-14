
suppressPackageStartupMessages(library(genefu))
suppressWarnings(try(data(pam50.robust, package = "genefu"), silent = TRUE))
if (!exists("pam50.robust")) data(pam50, package = "genefu")
stopifnot(exists("pam50.robust"))
d   <- Sys.getenv("GDIR")
X   <- read.delim(file.path(d, "pam50_expr_all_samples.tsv"), check.names = FALSE)
rownames(X) <- X$SAMPLE_ID; X$SAMPLE_ID <- NULL; X <- as.matrix(X)
stopifnot(all(rownames(pam50.robust$centroids) %in% colnames(X)))
annot <- data.frame(probe = colnames(X), EntrezGene.ID = NA, row.names = colnames(X))
run <- function(M, tag) {
  r <- intrinsic.cluster.predict(sbt.model = pam50.robust, data = M, annot = annot,
                                 do.mapping = FALSE, do.prediction.strength = FALSE, verbose = FALSE)
  out <- data.frame(SAMPLE_ID = rownames(M), run = tag, call = as.character(r$subtype),
                    r$subtype.proba, check.names = FALSE)
  out
}
tests <- read.delim(file.path(d, "test_sets.tsv"))
res <- run(X, "full_cohort")
for (s in unique(tests$seed)) {
  ids <- as.character(tests$SAMPLE_ID[tests$seed == s])
  res <- rbind(res, run(X[ids, , drop = FALSE], paste0("seed_", s)))
}
write.table(res, file.path(d, "genefu_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("genefu calls written:", nrow(res), "\n")
