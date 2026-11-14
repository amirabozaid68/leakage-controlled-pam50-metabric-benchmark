
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("genefu", quietly = TRUE)) BiocManager::install("genefu", update = FALSE, ask = FALSE)
suppressPackageStartupMessages(library(genefu))
suppressWarnings(try(data(pam50.robust, package = "genefu"), silent = TRUE))
if (!exists("pam50.robust")) data(pam50, package = "genefu")
stopifnot(exists("pam50.robust"))
cent <- pam50.robust$centroids
write.table(data.frame(GENE = toupper(rownames(cent)), cent, check.names = FALSE),
            file.path(Sys.getenv("OUT"), "genefu", "pam50_centroids.tsv"), sep = "	", quote = FALSE, row.names = FALSE)
cat("std =", pam50.robust$std, "| method.cor =", pam50.robust$method.cor, "| genefu", as.character(packageVersion("genefu")), "
")
