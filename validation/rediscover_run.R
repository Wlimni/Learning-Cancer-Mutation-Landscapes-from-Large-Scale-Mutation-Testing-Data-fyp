# Script Name: rediscover_run.R
# Purpose : Run Rediscover (Ferrer-Bonsoms et al., Bioinformatics 2022) on the
#           MSK-IMPACT468 benchmark: its rate fit (getPM) and its exact
#           co-occurrence test (getMutex), for rediscover_compare.py.
# Inputs  : data/validation/msk468_genes_by_samples.csv (genes x samples, 0/1)
# Outputs : data/validation/rediscover_pm.csv  (genes x samples)
#           data/validation/rediscover_p.csv   (i, j, one-sided upper-tail p)
# Notes   : install.packages("Rediscover") needs Bioconductor's maftools and a
#           working C++ toolchain (PoissonBinomial links against fftw).
#           Run from the repo root: Rscript validation/rediscover_run.R [null]

suppressMessages(library(Rediscover))
out <- "data/validation"

# optional argument "null": run on the no-interaction matrix from null_calibration.py
args <- commandArgs(trailingOnly = TRUE)
null_run <- length(args) > 0 && args[1] == "null"
A <- as.matrix(read.csv(file.path(out, if (null_run) "null_genes_by_samples.csv" else "msk468_genes_by_samples.csv"), header = FALSE))
storage.mode(A) <- "integer"
cat(sprintf("A (genes x samples): %d x %d, altered cells %d\n", nrow(A), ncol(A), sum(A)))

PM <- getPM(A)                                   # S4 PMatrix, kept as-is for getMutex
if (null_run) {
  for (lt in c(TRUE, FALSE)) {
    P <- as.matrix(getMutex(A, PM = PM, lower.tail = lt, method = "Exact", mixed = FALSE)); p <- P[upper.tri(P)]
    cat(sprintf("%s: frac p<0.05 = %.3f | BH q<0.05: %d of %d\n", ifelse(lt, "exclusivity", "co-occurrence"),
                mean(p < 0.05), sum(p.adjust(p, "BH") < 0.05), length(p)))
  }
  quit(save = "no")
}
write.csv(as.matrix(PM), file.path(out, "rediscover_pm.csv"), row.names = FALSE)

# co-occurrence = upper tail; exact method, no mixed approximation
P <- as.matrix(getMutex(A, PM = PM, lower.tail = FALSE, method = "Exact", mixed = FALSE))
idx <- which(upper.tri(P), arr.ind = TRUE)
write.csv(data.frame(i = idx[, 1] - 1, j = idx[, 2] - 1, p_rediscover = P[idx]),
          file.path(out, "rediscover_p.csv"), row.names = FALSE)
cat(sprintf("wrote %d pairs\n", nrow(idx)))
