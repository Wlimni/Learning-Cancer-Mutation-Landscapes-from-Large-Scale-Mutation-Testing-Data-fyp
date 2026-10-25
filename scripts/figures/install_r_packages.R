# Script Name: install_r_packages.R
# Purpose : Install the R packages the figure scripts (and scripts/validation/rediscover_run.R) need.
# Usage   : Rscript scripts/figures/install_r_packages.R
# Notes   : tested with R 4.5.1. ComplexHeatmap comes from Bioconductor. No package here needs a C++
#           compiler, which matters on machines where source builds fail (e.g. arrow, ggraph).
options(repos = "https://cloud.r-project.org")
cran <- c("ggplot2", "dplyr", "tidyr", "patchwork", "scales", "ggrepel", "igraph", "Rediscover")
missing <- setdiff(cran, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
for (p in c(cran, "ComplexHeatmap")) cat(sprintf("%-15s %s\n", p, if (requireNamespace(p, quietly = TRUE)) "ok" else "MISSING"))
