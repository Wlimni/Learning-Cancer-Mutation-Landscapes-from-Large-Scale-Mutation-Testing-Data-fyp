# Script Name: _theme.R
# Purpose : Shared style for all figure scripts -- one theme, one colour convention, one saver.
# Notes   : colours follow the co-mutation literature (Sanchez-Vega 2018, maftools, SelectSim):
#           green = co-occurring, purple = mutually exclusive (ColorBrewer PRGn, colour-blind safe);
#           methods use the Okabe-Ito palette; oncoprints use cBioPortal's alteration colours.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales); library(ggrepel)
})

D <- "results/figures/data/"
FIG <- "results/figures/"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(paste0(D, f), check.names = FALSE, stringsAsFactors = FALSE)

COL_DIR <- c("co-occurring" = "#1B7837", "exclusive" = "#762A83", "shared copy-number event" = "grey45",
             "not significant" = "grey80")
PRGN <- c("#762A83", "#C2A5CF", "#F7F7F7", "#A6DBA0", "#1B7837")   # purple - white - green
COL_METHOD <- c("Fisher" = "#E69F00", "DISCOVER-style" = "#56B4E9", "Conditional" = "#009E73",
                "Conditional + subtype" = "#0072B2", "SELECT" = "#CC79A7")
COL_ALT <- c(MUT = "#008000", AMP = "#E41A1C", DEL = "#2166AC", FUSION = "#8B00FF")

theme_fyp <- function(base = 8) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          axis.text = element_text(colour = "black"), strip.background = element_blank(),
          strip.text = element_text(face = "bold"), legend.key.size = unit(3, "mm"),
          plot.title = element_text(face = "bold", size = base + 1), plot.tag = element_text(face = "bold", size = base + 4),
          panel.grid = element_blank(), plot.title.position = "plot", plot.caption.position = "plot",
          plot.caption = element_text(hjust = 0, size = base - 1.5, colour = "grey30"))
}
theme_set(theme_fyp())

short_ct <- function(x) {
  x <- sub("Non-Small Cell Lung Cancer", "NSCLC", x); x <- sub(" Cancer$", "", x)
  x <- sub("Myeloproliferative Neoplasms", "MPN", x); x <- sub("Esophagogastric", "Oesophagogastric", x); x
}

# Save as vector PDF (for the thesis/paper) and 300-dpi PNG (for slides and notebooks), one folder per phase
save_fig <- function(p, name, w, h, folder = "") {
  out <- paste0(FIG, folder, if (nzchar(folder)) "/" else "")
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0(out, name, ".pdf"), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(paste0(out, name, ".png"), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  message("saved ", folder, "/", name)
}
