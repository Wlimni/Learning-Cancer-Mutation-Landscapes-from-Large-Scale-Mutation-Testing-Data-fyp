#!/bin/bash
# Script Name: run_all.sh
# Purpose : Rebuild every figure from the processed data (run from the repo root).
# Outputs : results/figures/fig1-6_*.{pdf,png}, figS1-S6_*.{pdf,png}
# Note    : needs R with ggplot2, patchwork, ggrepel, igraph and ComplexHeatmap (Bioconductor).
set -e
.venv/bin/python scripts/figures/00_export_figure_data.py
for f in scripts/figures/fig*.R; do Rscript "$f"; done
