#!/bin/bash
# Script Name: run_all.sh
# Purpose : Rebuild every figure from the processed data (run from the repo root).
# Outputs : results/figures/<phase folder>/figure_<phase>.{pdf,png}, one combined figure per phase
# Note    : needs R with ggplot2, patchwork, ggrepel and ComplexHeatmap (Bioconductor).
set -e
.venv/bin/python scripts/figures/00_export_figure_data.py
for f in phase1_data phase2_pair_test phase2b_validation phase2c_extensions phase2c_methods phase3_lookup_tool phase4_biology; do Rscript "scripts/figures/$f.R"; done
