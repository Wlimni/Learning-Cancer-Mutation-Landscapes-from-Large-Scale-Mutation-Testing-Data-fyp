#!/bin/bash
# Script Name: run_all.sh
# Purpose : Rebuild every figure from the processed data (run from the repo root).
# Outputs : results/figures/<phase folder>/<step>_<name>.{pdf,png}, one figure per step
# Note    : needs R with ggplot2, patchwork, ggrepel and ComplexHeatmap (Bioconductor).
set -e
.venv/bin/python scripts/figures/00_export_figure_data.py
for f in phase1_data phase2_pair_test phase3_lookup_tool phase4_biology validation; do Rscript "scripts/figures/$f.R"; done
