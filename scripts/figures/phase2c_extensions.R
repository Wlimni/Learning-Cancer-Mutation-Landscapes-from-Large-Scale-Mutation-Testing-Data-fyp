# Script Name: phase2c_extensions.R
# Purpose : Phase 2c figures -- side analyses of the pair test (notebook 02c): the fusion screen.
# Inputs  : results/figures/data/fig2e*.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase2c_extensions/figure_phase2c.{pdf,png} -- one combined figure (2 panels)
source("scripts/figures/_theme.R")
suppressPackageStartupMessages(library(ComplexHeatmap))
OUT_DIR <- "phase2c_extensions"

# a -- fusion partners
fu <- rd("fig2e_fusion_pairs.csv") %>% filter(q_value < 0.05) %>%
  mutate(fusion = ifelse(grepl("_fusion", Gene_A), Gene_A, Gene_B), partner = ifelse(fusion == Gene_A, Gene_B, Gene_A),
         fusion = sub("_fusion", "", fusion), ct = short_ct(Cancer_Type))
top <- fu %>% filter(ct %in% c("NSCLC", "Prostate", "Hepatobiliary", "Breast")) %>%
  group_by(ct) %>% slice_min(q_value, n = 8, with_ties = FALSE) %>% ungroup() %>%
  mutate(lab = paste0(fusion, " fusion + ", sub("_fusion", " fusion", partner)),
         key = factor(paste(ct, lab, sep = "|"), levels = paste(ct, lab, sep = "|")[order(ct, log2_fold_enrichment)]))
p25 <- ggplot(top, aes(log2_fold_enrichment, key, colour = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2_fold_enrichment, yend = key), linewidth = 0.4) +
  geom_point(aes(size = n_both_altered)) +
  facet_grid(ct ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*\\|", "", x)) +
  scale_colour_manual(values = COL_DIR[1:2]) +
  scale_size(range = c(0.6, 2.5), name = "Patients\nwith both") +
  labs(x = "log2(observed / expected)", y = NULL, colour = NULL, title = "Fusion partners (q < 0.05)",
       caption = sprintf("%d of %d fusion pairs significant; 8 strongest per cancer type shown.\nBreast FGFR1 fusion + FGFR1-region genes: same-locus artefact (flagged in Phase 2).",
                         nrow(fu), nrow(rd("fig2e_fusion_pairs.csv")))) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), axis.text.y = element_text(size = 6), legend.position = "right")

# b -- lung oncoprint (patients with >= 1 alteration in these genes; cBioPortal colours)
on <- rd("fig2e_lung_oncoprint.csv")
n_lung <- rd("fig2e_lung_n.csv")$n_patients
GENES <- c("EGFR", "KRAS", "ALK_fusion", "ROS1_fusion", "RET_fusion", "BRAF", "MET", "ERBB2", "STK11", "KEAP1", "TP53")
m <- on %>% pivot_wider(names_from = Sample_ID, values_from = type, values_fill = "") %>%
  tibble::column_to_rownames("Hugo_Symbol") %>% as.matrix()
m <- m[intersect(GENES, rownames(m)), ]
rownames(m) <- sub("_fusion", " fusion", rownames(m))
alter_fun <- list(
  background = function(x, y, w, h) grid.rect(x, y, w, h * 0.9, gp = gpar(fill = "grey93", col = NA)),
  MUT = function(x, y, w, h) grid.rect(x, y, w, h * 0.9, gp = gpar(fill = COL_ALT[["MUT"]], col = NA)),
  AMP = function(x, y, w, h) grid.rect(x, y, w, h * 0.9, gp = gpar(fill = COL_ALT[["AMP"]], col = NA)),
  DEL = function(x, y, w, h) grid.rect(x, y, w, h * 0.9, gp = gpar(fill = COL_ALT[["DEL"]], col = NA)),
  FUSION = function(x, y, w, h) grid.rect(x, y, w, h * 0.45, gp = gpar(fill = COL_ALT[["FUSION"]], col = NA)))
# column order: driver groups first (EGFR, KRAS, fusions, ...), so exclusivity is visible as a staircase
ord <- do.call(order, c(lapply(rownames(m), function(g) -(m[g, ] != "")), list(decreasing = FALSE)))
op <- oncoPrint(m, alter_fun = alter_fun, col = COL_ALT, column_order = ord, row_order = seq_len(nrow(m)),
                remove_empty_columns = TRUE, show_column_names = FALSE, top_annotation = NULL,
                right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(border = FALSE, width = unit(10, "mm"),
                                                                               axis_param = list(gp = gpar(fontsize = 5)))),
                row_names_gp = gpar(fontsize = 7, fontface = "italic"), pct_gp = gpar(fontsize = 6),
                heatmap_legend_param = list(title = "Alteration", at = names(COL_ALT),
                                            labels = c("Mutation", "Amplification", "Deep deletion", "Fusion"), title_gp = gpar(fontsize = 7, fontface = "bold"),
                                            labels_gp = gpar(fontsize = 6), nrow = 1),
                column_title = sprintf("NSCLC: %s of %s patients tested for mutations, copy number and fusions",
                                       format(ncol(m), big.mark = ","), format(n_lung, big.mark = ",")),
                column_title_gp = gpar(fontsize = 8, fontface = "bold"))
p26 <- wrap_elements(full = grid.grabExpr(draw(op, heatmap_legend_side = "bottom")))

# One combined figure for the phase; the title states the main message
fig <- p26 / p25 + plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = "Phase 2c: fusions avoid the other drivers -- most lung cancers carry only one of them")
save_fig(fig, "figure_phase2c", 183, 200, OUT_DIR)
