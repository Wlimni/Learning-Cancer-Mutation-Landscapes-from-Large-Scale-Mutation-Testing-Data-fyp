# Script Name: fig5_fusions_alleles.R
# Purpose : Figure 5 -- fusions and variant-level partners.
#           a lung cancer oncoprint of drivers incl. ALK/ROS1/RET fusions (Skoulidis 2015 Fig 2 style);
#           b fusion partners: significant fusion pairs per cancer type (lollipop of effect size);
#           c KRAS/EGFR variant partners before vs after adjusting for the smoking signature.
# Inputs  : results/figures/data/fig2e_lung_oncoprint.csv, fig2e_lung_n.csv, fig2e_fusion_pairs.csv, fig4d_smoking.csv
# Outputs : results/figures/fig5_fusions_alleles.{pdf,png}
source("scripts/figures/_theme.R")
suppressPackageStartupMessages(library(ComplexHeatmap))

# a -- oncoprint (patients with >= 1 alteration in these genes; cBioPortal colours)
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
pa <- wrap_elements(full = grid.grabExpr(draw(op, heatmap_legend_side = "bottom")))

# b -- fusion partners
fu <- rd("fig2e_fusion_pairs.csv") %>% filter(q_value < 0.05) %>%
  mutate(fusion = ifelse(grepl("_fusion", Gene_A), Gene_A, Gene_B), partner = ifelse(fusion == Gene_A, Gene_B, Gene_A),
         fusion = sub("_fusion", "", fusion), ct = short_ct(Cancer_Type))
top <- fu %>% filter(ct %in% c("NSCLC", "Prostate", "Hepatobiliary", "Breast")) %>%
  group_by(ct) %>% slice_min(q_value, n = 8, with_ties = FALSE) %>% ungroup() %>%
  mutate(lab = paste0(fusion, " fusion + ", sub("_fusion", " fusion", partner)),
         key = factor(paste(ct, lab, sep = "|"), levels = paste(ct, lab, sep = "|")[order(ct, log2_fold_enrichment)]))
pb <- ggplot(top, aes(log2_fold_enrichment, key, colour = direction)) +
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

# c -- variant partners: crude vs smoking-adjusted fold
sm <- rd("fig4d_smoking.csv") %>%
  mutate(lab = paste0(sub("_", " ", allele), " + ", partner), lab = reorder(lab, fold),
         lost = q_smoking_adjusted >= 0.05)
pc <- ggplot(sm, aes(y = lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_segment(aes(x = fold, xend = fold_smoking_adjusted, yend = lab), colour = "grey60", linewidth = 0.4) +
  geom_point(aes(x = fold), shape = 21, fill = "white", size = 1.4, stroke = 0.4) +
  geom_point(aes(x = fold_smoking_adjusted, colour = lost), size = 1.6) +
  scale_colour_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
                      labels = c("Still significant", "Not significant after adjustment")) +
  scale_x_log10() +
  labs(x = "Fold enrichment (log scale)", y = NULL, colour = NULL, title = "Variant partners and smoking (NSCLC)",
       caption = "Open: crude; filled: adjusted for the tobacco mutational signature.") +
  theme(axis.text.y = element_text(size = 6), legend.position = "bottom", legend.direction = "vertical")

fig <- pa / (pb | pc) + plot_layout(heights = c(0.9, 1.3)) + plot_annotation(tag_levels = "a")
save_fig(fig, "fig5_fusions_alleles", 183, 190)
