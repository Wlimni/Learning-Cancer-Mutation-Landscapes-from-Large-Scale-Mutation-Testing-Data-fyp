# Script Name: phase2_pair_test.R
# Purpose : Phase 2 figures -- testing gene pairs fairly, and what it finds (notebook 02).
# Inputs  : results/figures/data/fig2*.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase2_pair_test/<step>_<name>.{pdf,png} -- one figure per step, shown in the phase's notebook
source("scripts/figures/_theme.R")
suppressPackageStartupMessages(library(ComplexHeatmap))
F <- "phase2_pair_test"

# 2.1 -- why not Fisher: p-values on data with no real gene-gene links (flat = honest)
M2 <- c(p_fisher = "Fisher's exact test (usual)", p_cond_subtype = "Conditional test (ours)")
nl <- rd("fig2a_null_pvalues.csv") %>% select(all_of(names(M2))) %>%
  pivot_longer(everything(), names_to = "m", values_to = "p") %>% mutate(method = factor(M2[m], M2))
lab <- nl %>% group_by(method) %>% summarise(txt = sprintf("p < 0.05: %.1f%% (should be 5%%)", 100 * mean(p < 0.05)))
p21 <- ggplot(nl, aes(p, fill = method)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.05), colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_text(data = lab, aes(x = 0.5, y = Inf, label = txt), vjust = 1.5, size = 2.4, inherit.aes = FALSE) +
  facet_wrap(~method, nrow = 1) +
  scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "p-value on data with no real gene-gene links", y = "Density",
       title = "Why not Fisher's test: it calls chance patterns significant",
       caption = "Same patients, links removed by shuffling. Flat = honest; spike at 0 = false positives.\nFull comparison of five methods: validation, Figure 5.1.")
save_fig(p21, "2.1_why_not_fisher", 130, 65, F)

# 2.2 -- filter funnel: from tested pairs to interaction candidates
fn0 <- rd("fig2f_funnel.csv") %>% filter(!grepl("replicated", step)) %>% mutate(step = factor(step, rev(step)))
fn <- fn0 %>% pivot_longer(c(co_occurring, exclusive), names_to = "direction") %>% mutate(direction = sub("_", "-", direction))
p22 <- ggplot(fn, aes(value, step, fill = direction)) +
  geom_col(width = 0.65) +
  geom_text(data = fn0, aes(x = pairs, y = step, label = sprintf("%s (%.1f%%)", comma(pairs), pct_of_tested)),
            inherit.aes = FALSE, hjust = -0.1, size = 2.4) +
  scale_fill_manual(values = COL_DIR[1:2]) +
  scale_x_continuous(labels = comma, expand = expansion(c(0, 0.25))) +
  labs(x = "Gene pairs (summed over cancer types)", y = NULL, fill = NULL, title = "From tested pairs to interaction candidates") +
  theme(legend.position = "bottom")
save_fig(p22, "2.2_filter_funnel", 120, 55, F)

pr <- rd("fig2_pairs.csv") %>%
  mutate(across(c(interaction_candidate, replicated), ~ .x == "True"),
         l2 = pmax(pmin(log2_fold_enrichment, 2), -2),
         cls = case_when(interaction_candidate ~ direction, q_value < 0.05 ~ "shared copy-number event",
                         TRUE ~ "not significant"))

# 2.3 -- pair heatmaps for three common cancers
tri <- function(ct, n_genes = 16) {
  d <- pr %>% filter(Cancer_Type == ct)
  g <- d %>% filter(interaction_candidate) %>% select(Gene_A, Gene_B) %>% unlist() %>% table() %>%
    sort(decreasing = TRUE) %>% names() %>% head(n_genes)
  d <- d %>% filter(Gene_A %in% g, Gene_B %in% g) %>%
    mutate(i = match(Gene_A, g), j = match(Gene_B, g), x = pmin(i, j), y = pmax(i, j))
  ggplot(d, aes(x, y, fill = l2)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, interaction_candidate), size = 0.7) +
    scale_fill_gradientn(colours = PRGN, limits = c(-2, 2), breaks = -2:2,
                         name = "log2(observed / expected)\npurple = exclusive, green = co-occurring") +
    guides(fill = guide_colourbar(barwidth = unit(45, "mm"), barheight = unit(2.5, "mm"), title.position = "top")) +
    scale_x_continuous(breaks = seq_along(g), labels = g, expand = c(0, 0), position = "top") +
    scale_y_reverse(breaks = seq_along(g), labels = g, expand = c(0, 0)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = short_ct(ct)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 5.5, face = "italic"),
          axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5))
}
p23 <- (tri("Non-Small Cell Lung Cancer") +
         labs(caption = "Dot: interaction candidate (q < 0.05, not one shared copy-number event).") |
       tri("Breast Cancer") | tri("Colorectal Cancer")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# 2.4 -- all tested pairs (volcano) (pseudo-log y axis so the very strong pairs do not pile up at a cap)
vb <- pr %>% mutate(y = -log10(pmax(q_value, 1e-300)), lab = paste0(Gene_A, "/", Gene_B, " (", short_ct(Cancer_Type), ")"),
                    cls = factor(cls, names(COL_DIR)))
top <- vb %>% filter(interaction_candidate, q_value < 1e-10) %>% group_by(direction) %>%
  slice_max(abs(log2_fold_enrichment), n = 5, with_ties = FALSE)
cnt <- vb %>% count(cls)
p24 <- ggplot(vb %>% arrange(desc(cls)), aes(log2_fold_enrichment, y, colour = cls)) +
  geom_point(size = 0.5, alpha = 0.7, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  geom_text_repel(data = top, aes(label = lab), size = 2, colour = "black", max.overlaps = 30,
                  min.segment.length = 0, segment.size = 0.2, seed = 2, force = 2, max.iter = 20000) +
  scale_colour_manual(values = COL_DIR, labels = function(x) paste0(x, " (n = ", comma(cnt$n[match(x, cnt$cls)]), ")")) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), breaks = c(0, 1, 3, 10, 30, 100, 300),
                     limits = c(0, 320), expand = expansion(c(0, 0.02))) +
  scale_x_continuous(limits = c(-7.2, 4.5)) +
  labs(x = "log2(observed / expected co-altered patients)", y = "-log10 q", colour = NULL,
       title = sprintf("%s pairs in %d cancer types", comma(nrow(vb)), n_distinct(vb$Cancer_Type))) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = c(0.27, 0.86))
save_fig(p23, "2.3_pair_heatmaps", 183, 85, F)
save_fig(p24, "2.4_all_pairs_volcano", 120, 100, F)

# 2.5 -- fusion partners
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
save_fig(p25, "2.5_fusion_partners", 120, 125, F)

# 2.6 -- lung oncoprint (patients with >= 1 alteration in these genes; cBioPortal colours)
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
save_fig(p26, "2.6_lung_oncoprint", 183, 80, F)
