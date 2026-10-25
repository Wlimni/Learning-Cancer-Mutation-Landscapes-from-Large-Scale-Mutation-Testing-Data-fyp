# Script Name: fig3_landscape.R
# Purpose : Figure 3 -- the co-occurrence / exclusivity landscape.
#           a triangular pair heatmaps per cancer type (maftools somaticInteractions; Sanchez-Vega 2018 Fig 6);
#           b volcano of all tested pairs (log fold vs significance; e.g. PMC12482361 Fig 7A);
#           c interaction network of replicated candidates in lung cancer (Park & Lehner 2015 Fig 3).
# Inputs  : results/figures/data/fig2_pairs.csv
# Outputs : results/figures/fig3_landscape.{pdf,png}
# Notes   : heatmap genes = the 16 genes in most interaction candidates in that cancer type.
source("scripts/figures/_theme.R")
suppressPackageStartupMessages(library(igraph))

pr <- rd("fig2_pairs.csv") %>%
  mutate(across(c(interaction_candidate, replicated), ~ .x == "True"),
         l2 = pmax(pmin(log2_fold_enrichment, 2), -2),
         cls = case_when(interaction_candidate ~ direction, q_value < 0.05 ~ "shared copy-number event",
                         TRUE ~ "not significant"))

# a -- triangle heatmaps
tri <- function(ct, n_genes = 16) {
  d <- pr %>% filter(Cancer_Type == ct)
  g <- d %>% filter(interaction_candidate) %>% select(Gene_A, Gene_B) %>% unlist() %>% table() %>%
    sort(decreasing = TRUE) %>% names() %>% head(n_genes)
  d <- d %>% filter(Gene_A %in% g, Gene_B %in% g) %>%
    mutate(i = match(Gene_A, g), j = match(Gene_B, g), x = pmin(i, j), y = pmax(i, j))
  ggplot(d, aes(x, y, fill = l2)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, interaction_candidate & replicated), size = 0.7) +
    geom_point(data = ~ filter(.x, interaction_candidate & !replicated), size = 0.7, shape = 1, stroke = 0.3) +
    scale_fill_gradientn(colours = PRGN, limits = c(-2, 2), breaks = -2:2,
                         name = "log2(observed / expected)\npurple = exclusive, green = co-occurring") +
    guides(fill = guide_colourbar(barwidth = unit(35, "mm"), barheight = unit(2.5, "mm"), title.position = "left")) +
    scale_x_continuous(breaks = seq_along(g), labels = g, expand = c(0, 0), position = "top") +
    scale_y_reverse(breaks = seq_along(g), labels = g, expand = c(0, 0)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = short_ct(ct)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 5.5, face = "italic"),
          axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5))
}
pa <- (tri("Non-Small Cell Lung Cancer") +
         labs(caption = "Open circle: interaction candidate (q < 0.05, not one shared copy-number event). Filled: also replicated outside MSK.") |
       tri("Breast Cancer") | tri("Colorectal Cancer")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# b -- volcano (pseudo-log y axis so the very strong pairs do not pile up at a cap)
vb <- pr %>% mutate(y = -log10(pmax(q_value, 1e-300)), lab = paste0(Gene_A, "/", Gene_B, " (", short_ct(Cancer_Type), ")"),
                    cls = factor(cls, names(COL_DIR)))
top <- vb %>% filter(interaction_candidate, replicated, q_value < 1e-10) %>% group_by(direction) %>%
  slice_max(abs(log2_fold_enrichment), n = 5, with_ties = FALSE)
cnt <- vb %>% count(cls)
pb <- ggplot(vb %>% arrange(desc(cls)), aes(log2_fold_enrichment, y, colour = cls)) +
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

# c -- network: replicated candidates in lung cancer
ed <- pr %>% filter(Cancer_Type == "Non-Small Cell Lung Cancer", interaction_candidate, replicated)
g <- graph_from_data_frame(ed %>% select(Gene_A, Gene_B), directed = FALSE)
set.seed(3); xy <- layout_with_kk(g)
nodes <- data.frame(gene = V(g)$name, x = xy[, 1], y = xy[, 2])
e2 <- ed %>% mutate(x = nodes$x[match(Gene_A, nodes$gene)], y = nodes$y[match(Gene_A, nodes$gene)],
                    xend = nodes$x[match(Gene_B, nodes$gene)], yend = nodes$y[match(Gene_B, nodes$gene)])
pc <- ggplot() +
  geom_segment(data = e2, aes(x, y, xend = xend, yend = yend, colour = direction, linewidth = abs(log2_fold_enrichment)),
               alpha = 0.6) +
  geom_label(data = nodes, aes(x, y, label = gene), size = 2, fontface = "italic", linewidth = 0.15,
             label.padding = unit(0.6, "mm"), fill = "white") +
  scale_colour_manual(values = COL_DIR[1:2]) +
  scale_linewidth(range = c(0.15, 1.2), name = "Effect size\n(abs. log2)") +
  scale_x_continuous(expand = expansion(0.12)) + scale_y_continuous(expand = expansion(0.05)) +
  labs(colour = NULL, title = sprintf("NSCLC: %d replicated interactions", nrow(ed))) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9), plot.title.position = "plot", legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 12), legend.key.size = unit(3, "mm"))

fig <- pa / (pb | pc) + plot_layout(heights = c(0.85, 1)) + plot_annotation(tag_levels = list(c("a", "", "", "b", "c")))
save_fig(fig, "fig3_landscape", 183, 190)
