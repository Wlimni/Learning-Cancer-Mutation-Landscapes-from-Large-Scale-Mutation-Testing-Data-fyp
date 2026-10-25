# Script Name: fig4_biology.R
# Purpose : Figure 4 -- what the patterns say about biology.
#           a exclusivity / co-occurrence within and between the 10 Sanchez-Vega pathways (Sanchez-Vega 2018 Fig 6);
#           b the same gene pair across cancer types: universal vs reversing (Park & Lehner 2015);
#           c copy-number burden confound: effect before vs after adjusting (van de Haar 2019 Fig 3 idea).
# Inputs  : results/figures/data/fig2_pairs.csv, fig4a_pathways.csv, fig4b_heterogeneity.csv, fig4e_burden.csv
# Outputs : results/figures/fig4_biology.{pdf,png}
source("scripts/figures/_theme.R")

pr <- rd("fig2_pairs.csv") %>% mutate(across(c(interaction_candidate, replicated), ~ .x == "True"))
pw <- rd("fig4a_pathways.csv")

# a1 -- pathway x pathway: share of tested pairs that are exclusive (lower) / co-occurring (upper) candidates
pp <- pr %>% inner_join(pw, by = c(Gene_A = "gene"), relationship = "many-to-many") %>% rename(pA = pathway) %>%
  inner_join(pw, by = c(Gene_B = "gene"), relationship = "many-to-many") %>% rename(pB = pathway) %>%
  mutate(p1 = pmin(pA, pB), p2 = pmax(pA, pB))
lv <- sort(unique(pw$pathway))
cell <- pp %>% group_by(p1, p2) %>%
  summarise(n = n(), ex = mean(interaction_candidate & direction == "exclusive"),
            co = mean(interaction_candidate & direction == "co-occurring"), .groups = "drop") %>%
  filter(n >= 10)
mat <- bind_rows(cell %>% transmute(x = p1, y = p2, v = -ex, lab = ex, n),     # lower triangle + diagonal: exclusive
                 cell %>% filter(p1 != p2) %>% transmute(x = p2, y = p1, v = co, lab = co, n)) %>%
  mutate(x = factor(x, lv), y = factor(y, rev(lv)))
pa1 <- ggplot(mat, aes(x, y, fill = v)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(100 * lab)), size = 1.9) +
  geom_tile(data = data.frame(x = factor(lv, lv), y = factor(lv, rev(lv))), aes(x, y), inherit.aes = FALSE,
            fill = NA, colour = "black", linewidth = 0.4) +
  scale_fill_gradientn(colours = PRGN, limits = c(-0.4, 0.4), oob = squish,
                       breaks = c(-0.4, 0, 0.4), labels = c("40% exclusive", "0", "40% co-occurring"), name = NULL) +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  labs(x = NULL, y = NULL, title = "Pathway structure",
       caption = "% of tested pairs that are candidates. Lower/diagonal (boxed): exclusive;\nupper: co-occurring. Cells with >= 10 tested pairs.") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x.top = element_text(angle = 45, hjust = 0),
        legend.position = "bottom", legend.key.width = unit(8, "mm"))

# a2 -- within vs between pathways
gp <- split(pw$pathway, pw$gene)
wb <- pr %>% filter(Gene_A %in% names(gp), Gene_B %in% names(gp)) %>%
  mutate(same = mapply(function(a, b) length(intersect(gp[[a]], gp[[b]])) > 0, Gene_A, Gene_B)) %>%
  group_by(same) %>%
  summarise(exclusive = mean(interaction_candidate & direction == "exclusive"),
            `co-occurring` = mean(interaction_candidate & direction == "co-occurring"), n = n()) %>%
  pivot_longer(c(exclusive, `co-occurring`), names_to = "direction") %>%
  mutate(where = ifelse(same, "Same pathway", "Different pathways"))
pa2 <- ggplot(wb, aes(where, value, fill = direction)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = percent(value, 1)), position = position_dodge(0.75), vjust = -0.3, size = 2.2) +
  scale_fill_manual(values = COL_DIR[1:2]) +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.15))) +
  labs(x = NULL, y = "Tested pairs that are candidates", fill = NULL, title = "Exclusivity is pathway-bound") +
  theme(legend.position = c(0.25, 0.9))

# b -- pair x cancer type heatmap: reversing vs universal pairs
het <- rd("fig4b_heterogeneity.csv")
sel <- bind_rows(het %>% filter(pattern == "opposite directions") %>% slice_min(q, n = 10) %>% mutate(grp = "Reverses direction"),
                 het %>% filter(pattern == "universal") %>% slice_max(n_cancer_types, n = 6, with_ties = FALSE) %>% mutate(grp = "Universal"))
hb <- pr %>% mutate(pair = paste(Gene_A, "/", Gene_B)) %>% filter(pair %in% sel$pair) %>%
  left_join(sel %>% select(pair, grp), by = "pair") %>%
  mutate(ct = short_ct(Cancer_Type), l2 = pmax(pmin(log2_fold_enrichment, 2), -2))
ct_keep <- hb %>% filter(interaction_candidate) %>% count(ct, sort = TRUE) %>% filter(n >= 2) %>% pull(ct)
hb <- hb %>% filter(ct %in% ct_keep) %>% mutate(pair = factor(pair, rev(sel$pair)))
pb <- ggplot(hb, aes(ct, pair, fill = l2)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = ~ filter(.x, interaction_candidate), size = 0.6) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(colours = PRGN, limits = c(-2, 2), name = "log2 obs/exp") +
  labs(x = NULL, y = NULL, title = "Same pair, different cancer types",
       caption = "Dot: interaction candidate in that cancer type. Blank: not tested there.") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic", size = 6), legend.key.height = unit(5, "mm"))

# c -- burden confound: dumbbell before / after adjusting for copy-number burden
bu <- rd("fig4e_burden.csv") %>%
  mutate(lab = paste0(Gene_A, "/", Gene_B, " (", short_ct(Cancer_Type), ")"),
         verdict = ifelse(grepl("^survives", verdict), "Survives adjustment", "Explained by burden"),
         lab = reorder(lab, burden_adjusted_OR))
pc <- ggplot(bu, aes(y = lab)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_segment(aes(x = phase2_fold, xend = burden_adjusted_OR, yend = lab), colour = "grey60", linewidth = 0.4) +
  geom_point(aes(x = phase2_fold), shape = 21, fill = "white", size = 1.4, stroke = 0.4) +
  geom_point(aes(x = burden_adjusted_OR, colour = verdict), size = 1.6) +
  scale_colour_manual(values = c("Survives adjustment" = "#1B7837", "Explained by burden" = "#D55E00")) +
  scale_x_log10() +
  labs(x = "Effect size (log scale)", y = NULL, colour = NULL,
       title = "Co-amplification and genome instability",
       caption = sprintf("Open: Phase 2 fold; filled: odds ratio adjusted for copy-number burden.\n%d copy-number-driven co-occurring candidates; %d explained by burden.",
                         nrow(bu), sum(bu$verdict == "Explained by burden"))) +
  theme(axis.text.y = element_text(size = 5.5), legend.position = "bottom")

fig <- ((pa1 | pa2) + plot_layout(widths = c(1.3, 1))) / (pb | pc) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a")
save_fig(fig, "fig4_biology", 183, 200)
