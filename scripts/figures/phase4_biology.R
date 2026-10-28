# Script Name: phase4_biology.R
# Purpose : Phase 4 figures -- what the patterns say about biology (notebook 04).
# Inputs  : results/figures/data/fig2_pairs.csv, fig4*.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase4_biology/<step>_<name>.{pdf,png} -- one figure per step, shown in the phase's notebook
source("scripts/figures/_theme.R")
F <- "phase4_biology"

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
save_fig((pa1 | pa2) + plot_layout(widths = c(1.3, 1)), "4.1_pathway_structure", 183, 95, F)

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
save_fig(pb, "4.2_pairs_across_cancer_types", 140, 110, F)

# S3 -- forest plot: KRAS/TP53 across cancer types, raw vs shrunk (log2 scale)
sh <- rd("fig4c_shrunk.csv") %>% filter(key == "KRAS/TP53") %>%
  mutate(raw = theta / log(2), se = sqrt(var) / log(2), shrunk = theta_shrunk / log(2),
         ct = reorder(short_ct(Cancer_Type), raw), sig = ifelse(q_value < 0.05, direction, "not significant"))
s3 <- ggplot(sh, aes(y = ct)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(xmin = raw - 1.96 * se, xmax = raw + 1.96 * se), width = 0, orientation = "y", linewidth = 0.3, colour = "grey50") +
  geom_point(aes(x = raw, colour = sig), size = 1.6) +
  geom_point(aes(x = shrunk), shape = 4, size = 1.4) +
  scale_colour_manual(values = COL_DIR[c("co-occurring", "exclusive", "not significant")]) +
  labs(x = "log2(observed / expected) with 95% CI", y = NULL, colour = NULL,
       title = "KRAS/TP53: opposite directions in different cancer types",
       caption = "Dot: per-cancer-type estimate; cross: shrunk towards the across-cancer mean (DerSimonian-Laird).") +
  theme(legend.position = "bottom")
save_fig(s3, "4.3_example_kras_tp53", 110, 110, F)

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
save_fig(pc, "4.4_variant_partners_and_smoking", 110, 105, F)

# S4 -- rare combinations
rc <- rd("fig4f_rare.csv") %>% mutate(civic = civic == "True")
s4 <- ggplot(rc, aes(expected_both, patients_with_both + 0.5)) +
  geom_abline(linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(colour = civic), size = 1, alpha = 0.8, stroke = 0) +
  geom_text_repel(data = rc %>% filter(civic) %>% slice_max(times_rarer_than_expected, n = 8),
                  aes(label = paste0(Gene_A, "/", Gene_B, " (", short_ct(Cancer_Type), ")")),
                  size = 2, min.segment.length = 0, seed = 5) +
  scale_x_log10() + scale_y_log10(breaks = c(0.5, 1.5, 5.5, 10.5, 50.5, 100.5), labels = c(0, 1, 5, 10, 50, 100)) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#D55E00"), labels = c("No", "Yes"),
                      name = "CIViC evidence\nfor the combination") +
  labs(x = "Expected patients with both", y = "Observed patients with both",
       caption = "Labelled: pairs with CIViC evidence for the combination (most depleted first).",
       title = sprintf("Rare combinations: %d exclusive pairs, observed vs expected", nrow(rc))) +
  theme(legend.position = c(0.2, 0.8))
save_fig(s4, "4.5_rare_combinations", 110, 95, F)

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
save_fig(pc, "4.6_copy_number_burden_check", 120, 110, F)
