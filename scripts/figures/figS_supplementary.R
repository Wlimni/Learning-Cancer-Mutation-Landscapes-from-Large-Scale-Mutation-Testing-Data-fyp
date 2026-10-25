# Script Name: figS_supplementary.R
# Purpose : Supplementary figures.
#           S1 Phase 2 filter funnel; S2 exact driver combinations in NSCLC (UpSet, Lex 2014);
#           S3 one pair across cancer types, raw vs shrunk (forest plot, maftools Fig 4A style);
#           S4 rare combinations: observed vs expected; S5 dN/dS route vs Phase 2 effect;
#           S6 exclusive gene sets: Dendrix vs load-aware pathway coherence (CoMEt Fig 3 style bars).
# Inputs  : results/figures/data/*.csv
# Outputs : results/figures/figS1..S6_*.{pdf,png}
source("scripts/figures/_theme.R")
suppressPackageStartupMessages(library(ComplexHeatmap))

# S1 -- funnel
fn <- rd("fig2f_funnel.csv") %>% mutate(step = factor(step, rev(step))) %>%
  pivot_longer(c(co_occurring, exclusive), names_to = "direction") %>%
  mutate(direction = sub("_", "-", direction))
s1 <- ggplot(fn, aes(value, step, fill = direction)) +
  geom_col(width = 0.65) +
  geom_text(data = rd("fig2f_funnel.csv") %>% mutate(step = factor(step, rev(step))),
            aes(x = pairs, y = step, label = sprintf("%s (%.1f%%)", comma(pairs), pct_of_tested)),
            inherit.aes = FALSE, hjust = -0.1, size = 2.4) +
  scale_fill_manual(values = COL_DIR[1:2]) +
  scale_x_continuous(labels = comma, expand = expansion(c(0, 0.25))) +
  labs(x = "Gene pairs (summed over cancer types)", y = NULL, fill = NULL, title = "Phase 2: from tested pairs to replicated candidates") +
  theme(legend.position = "bottom")
save_fig(s1, "figS1_funnel", 120, 60)

# S2 -- UpSet of exact combinations among NSCLC drivers
cm <- rd("fig3b_lung_combinations.csv")
mat <- as.matrix(cm[, -1]) > 0
colnames(mat) <- sub("_fusion", " fusion", colnames(mat))
m <- make_comb_mat(mat, mode = "distinct")
m <- m[comb_size(m) >= 40]
us <- UpSet(m, comb_order = order(-comb_size(m)), set_order = order(-set_size(m)),
            top_annotation = upset_top_annotation(m, add_numbers = TRUE, numbers_gp = gpar(fontsize = 5),
                                                  annotation_name_gp = gpar(fontsize = 7)),
            right_annotation = upset_right_annotation(m, add_numbers = TRUE, numbers_gp = gpar(fontsize = 5),
                                                      annotation_name_gp = gpar(fontsize = 7)),
            row_names_gp = gpar(fontsize = 7, fontface = "italic"),
            column_title = "NSCLC: exact combinations of 8 drivers (combinations with >= 40 patients)",
            column_title_gp = gpar(fontsize = 8, fontface = "bold"))
s2 <- wrap_elements(full = grid.grabExpr(draw(us)))
save_fig(s2, "figS2_upset_lung", 183, 90)

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
save_fig(s3, "figS3_forest_kras_tp53", 110, 110)

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
       title = sprintf("%d replicated exclusive pairs: how rare is the combination?", nrow(rc))) +
  theme(legend.position = c(0.2, 0.8))
save_fig(s4, "figS4_rare_combinations", 110, 95)

# S5 -- dN/dS route vs Phase 2
dn <- rd("fig5c_dnds.csv") %>% mutate(sig = p < 0.05)
s5 <- ggplot(dn, aes(log2_fold_enrichment, log2(selection_ratio))) +
  geom_hline(yintercept = 0, linewidth = 0.2) + geom_vline(xintercept = 0, linewidth = 0.2) +
  geom_point(aes(colour = sig), size = 1.2, alpha = 0.8, stroke = 0) +
  geom_text_repel(data = dn %>% slice_min(p, n = 8), aes(label = paste0(A, " -> ", B, " (", short_ct(Cancer_Type), ")")),
                  size = 1.9, min.segment.length = 0, seed = 6) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#0072B2"), labels = c("p >= 0.05", "p < 0.05"), name = NULL) +
  labs(x = "Phase 2 effect, log2(observed / expected)", y = "log2 selection ratio of B\n(protein-changing : silent, A-mutant vs A-wild-type)",
       title = "An independent route agrees with the pair test",
       caption = sprintf("%d testable directions: %.0f%% same sign; %d of %d significant agree.",
                         nrow(dn), 100 * mean(sign(dn$log2_fold_enrichment) == sign(log2(dn$selection_ratio))),
                         sum(dn$sig & sign(dn$log2_fold_enrichment) == sign(log2(dn$selection_ratio))), sum(dn$sig))) +
  theme(legend.position = c(0.15, 0.9))
save_fig(s5, "figS5_dnds", 110, 95)

# S6 -- gene sets
gs <- rd("fig5d_gene_sets.csv") %>% pivot_longer(c(all3, any2), names_to = "metric") %>%
  mutate(metric = recode(metric, all3 = "All 3 genes in one pathway", any2 = ">= 2 genes in one pathway"),
         ct = short_ct(Cancer_Type))
s6 <- ggplot(gs, aes(method, value, fill = method)) +
  stat_summary(fun = mean, geom = "col", width = 0.6) +
  geom_line(aes(group = ct), colour = "grey60", linewidth = 0.3) +
  geom_point(size = 0.9, colour = "grey20") +
  facet_wrap(~metric) +
  scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  scale_x_discrete(labels = function(x) sub(" \\(", "\n(", x)) +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.05))) +
  labs(x = NULL, y = "Top-20 sets per cancer type", title = "Exclusive gene sets: which score recovers pathways?",
       caption = "Bars: mean over 6 cancer types; lines: one cancer type.")
save_fig(s6, "figS6_gene_sets", 110, 80)
