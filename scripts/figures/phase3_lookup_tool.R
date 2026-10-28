# Script Name: phase3_lookup_tool.R
# Purpose : Phase 3 figures -- what the lookup tool answers (notebook 03).
# Inputs  : results/figures/data/fig3*.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase3_lookup_tool/figure_phase3.{pdf,png} -- one combined figure (4 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase3_lookup_tool"

# a -- KRAS G12C partners in lung cancer
kr <- rd("fig3a_kras_g12c.csv") %>%
  mutate(sig = ifelse(q_value < 0.05, direction, "not significant"),
         partner = reorder(partner, log2_fold_enrichment),
         lab = sprintf("%d / %.0f", n_both_altered, expected_both))
pa <- ggplot(kr, aes(log2_fold_enrichment, partner, colour = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2_fold_enrichment, yend = partner), linewidth = 0.4) +
  geom_point(size = 1.6) +
  geom_text(aes(x = 2.1, label = lab), colour = "grey30", size = 1.9, hjust = 1) +
  annotate("text", x = 2.1, y = nrow(kr) + 1, label = "observed / expected", size = 1.9, hjust = 1, fontface = "bold") +
  scale_colour_manual(values = COL_DIR[c("co-occurring", "exclusive", "not significant")]) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.6))) +
  coord_cartesian(xlim = c(-3.8, 2.1)) +
  labs(x = "log2(observed / expected)", y = NULL, colour = NULL, title = "lookup(): partners of KRAS G12C in NSCLC") +
  theme(axis.text.y = element_text(face = "italic", size = 6), legend.position = "bottom")

# 3.2 -- the 10 most common exact combinations of eight lung-cancer drivers (what combo_prevalence() counts)
cm <- rd("fig3b_lung_combinations.csv")
genes <- setdiff(names(cm), names(cm)[1])
combo <- apply(cm[, genes] > 0, 1, function(r) if (any(r)) paste(sub("_fusion", " fusion", genes[r]), collapse = " + ") else NA)
tc <- as.data.frame(table(combo), stringsAsFactors = FALSE) %>% slice_max(Freq, n = 10, with_ties = FALSE) %>%
  mutate(combo = factor(combo, rev(combo)))
s2 <- ggplot(tc, aes(Freq, combo)) +
  geom_col(fill = "grey45", width = 0.65) +
  geom_text(aes(label = comma(Freq)), hjust = -0.15, size = 2.5) +
  scale_x_continuous(expand = expansion(c(0, 0.15))) +
  labs(x = "Lung cancer patients", y = NULL, title = "combo_prevalence(): the most common exact combinations",
       caption = sprintf("Among %s lung cancers with at least one of 8 drivers; each patient counted once, by their exact set.",
                         comma(sum(!is.na(combo))))) +
  theme(axis.text.y = element_text(face = "italic"))

# b -- partner prediction AUROC
pp <- rd("fig3c_predict_partners.csv") %>%
  mutate(method = recode(method, baseline = "Gene base rates", subtype = "Base rates by subtype",
                         logistic = "Logistic model (whole profile)"),
         method = factor(method, c("Gene base rates", "Base rates by subtype", "Logistic model (whole profile)")),
         ct = reorder(short_ct(cancer_type), auroc, max))
pb <- ggplot(pp, aes(auroc, ct)) +
  geom_line(aes(group = ct), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(colour = method), size = 1.8) +
  scale_colour_manual(values = c("Gene base rates" = "grey55", "Base rates by subtype" = "#56B4E9",
                                 "Logistic model (whole profile)" = "#0072B2")) +
  labs(x = "AUROC on held-out patients", y = NULL, colour = NULL, title = "predict_partners(): beyond base rates") +
  theme(legend.position = "bottom", legend.direction = "vertical")

# c -- frequency agreement MSK vs other centres
cf <- rd("fig3d_centre_frequency.csv") %>% filter(MSK_n >= 200, other_n >= 200) %>% mutate(ct = short_ct(Cancer_Type))
r <- cf %>% group_by(ct) %>% summarise(r = cor(MSK, other), n = n())
pc <- ggplot(cf, aes(MSK, other)) +
  geom_abline(linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 0.8, alpha = 0.7, colour = "#0072B2", stroke = 0) +
  geom_text(data = r, aes(x = 0, y = Inf, label = sprintf("r = %.2f", r)), hjust = 0, vjust = 1.3, size = 2.2) +
  facet_wrap(~ct, nrow = 2) +
  scale_x_continuous(labels = percent, trans = "sqrt", breaks = c(0.01, 0.1, 0.3, 0.6)) +
  scale_y_continuous(labels = percent, trans = "sqrt", breaks = c(0.01, 0.1, 0.3, 0.6)) +
  coord_fixed() +
  labs(x = "Patients altered, MSK", y = "Patients altered, other centres",
       title = "Frequencies transfer across centres",
       caption = "One point per analysed gene; only patients whose panel covers the gene.")

# One combined figure for the phase; the title states the main message
fig <- (pa | (pb / pc + plot_layout(heights = c(1, 1.3)))) / s2 + plot_layout(heights = c(2.2, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = "Phase 3: the lookup tool's answers -- partners, exact combinations, predictions -- hold up")
save_fig(fig, "figure_phase3", 183, 235, OUT_DIR)
