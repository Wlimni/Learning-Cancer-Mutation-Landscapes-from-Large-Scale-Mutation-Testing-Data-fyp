# Script Name: fig6_tool.R
# Purpose : Figure 6 -- the lookup tool (Aim 1).
#           a lookup("KRAS G12C", NSCLC): partners, observed vs expected (the tool's main output);
#           b partner prediction: logistic model vs base rates per cancer type (SelectSim-style benchmark);
#           c do the tool's frequencies transfer between centres? MSK vs other centres (GENIE 2017 Fig 3 idea).
# Inputs  : results/figures/data/fig3a_kras_g12c.csv, fig3c_predict_partners.csv, fig3d_centre_frequency.csv
# Outputs : results/figures/fig6_tool.{pdf,png}
source("scripts/figures/_theme.R")

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

fig <- (pa | (pb / pc + plot_layout(heights = c(1, 1.3)))) + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "a")
save_fig(fig, "fig6_tool", 183, 170)
