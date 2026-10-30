# Script Name: phase2b_validation.R
# Purpose : Phase 2b figure -- which alteration set gives the better gene-pair list (notebook 02b):
#           a published pairs found per set; b precision and recall on data with planted links; c TCGA direction.
# Inputs  : results/figures/data/fig2b_known.csv, fig2b_sim.csv, fig2b_tcga.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase2b_validation/figure_phase2b.{pdf,png} -- one combined figure (3 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase2b_validation"
SETS <- c("SNV", "SNV+CNV", "SNV+CNV+SV")
COL_SET <- c(SNV = "#56B4E9", `SNV+CNV` = "#0072B2", `SNV+CNV+SV` = "#CC79A7")
VERD <- c("found" = "#009E73", "removed: one DNA event" = "grey45", "not significant" = "grey75",
          "not tested" = "grey90", "opposite direction" = "#D55E00")

# a -- published pairs: found, or why not
kn <- rd("fig2b_known.csv") %>% mutate(alteration_set = factor(alteration_set, rev(SETS)), verdict = factor(verdict, rev(names(VERD))))
nk <- n_distinct(paste(kn$Cancer_Type, kn$pair))
lab <- kn %>% group_by(alteration_set) %>% summarise(f = sum(verdict == "found"))
pa <- ggplot(kn, aes(y = alteration_set, fill = verdict)) +
  geom_bar(width = 0.65) +
  geom_text(data = lab, aes(x = nk + 0.5, y = alteration_set, label = sprintf("%d of %d found", f, nk)),
            inherit.aes = FALSE, hjust = 0, size = 2.6) +
  scale_fill_manual(values = VERD, breaks = names(VERD)) +
  scale_x_continuous(expand = expansion(c(0, 0.35))) +
  labs(x = "Published gene pairs", y = NULL, fill = NULL, title = "Known biology found",
       caption = "Pairs with a published direction, fixed before looking at the results.\nNot tested = a gene (or its copy-number change / fusion) is not in that set.") +
  guides(fill = guide_legend(nrow = 2)) + theme(legend.position = "bottom")

# b -- precision and recall on data with planted links
sm <- rd("fig2b_sim.csv") %>% mutate(alteration_set = factor(alteration_set, SETS))
pb <- ggplot(sm, aes(measure, value, fill = alteration_set)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = percent(value, 1)), position = position_dodge(0.75), vjust = -0.4, size = 2.4) +
  scale_fill_manual(values = COL_SET) +
  scale_y_continuous(labels = percent, limits = c(0, 1.08), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Planted links: precision and recall",
       caption = "Real data, links shuffled away, ~215 links of known strength planted (q < 0.05).") +
  theme(legend.position = "bottom")

# c -- TCGA: candidates keep their direction in independent patients
tc <- rd("fig2b_tcga.csv") %>% mutate(alteration_set = factor(alteration_set, SETS))
pc <- ggplot(tc, aes(alteration_set, same_direction, fill = group)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(n=%s)", percent(same_direction, 1), comma(n))), position = position_dodge(0.75),
            vjust = -0.25, size = 2.2, lineheight = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = c(candidates = "#0072B2", `control: not significant` = "grey70")) +
  scale_y_continuous(labels = percent, limits = c(0, 1.15), expand = c(0, 0)) +
  labs(x = NULL, y = "Same direction in TCGA", fill = NULL, title = "Independent patients (TCGA)",
       caption = "TCGA has no fusion calls; dashed line = chance.") +
  theme(legend.position = "bottom")

fig <- (pa / (pb | pc)) + plot_layout(heights = c(0.8, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = "Phase 2b: equally precise on all sets; copy number and fusions add known biology")
save_fig(fig, "figure_phase2b", 183, 150, OUT_DIR)
