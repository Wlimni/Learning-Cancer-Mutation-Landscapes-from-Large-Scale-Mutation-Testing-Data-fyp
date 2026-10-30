# Script Name: phase2c_methods.R
# Purpose : Phase 2c figure -- precision and recall of our test vs Fisher, DISCOVER-style and SELECT (notebook 02c, Section G)
#           a precision per method in three simulation scenarios; b recall by planted strength;
#           c precision-recall curves; d textbook pairs found in real data; e GENIE vs TCGA effect for every
#           candidate (independent patients).
# Inputs  : results/figures/data/fig5g_sim_precision_recall.csv, fig5h_sim_recall_by_strength.csv, fig5i_sim_pr.csv, fig5j_textbook.csv,
#           fig5a_tcga.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase2c_methods/figure_phase2c_methods.{pdf,png} -- one combined figure (5 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase2c_methods"
METHODS <- c("Fisher", "DISCOVER-style", "SELECT", "Ours")
COL_M <- c(Fisher = "#E69F00", `DISCOVER-style` = "#56B4E9", SELECT = "#CC79A7", Ours = "#0072B2")
SCEN <- c("S1 model-based" = "S1 model-based\n(our assumption)", "S2 real data, links removed" = "S2 real data,\nlinks removed",
          "S3 misspecified" = "S3 misspecified\n(stress test)")
REALISTIC <- "S2 real data, links removed"

# a -- precision at q < 0.05 (share of calls that are planted links, found in the right direction)
fd <- rd("fig5g_sim_precision_recall.csv") %>% mutate(method = factor(method, METHODS), scenario = factor(SCEN[scenario], SCEN))
pa <- ggplot(fd, aes(scenario, precision, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = percent(precision, 1)), position = position_dodge(0.8), angle = 90, hjust = -0.15, size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = COL_M) +
  scale_y_continuous(labels = percent, limits = c(0, 1.18), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "Precision\n(calls that are real links)", fill = NULL,
       title = "Precision: how many calls are real",
       caption = "Dashed line: 95% precision, what q < 0.05 aims for.") +
  theme(legend.position = "bottom")

# b -- recall: share of planted links found (q < 0.05, right direction), realistic scenario
pw <- rd("fig5h_sim_recall_by_strength.csv") %>% filter(scenario == REALISTIC) %>%
  mutate(method = factor(method, METHODS),
         direction = factor(direction, c("co-occurring", "exclusive")))
pb <- ggplot(pw, aes(fold, recall, colour = method)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  facet_wrap(~direction, labeller = as_labeller(c(`co-occurring` = "planted co-occurrence",
                                                  exclusive = "planted exclusivity"))) +
  scale_colour_manual(values = COL_M) +
  scale_x_continuous(breaks = c(1.5, 2, 3), labels = c("1.5x", "2x", "3x")) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "Strength of the planted link", y = "Recall\n(planted links found)", colour = NULL,
       title = "Recall: how many real links are found",
       caption = "Scenario S2 (real data, links removed); colours as in a.") +
  theme(legend.position = "none")

# c -- precision-recall over all thresholds (realistic scenario); dots = the q < 0.05 cut-off
pr <- rd("fig5i_sim_pr.csv") %>% filter(scenario == REALISTIC) %>% mutate(method = factor(method, METHODS))
op <- fd %>% filter(scenario == SCEN[[REALISTIC]])
pc <- ggplot(pr, aes(recall, precision, colour = method)) +
  geom_path(linewidth = 0.6) +
  geom_point(data = op, size = 2.4, shape = 21, fill = "white", stroke = 1) +
  scale_colour_manual(values = COL_M) +
  scale_x_continuous(labels = percent, limits = c(0, 1)) + scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "Recall (planted links found)", y = "Precision (calls that are real)", colour = NULL,
       title = "Precision and recall at every threshold",
       caption = "Circles: each method at its q < 0.05 cut-off. Scenario S2; colours as in a.") +
  theme(legend.position = "none")

# d -- textbook pairs (chosen from the literature beforehand) found in real data
tb <- rd("fig5j_textbook.csv") %>%
  mutate(found = found == "True", method = factor(method, METHODS),
         row = factor(paste0(pair, "  (", short_ct(Cancer_Type), ")"),
                      rev(unique(paste0(pair, "  (", short_ct(Cancer_Type), ")")))))
nf <- tb %>% group_by(method) %>% summarise(n = sum(found)) %>% mutate(lab = sprintf("%s\n%d/%d", method, n, n_distinct(tb$row)))
pd_ <- ggplot(tb, aes(method, row, fill = ifelse(found, expected, "not found"))) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = c(`co-occurring` = COL_DIR[[1]], exclusive = COL_DIR[[2]], `not found` = "grey88"),
                    breaks = c("co-occurring", "exclusive", "not found")) +
  scale_x_discrete(position = "top", labels = setNames(nf$lab, nf$method)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Recall of textbook pairs (real data)",
       caption = "29 pairs with a published direction (data/curated/textbook_pairs.csv).") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6, face = "italic"), axis.line = element_blank(),
        axis.ticks = element_blank())

# e -- independent patients: every candidate testable in TCGA, effect in GENIE vs TCGA
tv <- rd("fig5a_tcga.csv") %>% mutate(across(c(interaction_candidate, replicated), ~ .x == "True"))
kb <- tv %>% filter(interaction_candidate) %>%
  mutate(GENIE = log2((n_both_altered + 0.5) / (expected_both + 0.5)), TCGA = log2((tcga_n_both + 0.5) / (tcga_expected + 0.5)))
rho <- cor(kb$GENIE, kb$TCGA, method = "spearman")
same_c <- mean(kb$tcga_direction == kb$direction)
ctl <- tv %>% filter(q_value >= 0.05); same_0 <- mean(ctl$tcga_direction == ctl$direction)
pe <- ggplot(kb, aes(GENIE, TCGA, colour = direction)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_abline(linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 0.9, alpha = 0.45) +
  scale_colour_manual(values = COL_DIR[1:2]) +
  labs(x = "GENIE: log2(observed / expected)", y = "TCGA: log2(observed / expected)", colour = NULL,
       title = "Independent patients (TCGA) agree",
       caption = sprintf("%d candidates; Spearman rho = %.2f. Same direction: %.0f%% (control pairs: %.0f%%).",
                         nrow(kb), rho, 100 * same_c, 100 * same_0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(legend.position = "bottom")

# One combined figure: simulation on the left, real data on the right; the title states the main message
fig <- ((pa / pb / pc) | (pd_ / pe + plot_layout(heights = c(1.6, 1)))) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = "Methods compared: the highest precision, and the most real links found at any given precision")
save_fig(fig, "figure_phase2c_methods", 183, 240, OUT_DIR)
