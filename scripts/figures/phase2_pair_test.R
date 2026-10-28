# Script Name: phase2_pair_test.R
# Purpose : Phase 2 figures -- testing gene pairs fairly, and what it finds (notebook 02).
# Inputs  : results/figures/data/fig2*.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase2_pair_test/figure_phase2.{pdf,png} -- one combined figure (4 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase2_pair_test"

# 2.1 -- why not Fisher: on data with no real links, how many pairs does each test call significant?
nl <- rd("fig2a_null_pvalues.csv")
fp <- data.frame(test = factor(c("Fisher's test\n(usual)", "Conditional test\n(ours)"), c("Fisher's test\n(usual)", "Conditional test\n(ours)")),
                 fp = c(mean(nl$p_fisher < 0.05), mean(nl$p_cond_subtype < 0.05)))
p21 <- ggplot(fp, aes(test, fp, fill = test)) +
  geom_col(width = 0.55) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = 2.45, y = 0.05, label = "expected by chance: 5%", vjust = -0.5, hjust = 1, size = 2.4) +
  geom_text(aes(label = percent(fp, 0.1)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  scale_y_continuous(labels = percent, limits = c(0, 0.14), expand = c(0, 0)) +
  labs(x = NULL, y = "Pairs called significant (p < 0.05)", title = "Why not Fisher's test",
       caption = "Data with all real gene-gene links removed: every 'significant'\npair is a false alarm. Fisher over-calls; ours stays below 5%.")

# 2.2 -- filter funnel: from tested pairs to interaction candidates
fn0 <- rd("fig2f_funnel.csv") %>% filter(!grepl("replicated", step)) %>% mutate(step = factor(step, rev(step)))
fn <- fn0 %>% pivot_longer(c(co_occurring, exclusive), names_to = "direction") %>% mutate(direction = sub("_", "-", direction))
p22 <- ggplot(fn, aes(value, step, fill = direction)) +
  geom_col(width = 0.65) +
  geom_text(data = fn0, aes(x = pairs, y = step, label = ifelse(pct_of_tested >= 100, comma(pairs), sprintf("%s (%.1f%%)", comma(pairs), pct_of_tested))),
            inherit.aes = FALSE, hjust = -0.1, size = 2.4) +
  scale_fill_manual(values = COL_DIR[1:2]) +
  scale_x_continuous(labels = comma, expand = expansion(c(0, 0.45))) +
  labs(x = "Gene pairs (summed over cancer types)", y = NULL, fill = NULL, title = "From tested pairs to candidates") +
  theme(legend.position = "bottom")

pr <- rd("fig2_pairs.csv") %>%
  mutate(across(c(interaction_candidate, replicated), ~ .x == "True"),
         l2 = pmax(pmin(log2_fold_enrichment, 2), -2),
         cls = case_when(interaction_candidate ~ direction, q_value < 0.05 ~ "shared copy-number event",
                         TRUE ~ "not significant"))

# 2.5 -- the strongest interaction candidates (not weak): 10 most significant per direction
st <- pr %>% filter(interaction_candidate, weak_effect != "True") %>% group_by(direction) %>%
  slice_min(q_value, n = 10, with_ties = FALSE) %>% ungroup() %>%
  mutate(lab = paste0(Gene_A, " / ", Gene_B, "  (", short_ct(Cancer_Type), ")"),
         lab = factor(lab, lab[order(log2_fold_enrichment)]),
         obs = sprintf("%d / %.0f", n_both_altered, expected_both))
p25 <- ggplot(st, aes(log2_fold_enrichment, lab, colour = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2_fold_enrichment, yend = lab), linewidth = 0.5) +
  geom_point(size = 2) +
  geom_text(aes(x = max(log2_fold_enrichment) + 2.4, label = obs), colour = "grey30", size = 2.2, hjust = 1) +
  annotate("text", x = max(st$log2_fold_enrichment) + 2.4, y = nrow(st) + 0.9, label = "seen / expected",
           size = 2.2, hjust = 1, fontface = "bold") +
  scale_colour_manual(values = COL_DIR[1:2]) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.5))) +
  labs(x = "log2(observed / expected)", y = NULL, colour = NULL,
       title = "The strongest gene pairs",
       caption = "10 most significant candidates per direction (>= 1.5-fold).") +
  theme(axis.text.y = element_text(size = 6.5), legend.position = "bottom")

# 2.6 -- replication across centres: pairs significant in MSK, re-tested in all other centres
rp <- rd("fig2_pairs.csv") %>% filter(disc_q_value < 0.05, !is.na(rep_p_value)) %>%
  mutate(outcome = case_when(rep_p_value < 0.05 & rep_direction == disc_direction ~ "Confirmed (p < 0.05, same direction)",
                             rep_p_value < 0.05 ~ "Significant in the opposite direction",
                             rep_direction == disc_direction ~ "Same direction, not significant",
                             TRUE ~ "Opposite direction, not significant"),
         outcome = factor(outcome, c("Confirmed (p < 0.05, same direction)", "Same direction, not significant",
                                     "Opposite direction, not significant", "Significant in the opposite direction")))
rs <- rp %>% count(disc_direction, outcome) %>% group_by(disc_direction) %>% mutate(share = n / sum(n))
p26 <- ggplot(rs, aes(share, disc_direction, fill = outcome)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(share > 0.06, n, "")), position = position_stack(vjust = 0.5), size = 2.4, colour = "white") +
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "grey70", "#D55E00")) +
  scale_x_continuous(labels = percent, expand = expansion(c(0, 0.03))) +
  labs(x = "Pairs significant in MSK patients", y = NULL, fill = NULL,
       title = sprintf("Replication: %d of %d MSK findings confirmed in other centres (%.0f%%)",
                       sum(rp$outcome == levels(rp$outcome)[1]), nrow(rp), 100 * mean(rp$outcome == levels(rp$outcome)[1])),
       caption = sprintf("Each pair re-tested in patients from all other centres (different panels and pipelines); %d significant in the opposite direction.",
                         sum(rp$outcome == levels(rp$outcome)[4]))) +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")

# One combined figure for the phase; the title states the main message
fig <- ((p21 / p22) | p25) / p26 + plot_layout(heights = c(2.6, 0.62)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = "Phase 2: a test that accounts for each patient's alteration load finds known biology that replicates")
save_fig(fig, "figure_phase2", 183, 185, OUT_DIR)
