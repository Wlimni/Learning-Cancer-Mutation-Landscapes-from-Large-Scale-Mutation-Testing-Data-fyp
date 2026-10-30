# Script Name: phase2_pair_test.R
# Purpose : Phase 2 figure -- the gene-pair list and its size (notebook 02). Evaluation is in figure 2b.
#           a the three alteration sets: candidates by direction, and how many replicate outside MSK;
#           b from tested pairs to replicated candidates (main set); c candidates per cancer type (main set).
# Inputs  : results/figures/data/fig2a_sets.csv, fig2f_funnel.csv, fig2_pairs.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase2_pair_test/figure_phase2.{pdf,png} -- one combined figure (3 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase2_pair_test"
DIRS <- c("co-occurring", "exclusive")

# a -- the three alteration sets
st <- rd("fig2a_sets.csv")
st$set <- factor(st$`alteration set`, rev(st$`alteration set`))
sl <- st %>% select(set, `co-occurring` = `candidates co-occurring`, exclusive = `candidates exclusive`) %>%
  pivot_longer(-set, names_to = "direction", values_to = "n")
pa <- ggplot(sl, aes(n, set, fill = direction)) +
  geom_col(width = 0.6) +
  geom_text(data = st, aes(candidates, set, label = sprintf("%s (%s replicated)", comma(candidates), comma(`replicated candidates`))),
            inherit.aes = FALSE, hjust = -0.05, size = 2.4) +
  scale_fill_manual(values = COL_DIR[DIRS]) +
  scale_x_continuous(labels = comma, expand = expansion(c(0, 0.85))) +
  labs(x = "Interaction candidates (summed over cancer types)", y = NULL, fill = NULL,
       title = "Three definitions of 'altered'",
       caption = sprintf("Pairs tested: %s.\nMSK hits replicated in other centres: %s.",
                         paste(sprintf("%s %s", st$`alteration set`, comma(st$`pairs tested`)), collapse = ", "),
                         paste(sprintf("%s %.0f%%", st$`alteration set`, st$`MSK hits replicated %`), collapse = ", "))) +
  theme(legend.position = "bottom")

# b -- funnel, main set
fn0 <- rd("fig2f_funnel.csv") %>% mutate(step = factor(step, rev(step)))
fn <- fn0 %>% pivot_longer(c(co_occurring, exclusive), names_to = "direction") %>% mutate(direction = sub("_", "-", direction))
pb <- ggplot(fn, aes(value, step, fill = direction)) +
  geom_col(width = 0.65) +
  geom_text(data = fn0, aes(x = pairs, y = step, label = ifelse(pct_of_tested >= 100, comma(pairs), sprintf("%s (%.1f%%)", comma(pairs), pct_of_tested))),
            inherit.aes = FALSE, hjust = -0.1, size = 2.4) +
  scale_fill_manual(values = COL_DIR[DIRS], guide = "none") +
  scale_x_continuous(labels = comma, expand = expansion(c(0, 0.35))) +
  labs(x = "Gene pairs (summed over cancer types)", y = NULL, title = "Main set (SNV+CNV): from tested pairs to candidates")

# c -- candidates per cancer type, main set
pr <- rd("fig2_pairs.csv") %>% mutate(across(c(interaction_candidate, replicated), ~ .x == "True")) %>% filter(interaction_candidate)
top <- pr %>% count(Cancer_Type, sort = TRUE) %>% slice_head(n = 15)
pc_d <- pr %>% filter(Cancer_Type %in% top$Cancer_Type) %>% count(Cancer_Type, direction) %>%
  mutate(ct = factor(short_ct(Cancer_Type), short_ct(rev(top$Cancer_Type))))
rep_ct <- pr %>% filter(Cancer_Type %in% top$Cancer_Type) %>% group_by(Cancer_Type) %>%
  summarise(n = n(), rep = sum(replicated)) %>% mutate(ct = factor(short_ct(Cancer_Type), levels(pc_d$ct)))
pc <- ggplot(pc_d, aes(n, ct, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(data = rep_ct, aes(n, ct, label = sprintf("%d (%d replicated)", n, rep)), inherit.aes = FALSE, hjust = -0.08, size = 2.3) +
  scale_fill_manual(values = COL_DIR[DIRS], guide = "none") +
  scale_x_continuous(expand = expansion(c(0, 0.7))) +
  labs(x = "Interaction candidates", y = NULL, title = "Where the candidates are (main set)",
       caption = sprintf("15 of %d cancer types with candidates.", n_distinct(pr$Cancer_Type))) +
  theme(axis.text.y = element_text(size = 6.5))

fig <- ((pa / pb) | pc) + plot_layout(widths = c(1.1, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = sprintf("Phase 2: %s interaction candidates in the main list; %s replicate outside MSK",
                                  comma(st$candidates[st$`alteration set` == "SNV+CNV"]),
                                  comma(st$`replicated candidates`[st$`alteration set` == "SNV+CNV"])))
save_fig(fig, "figure_phase2", 183, 120, OUT_DIR)
