# Script Name: fig2_method_validity.R
# Purpose : Figure 2 -- can the pair test be trusted?
#           a p-value histograms on null data with no gene-gene links (DISCOVER 2016 Fig 2);
#           b false-positive rate and TCGA replication per method (CoMEt / SelectSim benchmark bars);
#           c GENIE vs TCGA effect per pair (Park & Lehner 2015 EV3; SelectSim Fig 1j).
# Inputs  : results/figures/data/fig2a_null_pvalues.csv, fig5b_benchmark_real.csv, fig5a_tcga.csv
# Outputs : results/figures/fig2_method_validity.{pdf,png}
source("scripts/figures/_theme.R")

M <- c(p_fisher = "Fisher", p_discover = "DISCOVER-style", p_select = "SELECT", p_conditional = "Conditional",
       p_cond_subtype = "Conditional + subtype")
DIRS <- c(p_fisher = "dir_fisher", p_discover = "dir_discover", p_select = "dir_select", p_conditional = "dir_conditional",
          p_cond_subtype = "dir_cond_subtype")
# SELECT reports its own empirical FDR (q_select); the other methods use Benjamini-Hochberg
qv <- function(df, pc) if (pc == "p_select") df$q_select else p.adjust(pmax(df[[pc]], 1e-300), "BH")

# a -- null p-value histograms: a calibrated test is flat (or below the line near 0 if conservative)
nl <- rd("fig2a_null_pvalues.csv") %>% select(all_of(names(M))) %>%
  pivot_longer(everything(), names_to = "m", values_to = "p") %>%
  mutate(method = factor(M[m], M))
lab <- nl %>% group_by(method) %>% summarise(fp = mean(p < 0.05)) %>%
  mutate(txt = sprintf("p < 0.05: %.1f%%", 100 * fp))
pa <- ggplot(nl, aes(p, fill = method)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.05), colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_text(data = lab, aes(x = 0.5, y = Inf, label = txt), vjust = 1.5, size = 2.4, inherit.aes = FALSE) +
  facet_wrap(~method, nrow = 1) +
  theme(strip.text = element_text(size = 6.5)) +
  scale_fill_manual(values = COL_METHOD, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(x = "p-value on null data (no real gene-gene links)", y = "Density",
       title = "Calibration: same patients, links removed by curveball swaps (SELECT re-implemented from Mina 2017)",
       caption = sprintf("%s null pairs in 9 cancer types; dashed line = ideal flat distribution. Spike at 0 = false positives;\nmass at 1 = conservative (discrete) test.",
                         comma(nrow(nl) / length(M))))

# b -- false positives on null data vs replication of significant calls in TCGA
real <- rd("fig5b_benchmark_real.csv")
rep <- lapply(names(M), function(pc) {
  q <- qv(real, pc)
  s <- real[q < 0.05 & !is.na(real$tcga_p), ]
  data.frame(method = M[[pc]], n_sig = sum(q < 0.05),
             value = mean(s$tcga_p < 0.05 & s$tcga_direction == s[[DIRS[[pc]]]]))
}) %>% bind_rows()
bdat <- bind_rows(lab %>% transmute(method = as.character(method), value = fp, metric = "False positives on null data\n(p < 0.05; target 5%)"),
                  rep %>% transmute(method, value, metric = "Significant calls replicating\nin TCGA (same direction)")) %>%
  mutate(method = factor(method, M))
pb <- ggplot(bdat, aes(method, value, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = percent(value, 0.1)), vjust = -0.3, size = 1.9) +
  geom_hline(data = data.frame(metric = unique(bdat$metric)[1], y = 0.05), aes(yintercept = y),
             linetype = "dashed", linewidth = 0.3, colour = "grey55") +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = COL_METHOD, guide = "none") +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.15))) +
  labs(x = NULL, y = NULL, title = "Head-to-head on identical data") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# c -- external validation: GENIE effect vs TCGA effect
tv <- rd("fig5a_tcga.csv") %>%
  mutate(genie = log2((n_both_altered + 0.5) / (expected_both + 0.5)),
         tcga = log2((tcga_n_both + 0.5) / (tcga_expected + 0.5)),
         cls = case_when(interaction_candidate == "True" & replicated == "True" ~ "Replicated candidate",
                         q_value >= 0.05 ~ "Not significant in GENIE", TRUE ~ "Other"))
tv <- tv %>% filter(cls != "Other") %>% arrange(desc(cls))
rc <- tv %>% filter(cls == "Replicated candidate")
stats <- sprintf("Replicated candidates (n = %d): %.0f%% same direction,\n%.0f%% significant in TCGA, %d significant in the opposite direction",
                 nrow(rc), 100 * mean(rc$tcga_direction == rc$direction),
                 100 * mean(rc$tcga_p < 0.05 & rc$tcga_direction == rc$direction),
                 sum(rc$tcga_p < 0.05 & rc$tcga_direction != rc$direction))
known <- tv %>% filter(paste(Cancer_Type, Gene_A, Gene_B) %in%
                         c("Non-Small Cell Lung Cancer KEAP1 STK11", "Non-Small Cell Lung Cancer EGFR KRAS",
                           "Breast Cancer CCND1 FGFR1", "Glioma CDKN2A EGFR")) %>%
  mutate(lbl = paste0(Gene_A, "/", Gene_B, " (", short_ct(Cancer_Type), ")"))
pc <- ggplot(tv, aes(genie, tcga)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey50") +
  geom_point(aes(colour = cls, size = cls), alpha = 0.7, stroke = 0) +
  geom_text_repel(data = known, aes(label = lbl), size = 2.3, min.segment.length = 0, box.padding = 0.4, seed = 1) +
  scale_colour_manual(values = c("Replicated candidate" = "#0072B2", "Not significant in GENIE" = "grey75")) +
  scale_size_manual(values = c("Replicated candidate" = 1.3, "Not significant in GENIE" = 0.6), guide = "none") +
  coord_cartesian(xlim = c(-5, 4), ylim = c(-5, 4)) +
  labs(x = "GENIE effect, log2(observed / expected)", y = "TCGA effect, log2(observed / expected)", colour = NULL,
       title = "Independent patients: TCGA PanCancer Atlas", caption = stats) +
  theme(legend.position = c(0.3, 0.92))

fig <- pa / ((pb | pc) + plot_layout(widths = c(1.25, 1))) + plot_layout(heights = c(1, 1.6)) + plot_annotation(tag_levels = "a")
save_fig(fig, "fig2_method_validity", 183, 150)
