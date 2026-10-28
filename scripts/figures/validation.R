# Script Name: validation.R
# Purpose : Validation figures -- can the results be trusted? (notebook 05).
# Inputs  : results/figures/data/fig2a*, fig5*.csv (from 00_export_figure_data.py)
# Outputs : results/figures/validation/<step>_<name>.{pdf,png} -- one figure per step, shown in the phase's notebook
source("scripts/figures/_theme.R")
F <- "validation"

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
save_fig(pa / pb + plot_layout(heights = c(1, 1.3)), "5.1_method_benchmark", 183, 140, F)

# c -- external validation: GENIE effect vs TCGA effect
tv <- rd("fig5a_tcga.csv") %>%
  mutate(genie = log2((n_both_altered + 0.5) / (expected_both + 0.5)),
         tcga = log2((tcga_n_both + 0.5) / (tcga_expected + 0.5)),
         cls = case_when(interaction_candidate == "True" ~ "Interaction candidate",
                         q_value >= 0.05 ~ "Not significant in GENIE", TRUE ~ "Other"))
tv <- tv %>% filter(cls != "Other") %>% arrange(desc(cls))
rc <- tv %>% filter(cls == "Interaction candidate")
stats <- sprintf("Interaction candidates testable in TCGA (n = %d): %.0f%% same direction,\n%.0f%% significant in TCGA, %d significant in the opposite direction",
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
  scale_colour_manual(values = c("Interaction candidate" = "#0072B2", "Not significant in GENIE" = "grey75")) +
  scale_size_manual(values = c("Interaction candidate" = 1.3, "Not significant in GENIE" = 0.6), guide = "none") +
  coord_cartesian(xlim = c(-5, 4), ylim = c(-5, 4)) +
  labs(x = "GENIE effect, log2(observed / expected)", y = "TCGA effect, log2(observed / expected)", colour = NULL,
       title = "Independent patients: TCGA PanCancer Atlas", caption = stats) +
  theme(legend.position = c(0.3, 0.92))
save_fig(pc, "5.2_tcga_agreement", 110, 105, F)

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
save_fig(s5, "5.3_dnds_route", 110, 95, F)

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
save_fig(s6, "5.4_exclusive_gene_sets", 110, 80, F)

# 5.5 -- mutation-load confounding (van de Haar 2019 diagnostic): does a gene's load score (MLA)
# predict how many exclusive pairs a method calls for it? Strongly negative = confounded by load.
L <- c(Fisher = "Fisher", `DISCOVER-style` = "DISCOVER-style", SELECT = "SELECT", `Conditional + subtype` = "Conditional + subtype")
corr <- function(d, design) bind_rows(lapply(names(L), function(m) {
  col <- paste0("ex_", m)
  if (!col %in% names(d)) return(NULL)
  data.frame(method = L[[m]], design = design, r = cor(d$MLA, d[[col]]))
}))
ld <- bind_rows(corr(rd("fig5e_load.csv"), "Non-hypermutated only (our design)"),
                corr(rd("fig5e_load_hyper.csv"), "Hypermutated tumours included")) %>%
  mutate(method = factor(method, L), design = factor(design, c("Non-hypermutated only (our design)", "Hypermutated tumours included")))
p55 <- ggplot(ld, aes(method, r, fill = design)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.87, ymax = -0.60, fill = "grey85", alpha = 0.6) +
  annotate("text", x = 0.55, y = -0.735, label = "range reported by van de Haar 2019 (exome data)", hjust = 0, size = 2.2, colour = "grey30") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  scale_y_continuous(limits = c(-0.9, 0.3)) +
  labs(x = NULL, y = "Correlation: load score vs exclusive calls per gene", fill = NULL,
       title = "Mutation-load confounding is removed by the hypermutation split",
       caption = "Negative = genes mutated mostly in low-load tumours collect extra exclusive calls (confounding).\nSELECT was not re-run with hypermutated tumours.") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
save_fig(p55, "5.5_mutation_load_check", 120, 100, F)
