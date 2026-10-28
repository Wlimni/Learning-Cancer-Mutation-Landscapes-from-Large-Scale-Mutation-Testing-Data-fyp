# Script Name: phase2b_validation.R
# Purpose : Phase 2b figure -- can the Phase 2 results be trusted? (notebook 02b)
#           a TCGA: how many pairs replicate in independent patients; b GENIE vs TCGA effect for every
#           testable candidate; c the method benchmark; d agreement with the literature and the SelectSim atlas.
# Inputs  : results/figures/data/fig5a_tcga.csv, fig2a_null_pvalues.csv, fig5b_benchmark_real.csv,
#           fig5f_known.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase2b_validation/figure_phase2b.{pdf,png} -- one combined figure (4 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase2b_validation"

# a -- TCGA: replication in independent patients, with plain counts
tv <- rd("fig5a_tcga.csv") %>% mutate(across(c(interaction_candidate, replicated), ~ .x == "True"))
grp <- list("Replicated candidates" = tv$interaction_candidate & tv$replicated,
            "All candidates" = tv$interaction_candidate,
            "Not significant in GENIE\n(control)" = tv$q_value >= 0.05)
ta <- bind_rows(lapply(names(grp), function(g) {
  d <- tv[grp[[g]], ]; same <- d$tcga_direction == d$direction
  data.frame(group = g, n = nrow(d),
             metric = c("Same direction", "Significant in TCGA,\nsame direction"),
             value = c(mean(same), mean(d$tcga_p < 0.05 & same)))
})) %>% mutate(group = factor(sprintf("%s\n(n = %d)", group, n), unique(sprintf("%s\n(n = %d)", group, n))))
opp <- with(tv[grp[[1]], ], sum(tcga_p < 0.05 & tcga_direction != direction))
pa <- ggplot(ta, aes(group, value, fill = metric)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = percent(value, 1)), position = position_dodge(0.75), vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("#56B4E9", "#0072B2")) +
  scale_y_continuous(labels = percent, limits = c(0, 1.08), expand = c(0, 0)) +
  labs(x = NULL, y = "Pairs testable in TCGA", fill = NULL,
       title = "TCGA (independent patients) confirms the pairs",
       caption = sprintf("Replicated candidates significant in the opposite direction in TCGA: %d.", opp)) +
  theme(legend.position = "bottom")

# b -- every candidate testable in TCGA: effect in GENIE vs TCGA, eight well-known pairs labelled
KNOWN <- c("Non-Small Cell Lung Cancer|KEAP1|STK11", "Glioma|ATRX|IDH1", "Glioma|IDH1|TP53", "Glioma|CDKN2A|EGFR",
           "Breast Cancer|PIK3CA|PTEN", "Colorectal Cancer|BRAF|KRAS", "Melanoma|BRAF|NRAS", "Non-Small Cell Lung Cancer|EGFR|KRAS")
kb <- tv %>% filter(interaction_candidate) %>%
  mutate(GENIE = log2((n_both_altered + 0.5) / (expected_both + 0.5)), TCGA = log2((tcga_n_both + 0.5) / (tcga_expected + 0.5)),
         key = paste(Cancer_Type, Gene_A, Gene_B, sep = "|"),
         lab = ifelse(key %in% KNOWN, paste0(Gene_A, "+", Gene_B), NA))
rho <- cor(kb$GENIE, kb$TCGA, method = "spearman")
pb <- ggplot(kb, aes(GENIE, TCGA, colour = direction)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_abline(linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 0.9, alpha = 0.45) +
  geom_point(data = filter(kb, !is.na(lab)), size = 1.8, colour = "black", shape = 21, stroke = 0.6) +
  geom_label_repel(aes(label = lab), size = 2.1, colour = "black", fontface = "italic", na.rm = TRUE,
                   fill = alpha("white", 0.85), label.size = 0, label.padding = 0.1, min.segment.length = 0,
                   box.padding = 0.6, force = 4, max.overlaps = Inf, seed = 3) +
  scale_colour_manual(values = COL_DIR[1:2]) +
  scale_x_continuous(expand = expansion(c(0.05, 0.12))) + scale_y_continuous(expand = expansion(c(0.05, 0.25))) +
  labs(x = "GENIE: log2(observed / expected)", y = "TCGA: log2(observed / expected)", colour = NULL,
       title = "Effects agree across the two cohorts",
       caption = sprintf("%d candidates testable in both; Spearman rho = %.2f. Dashed line = equal effect.\nCircled: well-known pairs.", nrow(kb), rho)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(legend.position = "bottom")

# c -- method benchmark: false alarms on null data, and TCGA replication of each method's calls
M <- c(p_fisher = "Fisher", p_discover = "DISCOVER-style", p_select = "SELECT", p_cond_subtype = "Ours")
DIRS <- c(p_fisher = "dir_fisher", p_discover = "dir_discover", p_select = "dir_select", p_cond_subtype = "dir_cond_subtype")
qv <- function(df, pc) if (pc == "p_select") df$q_select else p.adjust(pmax(df[[pc]], 1e-300), "BH")
nl <- rd("fig2a_null_pvalues.csv"); real <- rd("fig5b_benchmark_real.csv")
bc <- bind_rows(lapply(names(M), function(pc) {
  s <- real[qv(real, pc) < 0.05 & !is.na(real$tcga_p), ]
  data.frame(method = M[[pc]],
             metric = c("False alarms on data with\nno real links (target 5%)", "Significant calls confirmed\nin TCGA"),
             value = c(mean(nl[[pc]] < 0.05), mean(s$tcga_p < 0.05 & s$tcga_direction == s[[DIRS[[pc]]]])))
})) %>% mutate(method = factor(method, M))
pc <- ggplot(bc, aes(method, value, fill = method)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = percent(value, 1)), vjust = -0.3, size = 2.4) +
  geom_hline(data = data.frame(metric = unique(bc$metric)[1], y = 0.05), aes(yintercept = y), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c(Fisher = "#E69F00", `DISCOVER-style` = "#56B4E9", SELECT = "#CC79A7", Ours = "#0072B2"), guide = "none") +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.18))) +
  labs(x = NULL, y = NULL, title = "Our test vs published methods, on identical data") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# d -- agreement with prior knowledge: literature (Europe PMC) and the SelectSim atlas, robust pairs
kn <- rd("fig5f_known.csv")
LV <- c("well known", "some reports", "rarely or never reported", "same type, same direction", "pan-cancer only", "not reported")
kd <- bind_rows(count(kn, status = novelty) %>% mutate(source = "Literature"),
                count(kn, status = selectsim) %>% mutate(source = "SelectSim atlas\n(Nat Genet 2026)")) %>%
  mutate(status = factor(status, LV), source = factor(source, rev(unique(source))))
pd_ <- ggplot(kd, aes(n, source, fill = status)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.5, colour = "white") +
  scale_fill_manual(values = c("well known" = "#08519C", "some reports" = "#4292C6", "rarely or never reported" = "grey65",
                               "same type, same direction" = "#08519C", "pan-cancer only" = "#4292C6", "not reported" = "grey65")) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = sprintf("Robust pairs (n = %d)", nrow(kn)), y = NULL, fill = NULL,
       title = "Most pairs are known biology",
       caption = "SelectSim: none in the opposite direction. It covers mutations only,\nso copy-number and fusion pairs cannot appear there.") +
  guides(fill = guide_legend(ncol = 2)) +
  theme(legend.position = "bottom")

# One combined figure; the title states the main message
fig <- (pa | pb) / (pc | pd_) + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = "Phase 2b: the pairs replicate in new patients and match known biology")
save_fig(fig, "figure_phase2b", 183, 175, OUT_DIR)
