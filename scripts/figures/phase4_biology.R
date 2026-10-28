# Script Name: phase4_biology.R
# Purpose : Phase 4 figure -- what the patterns say about biology (notebook 04).
#           a exclusivity follows pathways; b the same pair reverses between cancer types (CDKN2A/TP53);
#           c KRAS G12C partners before and after the smoking adjustment; d the rarest combinations.
# Inputs  : results/figures/data/fig2_pairs.csv, fig4a_pathways.csv, fig4c_shrunk.csv, fig4d_smoking.csv,
#           fig4f_rare.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase4_biology/figure_phase4.{pdf,png} -- one combined figure (4 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase4_biology"

# a -- share of tested pairs that are exclusive / co-occurring candidates, same vs different pathway
pr <- rd("fig2_pairs.csv") %>% mutate(interaction_candidate = interaction_candidate == "True")
pw <- rd("fig4a_pathways.csv")
gp <- split(pw$pathway, pw$gene)
wb <- pr %>% filter(Gene_A %in% names(gp), Gene_B %in% names(gp)) %>%
  mutate(same = mapply(function(a, b) length(intersect(gp[[a]], gp[[b]])) > 0, Gene_A, Gene_B)) %>%
  group_by(same) %>%
  summarise(exclusive = mean(interaction_candidate & direction == "exclusive"),
            `co-occurring` = mean(interaction_candidate & direction == "co-occurring")) %>%
  pivot_longer(c(exclusive, `co-occurring`), names_to = "direction") %>%
  mutate(where = ifelse(same, "Same pathway", "Different pathways"))
pa <- ggplot(wb, aes(where, value, fill = direction)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = percent(value, 1)), position = position_dodge(0.75), vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = COL_DIR[1:2]) +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.15))) +
  labs(x = NULL, y = "Tested pairs that are candidates", fill = NULL, title = "Exclusivity follows pathways",
       caption = "Genes in the same pathway avoid each other 2.5x more often:\none hit per pathway is enough.") +
  theme(legend.position = c(0.3, 0.9))

# b -- one pair, opposite directions in different cancer types (cancer types where it is significant)
sh <- rd("fig4c_shrunk.csv") %>% filter(key == "CDKN2A/TP53", q_value < 0.05) %>%
  mutate(l2 = theta / log(2), se = sqrt(var) / log(2), ct = reorder(short_ct(Cancer_Type), l2))
pb <- ggplot(sh, aes(l2, ct, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(xmin = l2 - 1.96 * se, xmax = l2 + 1.96 * se), width = 0, orientation = "y", linewidth = 0.4) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = COL_DIR[1:2]) +
  labs(x = "log2(observed / expected), 95% CI", y = NULL, colour = NULL,
       title = "CDKN2A + TP53: together in some cancers, apart in others",
       caption = "Cancer types where the pair is significant (q < 0.05).") +
  theme(legend.position = "bottom")

# c -- KRAS G12C partners: crude vs smoking-adjusted
sm <- rd("fig4d_smoking.csv") %>% filter(allele == "KRAS_G12C") %>%
  mutate(partner = reorder(partner, fold)) %>%
  pivot_longer(c(fold, fold_smoking_adjusted), names_to = "model", values_to = "f") %>%
  mutate(model = recode(model, fold = "Before adjustment", fold_smoking_adjusted = "After smoking adjustment"),
         model = factor(model, c("Before adjustment", "After smoking adjustment")))
pc <- ggplot(sm, aes(f, partner, fill = model)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2fx", f)), position = position_dodge(0.75), hjust = -0.15, size = 2.2) +
  scale_fill_manual(values = c("grey65", "#0072B2")) +
  scale_x_continuous(expand = expansion(c(0, 0.18))) +
  labs(x = "Fold enrichment with KRAS G12C (1 = chance)", y = NULL, fill = NULL,
       title = "KRAS G12C partners and smoking (lung cancer)",
       caption = "Smoking explains part of G12C + STK11/KEAP1 and all of G12C + RBM10.") +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "italic"))

# d -- the 10 rarest combinations: patients seen vs expected
rc <- rd("fig4f_rare.csv") %>% filter(expected_both >= 20) %>% slice_max(times_rarer_than_expected, n = 10, with_ties = FALSE) %>%
  mutate(lab = paste0(sub("_fusion", " fusion", Gene_A), " + ", sub("_fusion", " fusion", Gene_B), "  (", short_ct(Cancer_Type), ")"),
         lab = reorder(lab, expected_both)) %>%
  pivot_longer(c(patients_with_both, expected_both), names_to = "what", values_to = "n") %>%
  mutate(what = recode(what, patients_with_both = "Seen", expected_both = "Expected by chance"),
         what = factor(what, c("Expected by chance", "Seen")))
pd_ <- ggplot(rc, aes(n, lab, fill = what)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = round(n)), position = position_dodge(0.75), hjust = -0.2, size = 2.2) +
  scale_fill_manual(values = c("grey70", "#762A83")) +
  scale_x_continuous(expand = expansion(c(0, 0.12))) +
  labs(x = "Patients with both alterations", y = NULL, fill = NULL, title = "The rarest combinations",
       caption = "Combinations seen far less often than chance predicts (expected >= 20 patients).") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5))

# One combined figure for the phase; the title states the main message
fig <- (pa | pb) / (pc | pd_) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = "Phase 4: exclusivity follows pathways, and pairs change with the cancer type")
save_fig(fig, "figure_phase4", 183, 175, OUT_DIR)
