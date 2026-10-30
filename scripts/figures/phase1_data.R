# Script Name: phase1_data.R
# Purpose : Phase 1 figure -- what Phase 1 hands to Phase 2 (notebook 01, Section 9):
#           f the most altered genes; a patient groups by hypermutation; b the groups in the largest cancer types; c alteration counts of
#           normal vs hypermutated patients; d calls before vs after filtering; e altered genes per patient by cancer type.
# Inputs  : results/figures/data/fig1a_groups.csv, fig1b_tmb_groups.csv, fig1c_events_by_group.csv,
#           fig1d_before_after.csv, fig1e_alterations_per_patient.csv, fig1f_top_genes.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase1_data/figure_phase1.{pdf,png} -- one combined figure (6 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase1_data"
COL_TMB <- c(non_hypermutated = "#0072B2", hypermutated = "#E69F00", ultramutated = "#D55E00", unknown = "grey70")
LAB_TMB <- c(non_hypermutated = "Normal: Phase 2", hypermutated = "Hypermutated (TMB 10-100): Phase 2c",
             ultramutated = "Ultramutated (TMB >= 100)", unknown = "TMB unknown (panel < 1 Mb)")
COL_EVENT <- c(Mutation = "#009E73", Amplification = "#D55E00", `Deep deletion` = "#0072B2", Fusion = "#CC79A7")
cnt <- label_comma()

# a -- who the 138,206 patients tested for mutation and copy number are
ga <- rd("fig1a_groups.csv") %>%
  mutate(col = c("non_hypermutated", "hypermutated", "ultramutated", "unknown", "excluded"),
         group = factor(group, rev(group)))
pa <- ggplot(ga, aes(patients, group, fill = col)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = cnt(patients)), hjust = -0.1, size = 2.8) +
  scale_fill_manual(values = c(COL_TMB, excluded = "grey70"), guide = "none") +
  scale_x_continuous(labels = label_comma(scale = 1e-3, suffix = "k"), expand = expansion(c(0, 0.4))) +
  labs(x = "Patients (one sample each)", y = NULL, title = "Patient groups",
       caption = sprintf("%s patients tested for mutation and copy number.", cnt(sum(ga$patients))))

# b -- the same groups in the 15 largest cancer types
tg <- rd("fig1b_tmb_groups.csv") %>% group_by(CANCER_TYPE) %>%
  mutate(n = sum(patients), share = patients / n,
         hyper = sum(share[TMB_GROUP %in% c("hypermutated", "ultramutated")])) %>%
  ungroup() %>%
  mutate(ct = reorder(sprintf("%s (%s)", short_ct(CANCER_TYPE), cnt(n)), hyper),
         TMB_GROUP = factor(TMB_GROUP, rev(names(COL_TMB))))
pb <- ggplot(tg, aes(share, ct, fill = TMB_GROUP)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = COL_TMB, labels = LAB_TMB, breaks = names(COL_TMB)) +
  scale_x_continuous(labels = percent, expand = c(0, 0)) +
  labs(x = "Share of patients", y = NULL, fill = NULL, title = "Hypermutation depends on cancer type") +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5))

# c -- alteration composition by group: share of each type, with the load per patient written at the bar end
gc <- rd("fig1c_events_by_group.csv") %>%
  mutate(grp = factor(sprintf("%s\n(%s patients)", c(non_hypermutated = "Normal", hypermutated = "Hypermutated",
                                                      ultramutated = "Ultramutated")[TMB_GROUP], cnt(patients)),
                      sprintf("%s\n(%s patients)", c("Ultramutated", "Hypermutated", "Normal"), cnt(c(477, 15501, 113042)))),
         event = factor(event, rev(names(COL_EVENT)))) %>%
  group_by(grp) %>% mutate(share = events / sum(events), load = sum(events) / patients) %>% ungroup()
pc <- ggplot(gc, aes(share, grp, fill = event)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(share > 0.06, percent(share, accuracy = 1), "")), position = position_stack(vjust = 0.5),
            size = 2.4, colour = "white") +
  geom_text(data = distinct(gc, grp, load), aes(1.02, grp, label = sprintf("%.1f per\npatient", load), fill = NULL),
            hjust = 0, size = 2.5, lineheight = 0.9) +
  scale_fill_manual(values = COL_EVENT, breaks = names(COL_EVENT)) +
  scale_x_continuous(labels = percent, breaks = c(0, 0.5, 1), expand = expansion(c(0, 0.3))) +
  labs(x = "Share of alterations", y = NULL, fill = NULL, title = "Hypermutated tumours are mostly mutations") +
  guides(fill = guide_legend(nrow = 1)) + coord_cartesian(clip = "off") + theme(legend.position = "bottom")

# d -- raw calls vs calls kept (all GENIE samples)
bd <- rd("fig1d_before_after.csv") %>%
  pivot_longer(c(raw, kept), names_to = "stage", values_to = "n") %>%
  mutate(event = factor(event, names(COL_EVENT)), stage = factor(stage, c("kept", "raw")),
         fill = ifelse(stage == "raw", "raw", as.character(event)),
         lab = ifelse(stage == "raw", cnt(n), sprintf("%s (%.0f%%)", cnt(n), 100 * n / n[stage == "raw"])),
         .by = event)
pd_ <- ggplot(bd, aes(n, stage, fill = fill)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = lab), hjust = -0.05, size = 2.5) +
  facet_wrap(~event, ncol = 1, scales = "free_x", strip.position = "right") +
  scale_fill_manual(values = c(COL_EVENT, raw = "grey70"), guide = "none") +
  scale_y_discrete(labels = c(raw = "Raw", kept = "Kept")) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()), expand = expansion(c(0, 0.45))) +
  labs(x = "Calls, all GENIE samples", y = NULL, title = "Filtering keeps likely drivers") +
  theme(strip.text.y = element_text(angle = 0, face = "bold"))

# e -- altered genes per patient in the Phase 2 cohort, box = quartiles
be <- rd("fig1e_alterations_per_patient.csv") %>%
  mutate(all = ct == "All cancer types",
         lab = sprintf("%s (%s)", short_ct(ct), cnt(n)))
be$lab <- factor(be$lab, rev(c(be$lab[be$all], be$lab[!be$all][order(-be$median[!be$all], -be$mean[!be$all])])))
pe <- ggplot(be, aes(y = lab)) +
  geom_boxplot(aes(xmin = lo, xlower = q1, xmiddle = median, xupper = q3, xmax = hi, fill = all),
               stat = "identity", width = 0.6, linewidth = 0.35) +
  geom_point(aes(x = mean), shape = 23, size = 2, fill = "white", colour = "black") +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "grey80"), guide = "none") +
  labs(x = "Altered genes per patient", y = NULL, title = "Altered genes per patient",
       caption = "Box = quartiles with median, diamond = mean, whiskers = 1.5 x IQR.\nMedian 1 (B-cell) vs 6 (colorectal): one shared cut-off would\nmisjudge both, so Phase 2 tests within each cancer type.") +
  theme(axis.text.y = element_text(size = 6.5))

# f -- most altered genes in the Phase 2 cohort; share = patients altered / patients whose panel tested the gene
bf <- rd("fig1f_top_genes.csv") %>% mutate(share = patients / tested)
ord <- bf %>% group_by(Hugo_Symbol) %>% summarise(t = sum(share), .groups = "drop") %>% arrange(t)
bf <- bf %>% mutate(Hugo_Symbol = factor(Hugo_Symbol, ord$Hugo_Symbol), event = factor(event, rev(names(COL_EVENT))))
pf <- ggplot(bf, aes(share, Hugo_Symbol, fill = event)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = COL_EVENT, breaks = names(COL_EVENT)) +
  scale_x_continuous(labels = percent, expand = expansion(c(0, 0.05))) +
  labs(x = "Patients altered (of those tested for the gene)", y = NULL, fill = NULL, title = "Most altered genes",
       caption = "Genes tested in at least half the patients. A patient with\ntwo types counts in both. Fusions are not shown here.") +
  guides(fill = guide_legend(nrow = 1)) + theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5, face = "italic"))

fig <- (pa | pb) / (pc | pd_) / (pe | pf) + plot_layout(heights = c(1, 1.1, 1.2)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = sprintf("Phase 1: %s normal patients and their alterations, ready for Phase 2", cnt(ga$patients[1])))
save_fig(fig, "figure_phase1", 183, 262, OUT_DIR)
