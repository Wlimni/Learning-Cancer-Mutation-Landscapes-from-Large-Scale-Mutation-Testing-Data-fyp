# Script Name: phase1_data.R
# Purpose : Phase 1 figure -- what Phase 1 hands to Phase 2 (notebook 01, Section 9):
#           a patient flow, each filter step and what it removes; b TMB groups (hypermutated tumours are
#           analysed separately); c filters applied to each alteration type; d the events Phase 2 receives.
# Inputs  : results/figures/data/fig1a_flow.csv, fig1b_tmb_groups.csv, fig1c_filters.csv, fig1d_composition*.csv
#           (from 00_export_figure_data.py, which reads cohort.parquet and events.parquet)
# Outputs : results/figures/phase1_data/figure_phase1.{pdf,png} -- one combined figure (4 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase1_data"
COL_TMB <- c(non_hypermutated = "#0072B2", hypermutated = "#E69F00", ultramutated = "#D55E00", unknown = "grey70")
LAB_TMB <- c(non_hypermutated = "Not hypermutated: Phase 2", hypermutated = "Hypermutated (TMB 10-100): Phase 2c",
             ultramutated = "Ultramutated (TMB >= 100)", unknown = "TMB unknown (panel < 1 Mb)")
COL_EVENT <- c(Mutation = "#009E73", Amplification = "#D55E00", `Deep deletion` = "#0072B2", Fusion = "#CC79A7")

# a -- patient flow: every filter between GENIE and Phase 2, with what it removes
fl <- rd("fig1a_flow.csv") %>%
  mutate(lab = ifelse(is.na(removed) | removed == "", step, paste0(step, "\n- ", removed)),
         lab = factor(lab, rev(lab)), last = row_number() == n())
pa <- ggplot(fl, aes(n, lab, fill = last)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = comma(n)), hjust = -0.1, size = 2.6) +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "#0072B2"), guide = "none") +
  scale_x_continuous(labels = label_comma(scale = 1e-3, suffix = "k"), expand = expansion(c(0, 0.45))) +
  labs(x = "Patients (samples in the first row)", y = NULL, title = "From GENIE to the Phase 2 cohort") +
  theme(axis.text.y = element_text(size = 6.5, lineheight = 0.9))

# b -- TMB groups among patients tested for mutation and copy number, 15 largest cancer types
tg <- rd("fig1b_tmb_groups.csv") %>% group_by(CANCER_TYPE) %>%
  mutate(n = sum(patients), share = patients / n, hyper = sum(share[TMB_GROUP %in% c("hypermutated", "ultramutated")])) %>%
  ungroup() %>%
  mutate(ct = reorder(sprintf("%s (%s)", short_ct(CANCER_TYPE), comma(n)), hyper),
         TMB_GROUP = factor(TMB_GROUP, rev(names(COL_TMB))))
pb <- ggplot(tg, aes(share, ct, fill = TMB_GROUP)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = COL_TMB, labels = LAB_TMB, breaks = names(COL_TMB)) +
  scale_x_continuous(labels = percent, expand = c(0, 0)) +
  labs(x = "Share of patients", y = NULL, fill = NULL, title = "Hypermutated tumours are set aside") +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5))

# c -- filters applied to each alteration type: share of the first step's calls kept, counts as labels
fc <- rd("fig1c_filters.csv") %>%
  mutate(type = factor(type, c("Mutations", "Copy number", "Fusions")), step = factor(step, rev(unique(step))),
         kept = !duplicated(type, fromLast = TRUE)) %>%
  group_by(type) %>% mutate(share = n / first(n)) %>% ungroup()
pc <- ggplot(fc, aes(share, step, fill = ifelse(kept, as.character(type), "before"))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label_number(scale_cut = cut_short_scale(), accuracy = 0.1)(n)), hjust = -0.1, size = 2.4) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(before = "grey70", Mutations = COL_EVENT[["Mutation"]],
                               `Copy number` = COL_EVENT[["Amplification"]], Fusions = COL_EVENT[["Fusion"]]), guide = "none") +
  scale_x_continuous(labels = percent, breaks = c(0, 0.5, 1), expand = expansion(c(0, 0.35))) +
  labs(x = "Share of calls kept (labels: number of calls, all samples)", y = NULL, title = "Only likely drivers are kept",
       caption = "Coloured bar: what is kept. Mutations: AlphaMissense / PolyPhen-SIFT or hotspot evidence.\nCopy number: deep calls in the direction of the gene's role (OncoKB).") +
  theme(strip.text.y = element_text(angle = 0, face = "bold"), axis.text.y = element_text(size = 6.5))

# d -- what Phase 2 receives: driver events per patient by type, 15 largest cancer types
tot <- rd("fig1d_composition_total.csv")
cm <- rd("fig1d_composition.csv") %>% mutate(per_patient = events / patients) %>%
  group_by(ct) %>% mutate(total = sum(per_patient)) %>% ungroup() %>%
  mutate(ct = reorder(short_ct(ct), total), event = factor(event, rev(names(COL_EVENT))))
pd_ <- ggplot(cm, aes(per_patient, ct, fill = event)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = COL_EVENT, breaks = names(COL_EVENT)) +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  labs(x = "Driver events per patient", y = NULL, fill = NULL, title = "What Phase 2 receives",
       caption = sprintf("%s events for %s patients:\n%s.", comma(sum(tot$events)), comma(tot$patients[1]),
                         paste(sprintf("%s %.0f%%", tolower(tot$event), 100 * tot$events / sum(tot$events)), collapse = ", "))) +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5))

# One combined figure for the phase; the title states what is handed over
fig <- (pa | pb) / (pc | pd_) + plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)),
                  title = sprintf("Phase 1: %s patients and their driver events, ready for Phase 2", comma(fl$n[nrow(fl)])))
save_fig(fig, "figure_phase1", 183, 215, OUT_DIR)
