# Script Name: phase1_data.R
# Purpose : Phase 1 figures -- the cohort and why panel data need care (notebook 01).
# Inputs  : results/figures/data/fig1*.csv (from 00_export_figure_data.py)
# Outputs : results/figures/phase1_data/figure_phase1.{pdf,png} -- one combined figure (4 panels)
source("scripts/figures/_theme.R")
OUT_DIR <- "phase1_data"

# a -- patient flow
fl <- rd("fig1a_flow.csv") %>% mutate(step = factor(step, rev(step)))
pa <- ggplot(fl, aes(n, step)) +
  geom_col(fill = "grey45", width = 0.65) +
  geom_text(aes(label = comma(n)), hjust = -0.1, size = 2.6) +
  scale_x_continuous(labels = label_comma(scale = 1e-3, suffix = "k"), expand = expansion(c(0, 0.45))) +
  labs(x = "Samples / patients", y = NULL, title = "Patient flow")

# b -- cohort composition (top 25 cancer types by patients)
cb <- rd("fig1b_cohort.csv")
top <- cb %>% count(CANCER_TYPE, wt = patients, sort = TRUE) %>% slice_head(n = 25)
centres <- cb %>% count(centre, wt = patients, sort = TRUE) %>% pull(centre)
centres <- c(setdiff(centres, "Other"), "Other")
pb <- cb %>% filter(CANCER_TYPE %in% top$CANCER_TYPE) %>%
  mutate(ct = factor(short_ct(CANCER_TYPE), rev(short_ct(top$CANCER_TYPE))), centre = factor(centre, centres)) %>%
  ggplot(aes(patients, ct, fill = centre)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "grey70")) +
  scale_x_continuous(labels = label_comma(), expand = expansion(c(0, 0.05))) +
  labs(x = "Patients", y = NULL, fill = "Centre",
       title = sprintf("Cohort: top 25 of %d cancer types", n_distinct(cb$CANCER_TYPE))) +
  theme(legend.position = c(0.75, 0.3))

# c -- what panels claim vs what they release, per data type (all panels with patients)
pn <- rd("fig1c_panels.csv") %>%
  mutate(across(c(MUT_TESTED:claims_SV), ~ coalesce(as.logical(.x), FALSE)))
cnt <- function(type, claim, rel) data.frame(
  type = type, status = c("Claimed and released", "Claimed, never released", "Released, not claimed"),
  n = c(sum(pn[[claim]] & pn[[rel]]), sum(pn[[claim]] & !pn[[rel]]), sum(!pn[[claim]] & pn[[rel]])))
pcd <- bind_rows(cnt("Mutation", "claims_MUT", "MUT_TESTED"), cnt("Copy number", "claims_CNA", "CNA_TESTED"),
                 cnt("Fusion", "claims_SV", "SV_TESTED")) %>%
  mutate(type = factor(type, rev(c("Mutation", "Copy number", "Fusion"))),
         status = factor(status, c("Claimed and released", "Claimed, never released", "Released, not claimed")))
pc <- ggplot(pcd, aes(n, type, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5), size = 2.6, colour = "white") +
  scale_fill_manual(values = c("Claimed and released" = "grey40", "Claimed, never released" = "#D55E00",
                               "Released, not claimed" = "#56B4E9")) +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  labs(x = "Number of panels", y = NULL, fill = NULL, title = "Panel metadata over-claims what was tested",
       caption = "Orange: the panel says it tests this, but no data were ever released.\nSo testing is read from the data, not the metadata.") +
  guides(fill = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom")

# d -- alterations per patient (analysed genes), 20 largest cancer types
ad <- rd("fig1d_alterations_per_patient.csv")
top20 <- ad %>% count(CANCER_TYPE, sort = TRUE) %>% slice_head(n = 20) %>% pull(CANCER_TYPE)
ad20 <- ad %>% filter(CANCER_TYPE %in% top20) %>% mutate(ct = short_ct(CANCER_TYPE))
ord <- ad20 %>% group_by(ct) %>% summarise(m = mean(n_altered)) %>% arrange(m) %>% pull(ct)
pd_ <- ggplot(ad20, aes(n_altered + 1, factor(ct, ord))) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, fill = "grey90") +
  stat_summary(fun = mean, geom = "point", size = 1, colour = "#D55E00") +
  scale_x_log10(breaks = c(1, 2, 3, 5, 11, 21), labels = c(0, 1, 2, 4, 10, 20)) +
  coord_cartesian(xlim = c(1, 25)) +
  labs(y = NULL, x = "Analysed genes altered per patient", title = "Alteration load per patient",
       caption = "Box: median and IQR; orange: mean (20 largest cancer types).\nA pair test that ignores this load mistakes it for interaction.")

# One combined figure for the phase; the title states the main message
fig <- ((pa / pb + plot_layout(heights = c(1, 3.4))) | (pc / pd_ + plot_layout(heights = c(0.9, 1.6)))) +
  plot_annotation(tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 11)), title = "Phase 1: 112,779 non-hypermutated patients in 51 cancer types -- and panel data need care")
save_fig(fig, "figure_phase1", 183, 205, OUT_DIR)
