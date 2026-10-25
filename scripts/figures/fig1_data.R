# Script Name: fig1_data.R
# Purpose : Figure 1 -- the cohort and why the data need care.
#           a patient flow; b patients per cancer type and centre (GENIE 2017 Fig 2B);
#           c what panels claim vs what they release (GENIE 2017 Fig 2A/3 idea);
#           d alterations per patient by cancer type (GENIE 2017 Fig 2C idea).
# Inputs  : results/figures/data/fig1*.csv (from 00_export_figure_data.py)
# Outputs : results/figures/fig1_data.{pdf,png}
source("scripts/figures/_theme.R")

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

# c -- panel metadata vs released data (30 largest panels)
pn <- rd("fig1c_panels.csv") %>%
  mutate(across(c(MUT_TESTED:claims_SV), ~ coalesce(as.logical(.x), FALSE)))
status <- function(claim, rel) case_when(claim & rel ~ "Claimed and released", claim & !rel ~ "Claimed, never released",
                                         !claim & rel ~ "Released, not claimed", TRUE ~ "Not tested")
pc_long <- pn %>% slice_max(patients, n = 30) %>%
  transmute(panel = factor(panel, rev(panel)), patients,
            Mutation = status(claims_MUT, MUT_TESTED), `Copy number` = status(claims_CNA, CNA_TESTED),
            Fusion = status(claims_SV, SV_TESTED)) %>%
  pivot_longer(c(Mutation, `Copy number`, Fusion), names_to = "type", values_to = "status") %>%
  mutate(type = factor(type, c("Mutation", "Copy number", "Fusion")))
n_bad <- c(cna = sum(pn$claims_CNA & !pn$CNA_TESTED), sv = sum(pn$claims_SV & !pn$SV_TESTED),
           mut = sum(pn$claims_MUT & !pn$MUT_TESTED))
pc <- ggplot(pc_long, aes(type, panel, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Claimed and released" = "#4D4D4D", "Claimed, never released" = "#D55E00",
                               "Released, not claimed" = "#56B4E9", "Not tested" = "grey92")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Panel metadata vs released data",
       caption = sprintf("30 largest of %d panels. Across all panels, claimed but never released:\ncopy number %d, fusions %d, mutations %d",
                         nrow(pn), n_bad["cna"], n_bad["sv"], n_bad["mut"])) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 5.5),
        axis.text.x.top = element_text(angle = 35, hjust = 0, vjust = 0),
        legend.position = "right")

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

fig <- (pa / pb + plot_layout(heights = c(1, 3.4))) | (pc / pd_ + plot_layout(heights = c(1.5, 1)))
fig <- fig + plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = list(c("a", "b", "c", "d")))
save_fig(fig, "fig1_data", 183, 200)
