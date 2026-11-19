# Comparing the two copy-number filters

The established-driver filter reduces the number of Phase 2 candidates and
improves recovery in simulations with planted links. It also drops some
pairs that replicate, so the choice depends on how much coverage we want
to trade for a narrower candidate list.

| Measure | Current SNV+CNV | SNV+CNV with established CNV drivers |
|---|---:|---:|
| Tested cancer type and gene pair combinations | 6,809 | 5,641 |
| Interaction candidates | 1,569 | 1,232 |
| Candidates replicated outside MSK | 472 | 433 |
| Replication among all MSK hits eligible for replication | 533/814 (65.5%) | 482/704 (68.5%) |
| Published pairs recovered (58 total) | 44 | 43 |
| Precision / recall in planted-link simulation | 96.2% / 46.3% | 99.1% / 49.1% |
| TCGA same-direction candidates (control pairs) | 88% (69%) | 87% (64%) |

The current filter keeps deep calls whose direction is compatible with
known gene roles, including calls for genes without a reviewed role.
The narrower filter also requires an established driver tag. The Phase 1
[lookup table](../data/processed/cnv_driver_list.parquet) marks **445**
gene–direction entries as established and 468 as uncertain. The tag uses
evidence from at least one of COSMIC Cancer Gene Census v99, Vogelstein
2013, or CIViC, with a direction check. It is a gene-level rule. It does not prove that a copy-number
event drives a particular tumour or resolve which gene is causal in an amplicon.

The candidate lists overlap on **1,151** cancer type and gene pair
combinations, with the same direction in both analyses. The comparison
matches on cancer type and both genes in the two saved Phase 2 pair tables.

| Candidate group | Pairs | Eligible for replication | Replicated |
|---|---:|---:|---:|
| Shared, evaluated with the current filter | 1,151 | 646 | 449 |
| Shared, evaluated with the established-driver filter | 1,151 | 620 | 429 |
| Current filter only | 418 | 70 | 23 |
| Established-driver filter only | 81 | 13 | 4 |

For example, the current filter alone calls **FGF19–TP53**, **FGF4–TP53**
and **PIK3CA–WHSC1L1** exclusive in breast cancer. They did not replicate
in the saved results and are examples to review for passenger effects.
It also calls **BCL2L1–KRAS** exclusive in colorectal cancer, and that pair
did replicate. The narrower filter uniquely recovers **CBFB–PIK3CA**
co-occurrence in breast cancer, which also replicates.

Only 23 of the 418 current-only candidates had a recorded replication
(5.5% of all 418). Among the 70 eligible for that check, 23 replicated
(32.9%). Ineligible candidates do not contribute evidence either way.

I would use the established-driver list as the primary list if the priority
is a narrower set of candidate driver interactions, and keep the broader
list as a sensitivity analysis. The current main table is still SNV+CNV.
Before switching, compare the 23 current-only replicated pairs against
the evidence for their copy-number drivers.

Sources: [current pair table](../data/processed/comutation_pairs.parquet),
[established-driver pair table](../data/processed/comutation_pairs_snv_cnv_drivers.parquet),
and [Phase 2b validation](../notebooks/02_pair_analysis/02b_phase2_validation.ipynb).
