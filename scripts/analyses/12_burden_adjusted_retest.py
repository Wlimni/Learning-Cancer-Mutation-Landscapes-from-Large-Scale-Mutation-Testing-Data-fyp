"""
Script Name: 12_burden_adjusted_retest.py
Purpose : Systematic version of the CD79B check (11_novel_insights.py): for every replicated
          co-occurring candidate that is >= 80% driven by CNV+CNV (both genes amplified/deleted
          together, not mutated), re-test the pair after adjusting for each patient's overall
          copy-number "burden" in that cancer type (how many OTHER genes are altered by copy
          number) -- a proxy for a generally unstable "amplifier" tumour, which the Phase 2
          conditional test does not separate from a specific two-gene interaction.
Method  : logistic regression, gene B's altered status ~ gene A's altered status + log(1 +
          burden excluding A and B), on the same patients and gene set Phase 2 used (CNA-tested,
          non-hypermutated, that cancer type). Burden = count of OTHER genes (from that cancer
          type's full candidate gene list, not just the 222-gene screened set) altered by deep
          copy number in the same patient -- the broadest, most conservative proxy available.
Inputs  : data/processed/{comutation_pairs,mechanism_pairs,clinical_tidy,alterations_long}.parquet
Outputs : data/processed/burden_adjusted_retest.parquet; printed summary
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

P = "data/processed/"
cp = pd.read_parquet(P + "comutation_pairs.parquet")
cl = pd.read_parquet(P + "clinical_tidy.parquet")
al = pd.read_parquet(P + "alterations_long.parquet",
                     columns=["Sample_ID", "Hugo_Symbol", "Alteration_Type", "CANCER_TYPE"])

cand = cp[cp.interaction_candidate & cp.replicated & (cp.direction == "co-occurring")].copy()
cand["frac_cnv_cnv"] = cand.n_CNV_CNV / cand.n_both_altered.clip(lower=1)
targets = cand[cand.frac_cnv_cnv >= 0.8].copy()
print(f"{len(targets)} replicated co-occurring candidates are >= 80% CNV+CNV driven -- re-testing each")

rows = []
for ct in targets.Cancer_Type.unique():
    base = cl[(cl.CANCER_TYPE == ct) & cl.IS_REPRESENTATIVE & cl.CNA_TESTED
             & (cl.HYPERMUTATION_STATUS == "not_hypermutated")].set_index("SAMPLE_ID")
    cnv = al[(al.CANCER_TYPE == ct) & (al.Alteration_Type == "CNV")]
    cnv = cnv[cnv.Sample_ID.isin(base.index)]
    # Burden = how many distinct genes are copy-number-altered in this patient, within this
    # cancer type's full altered-gene pool (not restricted to the 222 genes Phase 2 screened).
    per_patient_genes = cnv.groupby("Sample_ID")["Hugo_Symbol"].apply(set)

    def altered(gene):
        s = set(cnv.loc[cnv.Hugo_Symbol == gene, "Sample_ID"])
        return base.index.isin(s).astype(int)

    for _, r in targets[targets.Cancer_Type == ct].iterrows():
        a, b = r.Gene_A, r.Gene_B
        xa, xb = altered(a), altered(b)
        burden = per_patient_genes.reindex(base.index).apply(lambda s: len(s - {a, b}) if isinstance(s, set) else 0)
        d = pd.DataFrame({"a": xa, "b": xb, "log_burden": np.log1p(burden.to_numpy())}, index=base.index)
        try:
            m = sm.Logit(d["b"], sm.add_constant(d[["a", "log_burden"]])).fit(disp=0)
            or_adj, p_adj = np.exp(m.params["a"]), m.pvalues["a"]
        except Exception as e:
            or_adj, p_adj = np.nan, np.nan
        rows.append(dict(Cancer_Type=ct, Gene_A=a, Gene_B=b, n_both_altered=r.n_both_altered,
                         phase2_fold=r.fold_enrichment, phase2_q=r.q_value,
                         burden_adjusted_OR=or_adj, burden_adjusted_p=p_adj,
                         median_burden_in_carriers=burden[xa.astype(bool) | xb.astype(bool)].median(),
                         median_burden_overall=burden.median()))
        print(f"  {ct:<28} {a:<8} {b:<8} raw fold {r.fold_enrichment:.2f} (q={r.q_value:.1e})  ->  "
              f"adjusted OR {or_adj:.2f} (p={p_adj:.3f})", flush=True)

out = pd.DataFrame(rows)
out["verdict"] = np.select(
    [out.burden_adjusted_p >= 0.05, (out.burden_adjusted_OR < 1) & (out.burden_adjusted_p < 0.05)],
    ["confound: not significant after adjustment", "confound: REVERSES direction after adjustment"],
    "survives: still co-occurring after adjustment")
out.to_parquet(P + "burden_adjusted_retest.parquet", index=False)
print("\n" + out.verdict.value_counts().to_string())
print(f"\nSurvives (likely real): {(out.verdict.str.startswith('survives')).sum()} of {len(out)}")
print(out.sort_values("burden_adjusted_p")[["Cancer_Type", "Gene_A", "Gene_B", "phase2_fold", "burden_adjusted_OR",
                                            "burden_adjusted_p", "verdict"]].to_string(index=False))
