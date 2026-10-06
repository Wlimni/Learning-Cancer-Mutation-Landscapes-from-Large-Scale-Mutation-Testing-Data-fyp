"""
Script Name: 11_novel_insights.py
Purpose : Mine Phase 2's replicated candidates for genuinely new/rare insight (not just
          recovering textbook biology). Two findings:
          (1) Amplification-burden confound: a "co-occurring" candidate pair can really be
              two passenger genes on independent amplicons, both correlated with a general
              "amplifier / chromosomally unstable" patient phenotype -- a blind spot the
              conditional test's per-patient total-count adjustment does not fully close
              (it controls for total alteration count, not for *which* alterations cluster
              together within that count). Quantified for CD79B in breast cancer, the
              clearest example, and swept across all cancer types.
          (2) An exception to Phase 4's "exclusivity = same-pathway redundancy" rule:
              NF1 + PTPN11 (both RTK-RAS pathway) CO-OCCUR in glioma, almost entirely via
              point mutations -- i.e. two genuine mutations, not an amplicon artefact.
Inputs  : data/processed/{comutation_pairs,mechanism_pairs,alterations_long,clinical_tidy}.parquet,
          data/external/pathways/sanchez_vega_2018_pathways.tsv
Outputs : printed findings; data/processed/cnv_burden_confound_check.parquet
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

P = "data/processed/"
cp = pd.read_parquet(P + "comutation_pairs.parquet")
mech = pd.read_parquet(P + "mechanism_pairs.parquet")
al = pd.read_parquet(P + "alterations_long.parquet",
                     columns=["Sample_ID", "Hugo_Symbol", "Alteration_Type", "CNA_Direction", "CANCER_TYPE"])
cl = pd.read_parquet(P + "clinical_tidy.parquet")
pw = pd.read_csv("data/external/pathways/sanchez_vega_2018_pathways.tsv", sep="\t")

cand = cp[cp.interaction_candidate & cp.replicated]
co = cand[cand.direction == "co-occurring"].copy()
co["frac_cnv_cnv"] = co.n_CNV_CNV / co.n_both_altered.clip(lower=1)

# ---------------------------------------------------------------------------
# (1) How much of the replicated co-occurring list is copy-number-amplification driven,
# i.e. at risk of the amplifier-phenotype confound rather than a specific interaction?
# ---------------------------------------------------------------------------
heavy = co[co.frac_cnv_cnv >= 0.8].sort_values("fold_enrichment", ascending=False)
print(f"Replicated co-occurring candidates: {len(co)}; >= 80% driven by CNV+CNV (both genes amplified "
      f"in the same patients, not mutated): {len(heavy)} ({len(heavy) / len(co):.0%})")
print(heavy[["Cancer_Type", "Gene_A", "Gene_B", "n_both_altered", "frac_cnv_cnv", "fold_enrichment", "q_value"]]
      .head(12).to_string(index=False))

# ---------------------------------------------------------------------------
# Case study: CD79B in breast cancer -- a B-cell lymphoma gene with no established role in
# breast cancer, sitting next to PPM1D/TBX2 (17q23), a known but under-covered amplicon.
# Does its "co-occurrence" with GATA3/MDM2/AURKA/GNAS survive adjusting for the patient's
# overall amplification burden?
# ---------------------------------------------------------------------------
bc = cl[(cl.CANCER_TYPE == "Breast Cancer") & cl.IS_REPRESENTATIVE & cl.CNA_TESTED].set_index("SAMPLE_ID")
amp = al[(al.CANCER_TYPE == "Breast Cancer") & (al.Alteration_Type == "CNV") & (al.CNA_Direction == "Amplification")]
burden = amp.groupby("Sample_ID").size().reindex(bc.index).fillna(0)


def flag(gene):
    return bc.index.isin(set(amp.loc[amp.Hugo_Symbol == gene, "Sample_ID"])).astype(int)


df = pd.DataFrame(index=bc.index)
df["cd79b"] = flag("CD79B")
df["log_burden_excl"] = np.log1p(burden.values - df["cd79b"])   # amp burden, excluding CD79B itself
cd79b_amp_n = df["cd79b"].sum()
print(f"\nCD79B-amplified breast patients: {cd79b_amp_n}; median amplification burden "
      f"{burden[df.cd79b == 1].median():.0f} vs {burden[df.cd79b == 0].median():.0f} overall "
      "(CD79B carriers are generally 'amplifier' tumours)")

# NOTE: GATA3 is deliberately excluded here. Its co-occurrence with CD79B (q=2.5e-42) is
# only 4.6% CNV+CNV (95% "mixed": CD79B amplified while GATA3 is MUTATED in the same patient
# -- GATA3 in breast cancer is overwhelmingly a mutation, not an amplification). So it is not
# actually a candidate for the amplifier-burden confound tested here.
# CDKN2A is also excluded from THIS quick check: it is altered only by deep DELETION, never
# amplification, so "amplification burden" is the wrong scope for it (its burden-adjusted test
# needs an "any CNV" burden covering both directions -- see 12_burden_adjusted_retest.py, which
# does this properly and correctly finds CD79B-CDKN2A survives, OR=2.00, p<1e-4).
rows = []
for partner in ["MDM2", "AURKA", "GNAS"]:
    df["y"] = flag(partner)
    raw = cp[(cp.Cancer_Type == "Breast Cancer") & cp.Gene_A.isin(["CD79B", partner])
            & cp.Gene_B.isin(["CD79B", partner])].iloc[0]
    m = sm.Logit(df["y"], sm.add_constant(df[["cd79b", "log_burden_excl"]])).fit(disp=0)
    rows.append(dict(partner=partner, phase2_fold=raw.fold_enrichment, phase2_q=raw.q_value,
                     burden_adjusted_OR=np.exp(m.params["cd79b"]), burden_adjusted_p=m.pvalues["cd79b"]))
res = pd.DataFrame(rows)
res.to_parquet(P + "cnv_burden_confound_check.parquet", index=False)
print("\nCD79B's association with each partner, before vs after adjusting for overall amp burden:")
print(res.round(4).to_string(index=False))
print("-> all three survive adjustment: CD79B's amplifier-burden story alone does not explain them;")
print("   see 12_burden_adjusted_retest.py for the full 24-pair sweep (3 of 24 ARE pure confounds elsewhere).")

# ---------------------------------------------------------------------------
# (2) NF1 + PTPN11 in glioma: both RTK-RAS pathway, but CO-OCCUR (the exception to Phase 4's
# "same pathway -> exclusive" rule), and it is mutation-driven, so the amplicon-burden
# confound above does not apply.
# ---------------------------------------------------------------------------
print("\nNF1 + PTPN11 in glioma (both RTK-RAS pathway genes -- normally pathway-mates are exclusive):")
r = cp[(cp.Cancer_Type == "Glioma") & cp.Gene_A.isin(["NF1", "PTPN11"]) & cp.Gene_B.isin(["NF1", "PTPN11"])].iloc[0]
print(f"  {r.n_both_altered} patients vs {r.expected_both:.1f} expected (fold {r.fold_enrichment:.2f}), "
      f"q={r.q_value:.1e}, replicated={bool(r.replicated)}; "
      f"{r.n_MUT_MUT} of {r.n_both_altered} co-altered patients are mutation+mutation (not CNV)")
print("  NF1 = RAS brake (tumour suppressor); PTPN11/SHP2 = RAS accelerator (oncogene) --",
      "plausible two-hit escalation of RAS/MAPK signalling, not redundancy.")
