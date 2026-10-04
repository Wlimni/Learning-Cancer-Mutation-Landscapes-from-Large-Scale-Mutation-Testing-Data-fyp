"""
Script Name: 08_allele_screen.py
Purpose : Variant-level co-mutation screen (Aim 2; Scharpf 2022, Cook 2021, Vaeyens 2023):
          hotspot genes are split into allele groups (e.g. KRAS G12C / G12D / G12V / other)
          and each allele is tested against every partner gene with the Phase 2
          conditional test (non-hypermutated patients, subtype strata). Then:
          (1) do alleles of the same gene differ in their partners? (Cochran's Q)
          (2) does KRAS G12C's partner pattern survive adjusting for the smoking
              mutational signature (per-patient C>A share, partner genes excluded)?
Inputs  : notebooks/02 (functions, via p2_variants.py), data/raw/data_mutations_extended.txt
Outputs : data/processed/allele_pairs.parquet, allele_heterogeneity.parquet,
          allele_smoking_adjustment.parquet (read by notebooks/04)
"""
import re, sys
import numpy as np, pandas as pd
from scipy.stats import chi2, false_discovery_control
sys.path.insert(0, "scripts/analyses")
import p2_variants as V

ns = V.ns
alts = ns["alterations"].copy()

def allele_group(gene, hgvs):
    """Allele label for a mutation in a hotspot gene; None = keep the plain gene name."""
    h = str(hgvs)
    if gene == "KRAS":
        for a in ["G12C", "G12D", "G12V", "G12A", "G12R", "G13D", "Q61H"]:
            if h == f"p.{a}": return f"KRAS_{a}"
        return "KRAS_other"
    if gene == "BRAF":
        return "BRAF_V600E" if h == "p.V600E" else "BRAF_nonV600"
    if gene == "EGFR":
        if h == "p.L858R": return "EGFR_L858R"
        m = re.match(r"p\.[A-Z](\d+)_[A-Z](\d+)(del|delins)", h)
        if m and 729 <= int(m.group(1)) <= 761: return "EGFR_ex19del"
        m = re.match(r"p\.[A-Z](\d+)_[A-Z](\d+)(ins|dup)", h) or re.match(r"p\.[A-Z](\d+)(dup)", h)
        if m and 762 <= int(m.group(1)) <= 823: return "EGFR_ex20ins"
        return "EGFR_other"
    if gene == "PIK3CA":
        if h in ("p.E545K", "p.E542K", "p.Q546K", "p.E545A", "p.E545G"): return "PIK3CA_helical"
        if h.startswith("p.H1047") or h.startswith("p.G1049"): return "PIK3CA_kinase"
        return "PIK3CA_other"
    if gene == "IDH1":
        return "IDH1_R132H" if h == "p.R132H" else "IDH1_other"
    return None

HOTSPOT = {"KRAS", "BRAF", "EGFR", "PIK3CA", "IDH1"}
mut = alts["Alteration_Type"] == "MUT"
hot = mut & alts["Hugo_Symbol"].isin(HOTSPOT)
alts.loc[hot, "Hugo_Symbol"] = [allele_group(g, h) for g, h in zip(alts.loc[hot, "Hugo_Symbol"], alts.loc[hot, "HGVSp_Short"])]
parent = {}
for lab in alts.loc[hot, "Hugo_Symbol"].unique():
    parent[lab] = lab.split("_")[0]
# Copy-number calls of hotspot genes keep the gene name (e.g. "EGFR" = EGFR amplification).

# Allele labels are tested wherever their parent gene is tested.
for key in ["panel_coverage", "cna_panel_coverage"]:
    cov = ns[key]
    extra = [cov[cov["Hugo_Symbol"] == p].assign(Hugo_Symbol=lab) for lab, p in parent.items()]
    ns[key] = pd.concat([cov] + extra, ignore_index=True)

# Smoking-signature proxy: share of C>A substitutions among a patient's SNVs, excluding the
# hotspot genes and partner genes tested here (so a G12C itself does not raise its own score).
maf = pd.read_csv("data/raw/data_mutations_extended.txt", sep="\t", low_memory=False,
                  usecols=["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2"])
snv = maf[(maf.Variant_Type == "SNP") & ~maf.Hugo_Symbol.isin(HOTSPOT | {"STK11", "KEAP1", "TP53", "SMARCA4", "NFE2L2"})]
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
ref, alt = snv.Reference_Allele.str.upper(), snv.Tumor_Seq_Allele2.str.upper()
pyr = ref.isin(["C", "T"])
is_ca = (ref.where(pyr, ref.map(comp)) == "C") & (alt.where(pyr, alt.map(comp)) == "A")
spec = is_ca.groupby(snv.Tumor_Sample_Barcode).agg(["size", "mean"])
spec = spec[spec["size"] >= 5]["mean"]

CTS = {"Non-Small Cell Lung Cancer": ["KRAS", "EGFR", "BRAF", "PIK3CA"], "Colorectal Cancer": ["KRAS", "BRAF", "PIK3CA"],
       "Pancreatic Cancer": ["KRAS"], "Breast Cancer": ["PIK3CA"], "Endometrial Cancer": ["KRAS", "PIK3CA"],
       "Glioma": ["IDH1"], "Melanoma": ["BRAF"]}


def run(ct, smoking=False):
    clin = ns["clinical"]
    if smoking:
        s = clin["SAMPLE_ID"].map(spec)
        in_ct = (clin["CANCER_TYPE"] == ct) & s.notna()        # tertiles within this cancer type
        tert = pd.Series("CA_unknown", index=clin.index)
        tert[in_ct] = pd.qcut(s[in_ct].rank(method="first"), 3, labels=["CA_low", "CA_mid", "CA_high"]).astype(str)
        clin = clin.assign(STRAT_SMOKE=clin["SUBTYPE_STRATUM"] + "|" + tert)
    saved = (ns["clinical"], ns["STRATUM_COLUMN"])
    ns["clinical"], ns["STRATUM_COLUMN"] = clin, ("STRAT_SMOKE" if smoking else "SUBTYPE_STRATUM")
    try:
        res = ns["test_cancer_type"](ct, alt_source=alts, require_cna_tested=True)
    finally:
        ns["clinical"], ns["STRATUM_COLUMN"] = saved
    if res.empty:
        return res
    res["q_value"] = ns["bh_fdr"](res["p_value"].to_numpy())
    res["Cancer_Type"] = ct
    return res


out = []
for ct, genes in CTS.items():
    res = run(ct)
    a_is = res.Gene_A.map(parent).isin(genes); b_is = res.Gene_B.map(parent).isin(genes)
    keep = (a_is ^ b_is)                                   # allele vs a non-hotspot partner
    r = res[keep].copy()
    r["allele"] = np.where(a_is[keep], r.Gene_A, r.Gene_B); r["partner"] = np.where(a_is[keep], r.Gene_B, r.Gene_A)
    r = r[~r.partner.map(parent).notna()]                  # partner must not be another allele label
    out.append(r)
    print(f"{ct:<28} allele-partner pairs tested: {len(r):>4}, significant: {(r.q_value < .05).sum()}", flush=True)
ap = pd.concat(out, ignore_index=True)
ap["gene"] = ap.allele.map(parent)
ap.to_parquet("data/processed/allele_pairs.parquet", index=False)

# (1) Do alleles of one gene differ in their partners?
ap["theta"] = np.log((ap.n_both_altered + 0.5) / (ap.expected_both + 0.5)); ap["var"] = 1 / (ap.n_both_altered + 0.5)
rows = []
for (ct, gene, partner), g in ap[~ap.allele.str.endswith("_other")].groupby(["Cancer_Type", "gene", "partner"]):
    if len(g) < 2: continue
    w = 1 / g["var"]; th = (w * g.theta).sum() / w.sum(); Q = (w * (g.theta - th) ** 2).sum()
    rows.append(dict(Cancer_Type=ct, gene=gene, partner=partner, n_alleles=len(g), p=chi2.sf(Q, len(g) - 1),
                     folds="; ".join(f"{a.split('_', 1)[1]} {np.exp(t):.2f}" for a, t in zip(g.allele, g.theta))))
het = pd.DataFrame(rows); het["q"] = false_discovery_control(het.p)
het.to_parquet("data/processed/allele_heterogeneity.parquet", index=False)
print(f"\n(gene, partner) combinations with >= 2 alleles tested: {len(het)}; alleles differ (q < 0.05): {(het.q < .05).sum()}")
print(het.sort_values("q").head(20)[["Cancer_Type", "gene", "partner", "folds", "q"]].to_string(index=False))

# (2) KRAS G12C in lung cancer, with and without the smoking-signature adjustment
lung = "Non-Small Cell Lung Cancer"
adj = run(lung, smoking=True)
sm = []
for allele in ["KRAS_G12C", "KRAS_G12D", "KRAS_G12V"]:
    for partner in ["STK11", "KEAP1", "SMARCA4", "TP53", "NTRK3", "ATM", "RBM10"]:
        sel = lambda d: d[((d.Gene_A == allele) & (d.Gene_B == partner)) | ((d.Gene_A == partner) & (d.Gene_B == allele))]
        a, b = sel(ap[ap.Cancer_Type == lung]), sel(adj)
        if len(a) and len(b):
            sm.append(dict(allele=allele, partner=partner, fold=a.fold_enrichment.iloc[0], q=a.q_value.iloc[0],
                           fold_smoking_adjusted=b.fold_enrichment.iloc[0], q_smoking_adjusted=b.q_value.iloc[0]))
sm = pd.DataFrame(sm); sm.to_parquet("data/processed/allele_smoking_adjustment.parquet", index=False)
print(sm.round(3).to_string(index=False))
