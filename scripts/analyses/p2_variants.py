"""
Script Name: p2_variants.py
Purpose : Run Phase 2 design variants with the notebook's OWN functions, to compare:
          (a) which patients: all (pooled) vs non-hypermutated only;
          (b) what counts as altered: mutation OR copy number (CNA-tested patients)
              vs mutation only (all mutation-tested patients);
          (c) gene filter: current (>= 3% of ALL patients) vs >= 3% of patients
              TESTED for the gene vs none (>= 5 patients only).
Inputs  : notebooks/02_phase2_comutation_matrix.ipynb (cells executed in a namespace)
Outputs : results/analyses/variants.parquet (one row per variant x cancer type)
"""
import hashlib, itertools, json, os, sys, time
import numpy as np, pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
_cwd = os.getcwd()
os.chdir(os.path.join(ROOT, "notebooks"))      # the notebook cells use paths relative to notebooks/
nb = json.load(open("02_phase2_comutation_matrix.ipynb"))
ns = {}
for cid in ["8455b24e", "57159198", "332f818c", "e3d5f851", "disc0v3rcd", "17d9c817", "056ff63e"]:
    src = "".join(next(c for c in nb["cells"] if c.get("id") == cid)["source"])
    exec(src, ns)
os.chdir(_cwd)
# The notebook now restricts `clinical` to its analysis group; the harness starts from all
# patients so that every variant (group="all" / "nonhyper") stays reproducible.
clinical = ns.get("clinical_all", ns["clinical"])
alterations = ns.get("alterations_all", ns["alterations"])
cov_mut, cov_cna = ns["panel_coverage"], ns["cna_panel_coverage"]
mask, Model, pair_test, bh = ns["_coverage_mask"], ns["ConditionalModel"], ns["pair_test"], ns["bh_fdr"]

# Direct (non-statistical) mutation-testability fixes: drop samples with no mutation
# data released at all, and genes covered only by intron/intergenic probes.
g = pd.read_csv(os.path.join(ROOT, "data/raw/genomic_information.txt"), sep="\t", low_memory=False,
                usecols=["Hugo_Symbol", "SEQ_ASSAY_ID", "Feature_Type", "includeInPanel"])
exon = g[(g.includeInPanel == True) & (g.Feature_Type == "exon")][["SEQ_ASSAY_ID", "Hugo_Symbol"]].drop_duplicates()
cov_mut_exon = cov_mut.merge(exon)
maf_samples = set(pd.read_csv(os.path.join(ROOT, "data/raw/data_mutations_extended.txt"), sep="\t", usecols=["Tumor_Sample_Barcode"],
                              low_memory=False).Tumor_Sample_Barcode)
panel_has_muts = clinical.assign(h=clinical.SAMPLE_ID.isin(maf_samples)).groupby("SEQ_ASSAY_ID").h.mean()
MUT_PANELS = set(panel_has_muts[panel_has_muts >= 0.5].index)     # e.g. drops VICC-02-XFV2, CHOP-FUSIP


def run(ct, mode="any", group="all", gene_rule="current", joint_cov=0.5, min_exp=5.0, subtype="no", center="all",
        genes_fixed=None):
    c = clinical[clinical.CANCER_TYPE == ct]
    if center == "MSK":
        c = c[c.CENTER == "MSK"]
    elif center == "nonMSK":
        c = c[c.CENTER != "MSK"]
    if mode == "any":
        c = c[c.CNA_TESTED]
        src = alterations
    else:                                                   # mutation only
        c = c[c.SEQ_ASSAY_ID.isin(MUT_PANELS)]
        src = alterations[alterations.Alteration_Type == "MUT"]
    if group == "nonhyper":
        c = c[c.HYPERMUTATION_STATUS == "not_hypermutated"]
    c = c.drop_duplicates("SAMPLE_ID").set_index("SAMPLE_ID")
    n = len(c)
    if n < 100:
        return None
    sub = src.loc[src.Sample_ID.isin(c.index), ["Sample_ID", "Hugo_Symbol"]].drop_duplicates()
    k = sub.groupby("Hugo_Symbol").Sample_ID.nunique()
    covm = cov_mut_exon if mode == "snv" else cov_mut
    tested_all = mask(covm, c, list(k.index)) & (mask(cov_cna, c, list(k.index)) if mode == "any" else True)
    n_tested = pd.Series(tested_all.sum(0), index=k.index)
    if gene_rule == "current":
        genes = k[k >= max(5, 0.03 * n)].index
    elif gene_rule == "of_tested":
        genes = k[(k >= 5) & (k >= 0.03 * n_tested)].index
    else:
        genes = k[k >= 5].index
    genes = sorted(genes) if genes_fixed is None else sorted(set(genes_fixed) & set(k.index))
    if len(genes) < 2:
        return None
    alt = (sub[sub.Hugo_Symbol.isin(genes)].assign(a=True).pivot(index="Sample_ID", columns="Hugo_Symbol", values="a")
           .reindex(index=c.index, columns=genes).fillna(False).astype(bool)).to_numpy()
    tested = mask(covm, c, genes) & (mask(cov_cna, c, genes) if mode == "any" else True)
    strata = c.HYPERMUTATION_STATUS.astype(str) if group == "all" else pd.Series("", index=c.index)
    if subtype == "yes":
        # detailed subtypes with >= 50 patients get their own gene rates; the rest share one
        det = c.CANCER_TYPE_DETAILED.astype(str)
        big = det.value_counts()
        det = det.where(det.map(big) >= 50, "other")
        strata = strata + "|" + det
    strata = strata.to_numpy() if (group == "all" or subtype == "yes") else None
    t0 = time.time()
    M = Model(alt, tested, strata)
    rows = []
    for i, j in itertools.combinations(range(len(genes)), 2):
        jt = tested[:, i] & tested[:, j]
        if jt.sum() / n < joint_cov:
            continue
        q = M.joint(i, j)[jt]
        if q.sum() < min_exp:
            continue
        both = int((alt[:, i] & alt[:, j] & jt).sum())
        p, d, e = pair_test(q, both)
        rows.append((genes[i], genes[j], int(jt.sum()), both, e, d, p))
    df = pd.DataFrame(rows, columns=["Gene_A", "Gene_B", "n_tested", "n_both", "expected", "direction", "p"])
    if len(df):
        df["q"] = bh(df.p.to_numpy())
    df.attrs.update(n_patients=n, n_genes=len(genes), secs=time.time() - t0)
    return df


if __name__ == "__main__":
    cts = sys.argv[1].split("|")
    variants = [dict(v.split("=") for v in a.split(",")) for a in sys.argv[2:]]
    out = []
    for ct in cts:
        for v in variants:
            kw = {k: (float(x) if k in ("joint_cov", "min_exp") else x) for k, x in v.items()}
            df = run(ct, **kw)
            if df is None:
                continue
            tag = ",".join(f"{k}={x}" for k, x in v.items())
            print(f"{ct:<28} {tag:<48} patients {df.attrs['n_patients']:>6,} genes {df.attrs['n_genes']:>4} "
                  f"pairs {len(df):>6,} sig {int((df.q < .05).sum()) if len(df) else 0:>5,} "
                  f"({(df.q < .05).mean() if len(df) else 0:.1%})  {df.attrs['secs']:.0f}s", flush=True)
            df.assign(Cancer_Type=ct, variant=tag).to_parquet(
                os.path.join(ROOT, "results/analyses/variants", f"var_{hashlib.md5((ct + tag).encode()).hexdigest()[:10]}.parquet"))
