"""
Script Name: 23_alteration_set_validation.py
Purpose : Precision and recall of the Phase 2 test for each alteration set of notebook 02 (SNV, SNV+CNV,
          SNV+CNV+SV), judged the same way as the main result in notebook 02b:
            (1) data with no real links  -- each set's real matrix shuffled within subtypes (every patient keeps
                its alteration count, every gene its frequency); every call is a false positive
            (2) planted links            -- the same shuffled data with links of known strength planted
                (script 22's scenario S2, its planting code); precision and recall at q < 0.05
            (3) TCGA                     -- each set's candidates re-tested in TCGA with the same test (script 14's
                rules and cached cBioPortal data); mutations only for the SNV set. TCGA has no fusion calls,
                so fusion pairs are not testable there
          Our test only (conditional, gene rates within subtypes); the method comparison stays in scripts 15 / 22.
Inputs  : data/processed/{cohort,events,testability}.parquet, comutation_pairs{_snv,,_snv_cnv_sv}.parquet,
          data/external/tcga/*.json (cached by script 14), scripts 14, 15 and 22 (shared code)
Outputs : data/processed/alteration_set_validation.parquet -- one row per pair and check (pair-level only,
          no patients); printed summary per set
Notes   : the 9 benchmark cancer types and patients on the panels covering 85% of each (as script 15), so every
          gene is tested in every patient; 6 replicates of the planted-link simulation. Seeded.
"""
import time

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

P = "data/processed/"
SETS = {"SNV": (["MUT"], "comutation_pairs_snv.parquet"),
        "SNV+CNV": (["MUT", "CNV"], "comutation_pairs.parquet"),
        "SNV+CNV+SV": (["MUT", "CNV", "SV"], "comutation_pairs_snv_cnv_sv.parquet")}
N_REP = 6

# Shared code: script 22's planting (which runs script 15's set-up: Phase 2 test functions, curveball, the 9 types)
src22 = open("scripts/analyses/22_planted_signal_simulation.py").read()
ns = {"__name__": "sim"}
exec(src22[:src22.index("rows, t0 = [], time.time()")], ns)
F, CTS, curveball_blocks = ns["F"], ns["CTS"], ns["curveball_blocks"]
ns["rng"] = np.random.default_rng(23)

co = pd.read_parquet(P + "cohort.parquet")
co = co[(co.TMB_GROUP == "non_hypermutated") & co.CNA_TESTED & co.MUT_TESTED]
ev = pd.read_parquet(P + "events.parquet", columns=["Sample_ID", "Hugo_Symbol", "Alteration_Type"])
tt = pd.read_parquet(P + "testability.parquet")
own = tt.MUT_TESTED | tt.SV_TESTED                     # as notebook 02, Step 1
cna_ok = tt.CNA_TESTED | ((tt.FEATURE_KIND != "gene") & own)


def matrix(set_name, ct):
    """Patients on the panels covering 85% of this cancer type, the set's genes tested on all of them."""
    types, f = SETS[set_name]
    cp = PAIRS[set_name]
    c = co[co.CANCER_TYPE == ct]
    covered = tt[own if set_name == "SNV" else own & cna_ok]
    per_panel = c.SEQ_ASSAY_ID.value_counts()
    big = per_panel.index[per_panel.cumsum() / per_panel.sum() <= 0.85].tolist() or per_panel.index[:1].tolist()
    common = set.intersection(*[set(covered.loc[covered.SEQ_ASSAY_ID == p, "Hugo_Symbol"]) for p in big])
    genes = sorted((set(cp.loc[cp.Cancer_Type == ct, "Gene_A"]) | set(cp.loc[cp.Cancer_Type == ct, "Gene_B"])) & common)
    c = c[c.SEQ_ASSAY_ID.isin(big)]
    e = ev[ev.Alteration_Type.isin(types) & ev.Sample_ID.isin(c.SAMPLE_ID) & ev.Hugo_Symbol.isin(genes)]
    X = (pd.crosstab(e.Sample_ID, e.Hugo_Symbol).reindex(index=c.SAMPLE_ID, columns=genes, fill_value=0) > 0).values
    keep = X.sum(0) >= 5; X = X[:, keep]; genes = [g for g, k in zip(genes, keep) if k]
    det = c.CANCER_TYPE_DETAILED.fillna("NA"); vc = det.value_counts()
    strata = det.where(det.map(vc) >= 50, "other").to_numpy()
    gi = {g: i for i, g in enumerate(genes)}
    pr = cp[(cp.Cancer_Type == ct) & cp.Gene_A.isin(gi) & cp.Gene_B.isin(gi)]
    return X, strata, genes, sorted((gi[a], gi[b]) for a, b in zip(pr.Gene_A, pr.Gene_B))


def ours(Y, strata, idx):
    """Our test (conditional, within subtypes) on matrix Y: p and direction per pair, BH q per cancer type."""
    m = F["ConditionalModel"](Y, np.ones_like(Y), strata)
    out = [F["pair_test"](m.joint(a, b), int((Y[:, a] & Y[:, b]).sum()))[:2] for a, b in idx]
    p = np.array([o[0] for o in out])
    return p, false_discovery_control(p.clip(1e-300, 1)), [o[1] for o in out]


PAIRS = {n: pd.read_parquet(P + f) for n, (_, f) in SETS.items()}
rows, t0 = [], time.time()
for set_name in SETS:
    for ct in CTS:
        X, strata, genes, idx = matrix(set_name, ct)
        if len(idx) < 2:
            continue
        name = lambda: dict(Gene_A=[genes[a] for a, _ in idx], Gene_B=[genes[b] for _, b in idx])
        bg = lambda: curveball_blocks(X, pd.factorize(strata)[0])
        p, q, d = ours(bg(), strata, idx)                                       # (1) no real links
        rows.append(pd.DataFrame(dict(check="no links", rep=0, p=p, q=q, direction=d, truth="none", fold=1.0, **name())))
        for rep in range(N_REP):                                               # (2) planted links
            Y, truth = ns["plant"](bg(), idx)
            p, q, d = ours(Y, strata, idx)
            rows.append(pd.DataFrame(dict(check="planted links", rep=rep, p=p, q=q, direction=d,
                                          truth=[truth.get(k, ("none", 1.0))[0] for k in range(len(idx))],
                                          fold=[truth.get(k, ("none", 1.0))[1] for k in range(len(idx))], **name())))
        for r in rows[-(N_REP + 1):]:
            r["alteration_set"], r["Cancer_Type"] = set_name, ct
        print(f"{set_name:<11} {ct:<28} patients {len(X):>6,} genes {len(genes):>3} pairs {len(idx):>4}  "
              f"{time.time() - t0:6.0f}s", flush=True)
sim = pd.concat(rows, ignore_index=True)

# (3) TCGA: script 14's set-up and cached data, re-tested per set
src14 = open("scripts/analyses/14_tcga_external_validation.py").read()
ns14 = {"__name__": "tcga"}
exec(src14[:src14.index("rows = []")], ns14)
cached, call, STUDIES = ns14["cached"], ns14["call"], ns14["STUDIES"]
TRUNC, INFRAME, kept_missense, role = (ns14[k] for k in ("TRUNC", "INFRAME", "kept_missense", "role"))


def missing(name):
    raise FileNotFoundError(f"TCGA cache {name} missing -- run script 14 first")


trows = []
for set_name, (types, _) in SETS.items():
    cp = PAIRS[set_name]
    for ct, studies in STUDIES.items():
        pairs = cp[cp.Cancer_Type == ct]
        if pairs.empty:
            continue
        genes = sorted(set(pairs.Gene_A) | set(pairs.Gene_B))
        alt_rows, samples, strata, tmb_all = [], [], [], {}
        for st in studies:
            sid = f"{st}_tcga_pan_can_atlas_2018"
            smp = cached(f"samples_{st}", lambda: missing(st))
            tmb_all.update({d["sampleId"]: float(d["value"]) for d in cached(f"tmb_{st}", lambda: missing(st))})
            for m in cached(f"mut_{st}_{ct}", lambda: missing(st)):
                g, mt, pc = m["gene"]["hugoGeneSymbol"], m.get("mutationType", ""), str(m.get("proteinChange", ""))
                if mt in TRUNC or mt in INFRAME or (mt == "Missense_Mutation" and (g, pc) in kept_missense):
                    alt_rows.append((m["sampleId"], g))
            if "CNV" in types:
                for c in cached(f"cna_{st}_{ct}", lambda: missing(st)):
                    g, v = c["gene"]["hugoGeneSymbol"], c["alteration"]
                    if {"TSG": v == -2, "ONC": v == 2}.get(role.get(g), True):
                        alt_rows.append((c["sampleId"], g))
            samples += smp; strata += [st] * len(smp)
        s = pd.Series(strata, index=samples)
        s = s[~s.index.duplicated()]
        s = s[[tmb_all.get(x, np.nan) < 10 for x in s.index]]                # non-hypermutated, as Phase 2
        alt = pd.DataFrame(alt_rows, columns=["s", "g"]).drop_duplicates()
        X = (pd.crosstab(alt.s, alt.g).reindex(index=s.index, columns=genes, fill_value=0) > 0).to_numpy()
        M = F["ConditionalModel"](X, np.ones_like(X, dtype=bool), s.to_numpy())
        gi = {g: i for i, g in enumerate(genes)}
        for a, b in sorted(zip(pairs.Gene_A.map(gi), pairs.Gene_B.map(gi))):
            jq = M.joint(a, b)
            if jq.sum() < 5:                                                    # same power gate
                continue
            p, d, e = F["pair_test"](jq, int((X[:, a] & X[:, b]).sum()))
            trows.append(dict(alteration_set=set_name, Cancer_Type=ct, Gene_A=genes[a], Gene_B=genes[b],
                              tcga_p=p, tcga_direction=d))
tcga = pd.DataFrame(trows)
tcga = tcga.merge(pd.concat([PAIRS[n].assign(alteration_set=n)[["alteration_set", "Cancer_Type", "Gene_A", "Gene_B",
                                                                 "direction", "q_value", "interaction_candidate"]]
                             for n in SETS]), on=["alteration_set", "Cancer_Type", "Gene_A", "Gene_B"], how="left")
tcga["check"] = "TCGA"

out = pd.concat([sim, tcga], ignore_index=True)
out.to_parquet(P + "alteration_set_validation.parquet", index=False)

print("\nPer alteration set (our test, q < 0.05)")
for set_name in SETS:
    n = sim[(sim.alteration_set == set_name) & (sim.check == "no links")]
    s = sim[(sim.alteration_set == set_name) & (sim.check == "planted links")]
    call = s.q < .05
    tp = call & (s.truth != "none") & (s.direction == s.truth)
    t = tcga[(tcga.alteration_set == set_name) & tcga.interaction_candidate]
    same = t.tcga_direction == t.direction
    print(f"  {set_name:<11} no links: p<0.05 {(n.p < .05).mean():5.1%}, q<0.05 {(n.q < .05).mean():5.1%} | planted: "
          f"recall {tp.sum() / max((s.truth != 'none').sum(), 1):5.1%}, precision {tp.sum() / max(call.sum(), 1):5.1%} | "
          f"TCGA: candidates testable {len(t):>4}, same direction {same.mean():5.1%}, "
          f"p<0.05 same direction {((t.tcga_p < .05) & same).mean():5.1%}, significant opposite {((t.tcga_p < .05) & ~same).sum()}")
