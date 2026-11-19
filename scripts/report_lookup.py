"""Script Name: report_lookup.py.
Purpose: search local GENIE data for observed subsets of report findings.
Inputs: cancer type, reviewed findings, local GENIE files. Outputs: aggregate tables.
Notes: Phase 2 evidence is reported separately at gene level.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = ROOT / "data" / "processed"


def as_feature(f) -> tuple[str, str | None]:
    return (f, None) if isinstance(f, str) else (f[0], f[1])


def label(f) -> str:
    g, d = as_feature(f)
    return g if d is None else f"{g} {d.replace('p.', '')}"


def carriers_of(sub_alt: pd.DataFrame, f) -> set:
    gene, detail = as_feature(f)
    hit = sub_alt["Hugo_Symbol"] == gene
    if detail in ("AMP", "DEL"):
        hit &= sub_alt["CNA_Direction"] == {"AMP": "Amplification", "DEL": "Deep_Deletion"}[detail]
    elif detail is not None:
        hit &= sub_alt["HGVSp_Short"] == detail
    return set(sub_alt.loc[hit, "Sample_ID"])


def scan_report(cancer_type, driver_findings, uncertain_findings, *, subtype=None, min_observed=1, max_combinations=None):
    """Return finding coverage, observed subsets of size >= 2, and Phase 2 gene-pair evidence."""

    if not isinstance(min_observed, int) or min_observed < 1:
        raise ValueError("min_observed must be a positive integer")

    cn_states = {"AMP": {2}, "DEL": {-2}, "GAIN": {1, 2}, "LOSS": {-1, -2}}
    groups = {}
    for group, findings in [("driver", driver_findings), ("uncertain", uncertain_findings)]:
        for feature in findings:
            gene, detail = as_feature(feature)
            gene = gene.strip()
            if detail is not None:
                detail = detail.strip().removeprefix("p.").replace("Ter", "*")
            groups.setdefault((gene, detail), set()).add(group)
    if not groups or any(not gene or detail == "" for gene, detail in groups):
        raise ValueError("Provide at least one valid gene or (gene, detail) finding")

    patients = pd.read_parquet(PROCESSED_DIR / "clinical_tidy.parquet")
    selected = patients["IS_REPRESENTATIVE"] & patients["CANCER_TYPE"].eq(cancer_type)
    if subtype is not None:
        selected &= patients["CANCER_TYPE_DETAILED"].eq(subtype)
    patients = patients.loc[selected].set_index("SAMPLE_ID")
    if patients.empty:
        raise ValueError(f"No representative samples for {cancer_type!r}, subtype {subtype!r}")

    raw = ROOT / "data" / "raw"
    mutation_targets = {f for f in groups if f[1] is not None and f[1] not in cn_states
                        and f[0] != "TERT_promoter"}
    cnv_genes = {gene for gene, detail in groups if detail in cn_states}
    events_needed = [f for f in groups if f not in mutation_targets and f[1] not in cn_states]
    carriers = {f: set() for f in groups}
    tested = {}
    if mutation_targets:
        genes = {gene for gene, _ in mutation_targets}
        cols = ["Hugo_Symbol", "HGVSp_Short", "Tumor_Sample_Barcode", "Mutation_Status"]
        # The raw somatic calls include protein changes excluded by the Phase 1 driver filter.
        for chunk in pd.read_csv(raw / "data_mutations_extended.txt", sep="\t", usecols=cols,
                                 dtype=str, chunksize=250_000):
            rows = chunk[chunk["Hugo_Symbol"].isin(genes)
                         & chunk["Tumor_Sample_Barcode"].isin(patients.index)
                         & chunk["Mutation_Status"].isin(["Somatic", "SOMATIC"])]
            variants = rows["HGVSp_Short"].str.removeprefix("p.").str.replace("Ter", "*", regex=False)
            for (gene, variant), hits in rows.groupby(["Hugo_Symbol", variants]):
                if (gene, variant) in mutation_targets:
                    carriers[gene, variant].update(hits["Tumor_Sample_Barcode"])

    if cnv_genes:
        with (raw / "data_CNA.txt").open() as handle:
            samples = handle.readline().rstrip("\r\n").split("\t")[1:]
            columns = [(i, sample) for i, sample in enumerate(samples, 1) if sample in patients.index]
            calls_by_gene = {}
            for row in handle:
                gene = row.split("\t", 1)[0]
                if gene in cnv_genes:
                    cells = row.rstrip("\r\n").split("\t")
                    calls_by_gene[gene] = pd.to_numeric(
                        pd.Series({s: cells[i] for i, s in columns}, dtype="object"), errors="coerce").dropna()
        for f in groups:
            gene, detail = f
            if detail in cn_states:
                calls = calls_by_gene.get(gene, pd.Series(dtype=float))
                tested[f] = set(calls.index)
                carriers[f] = set(calls.index[calls.isin(cn_states[detail])])

    def coverage_for(name):
        return pd.read_parquet(PROCESSED_DIR / f"{name}.parquet")

    mutation_cov = coverage_for("panel_gene_coverage")
    if events_needed:
        cols = ["Sample_ID", "Hugo_Symbol", "Alteration_Type", "HGVSp_Short", "CNA_Direction"]
        events = pd.read_parquet(PROCESSED_DIR / "alterations_long.parquet", columns=cols)
        events = events[events["Sample_ID"].isin(patients.index)]
        for f in events_needed:
            carriers[f] = carriers_of(events, f)

    fusion_cov = coverage_for("sv_gene_panel_coverage") if any(g.endswith("_fusion") for g, _ in groups) else None
    cna_cov = coverage_for("cna_gene_panel_coverage") if any(d is None and not g.endswith("_fusion") for g, d in groups) else None
    for f in groups:
        gene, detail = f
        if f not in tested:
            cov, flag = (fusion_cov, "SV_TESTED") if gene.endswith("_fusion") else (mutation_cov, "MUT_TESTED")
            panels = cov.loc[cov["Hugo_Symbol"].eq(gene), "SEQ_ASSAY_ID"]
            mask = patients[flag] & patients["SEQ_ASSAY_ID"].isin(panels)
            if detail is None and cna_cov is not None and not gene.endswith("_fusion"):
                cna_panels = cna_cov.loc[cna_cov["Hugo_Symbol"].eq(gene), "SEQ_ASSAY_ID"]
                mask &= patients["CNA_TESTED"] & patients["SEQ_ASSAY_ID"].isin(cna_panels)
            tested[f] = set(patients.index[mask])
        carriers[f] &= tested[f]

    names = {f: label(f) for f in groups}
    findings = pd.DataFrame([
        {"finding": names[f], "report_group": " + ".join(sorted(groups[f])),
         "n_observed": len(carriers[f]), "n_tested": len(tested[f])}
        for f in groups
    ])

    # Once an intersection falls below the count floor, none of its supersets can meet it.
    active = [f for f in groups if len(carriers[f]) >= min_observed]
    rows = []

    def visit(prefix, jointly_tested, observed, start):
        for i in range(start, len(active)):
            f = active[i]
            hits = observed & carriers[f]
            if len(hits) < min_observed:
                continue
            combo = (*prefix, f)
            eligible = jointly_tested & tested[f]
            if len(combo) >= 2:
                if max_combinations is not None and len(rows) >= max_combinations:
                    raise ValueError(
                        f"More than {max_combinations:,} combinations meet this count. "
                        "Increase the minimum patients observed and try again."
                    )
                rows.append((len(combo), tuple(names[x] for x in combo), len(hits), len(eligible),
                             100 * len(hits) / len(eligible)))
            visit(combo, eligible, hits, i + 1)

    all_patients = set(patients.index)
    visit((), all_patients, all_patients, 0)
    observed_combos = pd.DataFrame(rows, columns=["size", "combination", "n_observed", "n_tested", "pct"])
    observed_combos = observed_combos.sort_values(["size", "n_observed", "pct"], ascending=[True, False, False])

    genes = {gene for gene, _ in groups}
    pairs = pd.read_parquet(PROCESSED_DIR / "comutation_pairs.parquet")
    matched = (pairs["Cancer_Type"].eq(cancer_type) & pairs["interaction_candidate"]
               & pairs["Gene_A"].isin(genes) & pairs["Gene_B"].isin(genes))
    evidence = pairs.loc[matched, ["Gene_A", "Gene_B", "direction", "q_value", "replicated"]]
    return {"findings": findings, "combinations": observed_combos,
            "phase2_pairs": evidence.sort_values("q_value").reset_index(drop=True)}
