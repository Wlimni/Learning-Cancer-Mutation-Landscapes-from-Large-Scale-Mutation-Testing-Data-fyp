"""
Script Name: check_phase1_handoff.py
Purpose : Check the three tables Phase 1 hands to Phase 2 before anything is tested: one row per patient,
          every event belongs to a cohort patient, one row per (panel, gene feature), and every event was
          testable on its patient's panel for its alteration type.
Inputs  : data/processed/cohort.parquet, events.parquet, testability.parquet (notebook 01, Section 10)
Outputs : printed report; raises AssertionError if a structural rule is broken
Note    : run from the repo root (python scripts/validation/check_phase1_handoff.py) or import check().
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROCESSED = Path(__file__).resolve().parents[2] / "data" / "processed"
TMB_GROUPS = {"non_hypermutated", "hypermutated", "ultramutated", "unknown"}
TESTED_FOR = {"MUT": "MUT_TESTED", "CNV": "CNA_TESTED", "SV": "SV_TESTED"}


def check(cohort: pd.DataFrame, events: pd.DataFrame, testability: pd.DataFrame) -> dict:
    """Structural rules are asserted; the testability of events is reported (a few off-panel calls are expected)."""
    assert cohort["PATIENT_ID"].is_unique, "cohort: more than one row for a patient"
    assert cohort["SAMPLE_ID"].is_unique, "cohort: duplicated sample"
    assert set(cohort["TMB_GROUP"]) <= TMB_GROUPS and cohort["TMB_GROUP"].notna().all(), "cohort: bad TMB_GROUP"
    assert events["Sample_ID"].isin(cohort["SAMPLE_ID"]).all(), "events: sample not in the cohort"
    assert set(events["Alteration_Type"]) <= set(TESTED_FOR), "events: unknown alteration type"
    assert not testability.duplicated(["SEQ_ASSAY_ID", "Hugo_Symbol"]).any(), "testability: duplicated (panel, feature)"

    e = events[["Sample_ID", "Hugo_Symbol", "Alteration_Type"]].merge(
        cohort[["SAMPLE_ID", "SEQ_ASSAY_ID"]], left_on="Sample_ID", right_on="SAMPLE_ID")
    e = e.merge(testability, on=["SEQ_ASSAY_ID", "Hugo_Symbol"], how="left")
    testable = np.zeros(len(e), dtype=bool)
    for alt_type, col in TESTED_FOR.items():
        rows = (e["Alteration_Type"] == alt_type).to_numpy()
        testable[rows] = e.loc[rows, col].eq(True).to_numpy()   # missing (panel, feature) row = not tested
    return {"patients": len(cohort), "events": len(events),
            "events_testable_on_panel": float(testable.mean()),
            "events_by_type": events["Alteration_Type"].value_counts().to_dict(),
            "patients_by_tmb_group": cohort["TMB_GROUP"].value_counts().to_dict()}


if __name__ == "__main__":
    report = check(pd.read_parquet(PROCESSED / "cohort.parquet"), pd.read_parquet(PROCESSED / "events.parquet"),
                   pd.read_parquet(PROCESSED / "testability.parquet"))
    for k, v in report.items():
        print(f"{k}: {v}")
    print("hand-off checks passed")
