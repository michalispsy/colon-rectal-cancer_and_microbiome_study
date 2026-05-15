#!/usr/bin/env python3
"""
Step 5.2 — Cross-validation of SHAP features with Phase 3 differential abundance results.

Purpose
-------
Compare model-derived important species from Phase 5.1 SHAP analysis against
Phase 3 statistically supported differentially abundant species.

Implemented models only:
  - All        = previous Model A, universal CRC vs Control
  - OnlyMale   = previous Model B, male-only CRC vs Control
  - OnlyFemale = previous Model C, female-only CRC vs Control

This script intentionally ignores Adenoma / Models D,E,F.

Expected inputs
---------------
SHAP directory from Step 5.1, e.g.:
  OUTPUTS/phase_5/5.1_shap

Phase 3 directory, e.g.:
  OUTPUTS/phase_3
with files like:
  3.1/step3_1_blocked_wilcoxon_significant_results.csv
  3.2/step3_2_cmh_significant_results.csv
  3.3/step3_3_meta_analysis_significant_results.csv
  3.5/step3_5_biomarker_categories.csv      [optional]

Outputs
-------
  OUTPUTS/phase_5/5.2_phase3_validation/
    step5_2_summary_report.md
    step5_2_overlap_summary.csv
    <Model>_top<TOPK>_SHAP_vs_Phase3.csv
    <Model>_all_SHAP_species_vs_Phase3.csv
    phase3_reference_sets.csv
    figures/*.png
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODEL_ALIASES = {
    "All": ["All", "Model_A"],
    "OnlyMale": ["OnlyMale", "Model_B"],
    "OnlyFemale": ["OnlyFemale", "Model_C"],
}

# Per-model Phase 3 comparison filters.
PHASE3_FILTERS = {
    "All": {
        "3.1": {"comparison": {"CRC_vs_Control_all"}},
        "3.2": {"comparison": {"CRC_vs_Control__all"}},
        "3.3": {"comparison": {"CRC_vs_Control"}, "stratum": {"All"}},
    },
    "OnlyMale": {
        "3.1": {"comparison": {"CRC_vs_Control_Male"}},
        "3.2": {"comparison": {"CRC_vs_Control__Male"}},
        "3.3": {"comparison": {"CRC_vs_Control"}, "stratum": {"Male"}},
    },
    "OnlyFemale": {
        "3.1": {"comparison": {"CRC_vs_Control_Female"}},
        "3.2": {"comparison": {"CRC_vs_Control__Female"}},
        "3.3": {"comparison": {"CRC_vs_Control"}, "stratum": {"Female"}},
    },
}


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def clean_species_name(x: object) -> str:
    """Normalize feature/species names for robust set overlap."""
    if pd.isna(x):
        return ""
    s = str(x).strip()
    # Remove common feature prefixes used across scripts.
    s = re.sub(r"^(PA|rCLR|CLR|clr|pa)__", "", s)
    s = re.sub(r"^(PA|rCLR|CLR|clr|pa)_", "", s)
    s = s.replace("s__", "")
    # Normalize separators/whitespace.
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s


def canonical_key(x: object) -> str:
    """Case-insensitive key for set operations."""
    return clean_species_name(x).lower()


def find_existing_file(candidates: Iterable[Path], required: bool = True) -> Optional[Path]:
    for p in candidates:
        if p.exists():
            return p
    if required:
        msg = "Could not find any of:\n" + "\n".join(str(p) for p in candidates)
        raise FileNotFoundError(msg)
    return None


def find_model_species_consensus(shap_dir: Path, display_model: str) -> Path:
    """Find species-level SHAP consensus file for All / OnlyMale / OnlyFemale."""
    candidates: List[Path] = []
    for alias in MODEL_ALIASES[display_model]:
        model_dir = shap_dir / alias
        if model_dir.exists():
            candidates.extend(sorted(model_dir.glob("*species_level_consensus.csv")))
        candidates.extend(sorted(shap_dir.glob(f"**/{alias}*species_level_consensus.csv")))
    # Avoid duplicates while preserving order.
    seen = set()
    unique = []
    for p in candidates:
        if p not in seen:
            unique.append(p)
            seen.add(p)
    return find_existing_file(unique, required=True)  # type: ignore[arg-type]


def load_shap_species(shap_dir: Path, display_model: str) -> pd.DataFrame:
    path = find_model_species_consensus(shap_dir, display_model)
    df = pd.read_csv(path)
    if "Species" not in df.columns:
        raise ValueError(f"{path} does not contain a Species column")
    df = df.copy()
    df["DisplayModel"] = display_model
    df["SourceFile"] = str(path)
    df["Species_clean"] = df["Species"].map(clean_species_name)
    df["Species_key"] = df["Species"].map(canonical_key)

    # Make rank and importance columns predictable even if names change.
    if "SpeciesRank" not in df.columns:
        if "ConsensusRank" in df.columns:
            df["SpeciesRank"] = df["ConsensusRank"]
        else:
            df["SpeciesRank"] = np.arange(1, len(df) + 1)
    if "Species_MeanAbsSHAP" not in df.columns:
        for col in ["MeanAbsSHAP_Mean", "MeanAbsSHAP", "ConsensusScore"]:
            if col in df.columns:
                df["Species_MeanAbsSHAP"] = df[col]
                break
        else:
            df["Species_MeanAbsSHAP"] = np.nan
    return df.sort_values("SpeciesRank").reset_index(drop=True)


def phase3_file_candidates(phase3_dir: Path, subdir: str, filename: str) -> List[Path]:
    return [
        phase3_dir / subdir / filename,
        phase3_dir / filename,
        phase3_dir / f"phase_{subdir}" / filename,
        phase3_dir / f"step{subdir.replace('.', '_')}" / filename,
    ] + sorted(phase3_dir.glob(f"**/{filename}"))


def load_phase3_tables(phase3_dir: Path) -> Dict[str, pd.DataFrame]:
    files = {
        "3.1": ("3.1", "step3_1_blocked_wilcoxon_significant_results.csv"),
        "3.2": ("3.2", "step3_2_cmh_significant_results.csv"),
        "3.3": ("3.3", "step3_3_meta_analysis_significant_results.csv"),
        "3.5": ("3.5", "step3_5_biomarker_categories.csv"),
    }
    out: Dict[str, pd.DataFrame] = {}
    for key, (subdir, filename) in files.items():
        required = key in {"3.1", "3.2", "3.3"}
        path = find_existing_file(phase3_file_candidates(phase3_dir, subdir, filename), required=required)
        if path is None:
            continue
        df = pd.read_csv(path)
        if "species" not in df.columns and "Species" in df.columns:
            df = df.rename(columns={"Species": "species"})
        if "species" not in df.columns:
            raise ValueError(f"{path} does not contain a species/Species column")
        df = df.copy()
        df["Phase3Step"] = key
        df["Phase3SourceFile"] = str(path)
        df["Species_clean"] = df["species"].map(clean_species_name)
        df["Species_key"] = df["species"].map(canonical_key)
        out[key] = df
    return out


def apply_filters(df: pd.DataFrame, filters: Dict[str, Set[str]]) -> pd.DataFrame:
    out = df.copy()
    for col, allowed in filters.items():
        if col not in out.columns:
            # Missing stratum/comparison means no supported rows for this model-step.
            return out.iloc[0:0].copy()
        out = out[out[col].astype(str).isin(allowed)]
    return out.copy()


def build_reference_sets(phase3: Dict[str, pd.DataFrame]) -> Tuple[Dict[str, Dict[str, Set[str]]], pd.DataFrame]:
    """Create model-specific Phase 3 significant species sets."""
    sets: Dict[str, Dict[str, Set[str]]] = {}
    rows = []

    for model in ["All", "OnlyMale", "OnlyFemale"]:
        sets[model] = {}
        for step in ["3.1", "3.2", "3.3"]:
            df = phase3.get(step, pd.DataFrame())
            filtered = apply_filters(df, PHASE3_FILTERS[model][step]) if not df.empty else df
            species_set = set(filtered["Species_key"].dropna()) if not filtered.empty else set()
            sets[model][step] = species_set
            for _, r in filtered.iterrows():
                rows.append({
                    "DisplayModel": model,
                    "Phase3Step": step,
                    "Species": r.get("Species_clean", r.get("species", "")),
                    "Species_key": r.get("Species_key", ""),
                    "comparison": r.get("comparison", ""),
                    "stratum": r.get("stratum", ""),
                    "q_value": r.get("q_value", r.get("q_value_bh", r.get("p_adj_BH", np.nan))),
                    "direction": r.get("direction", ""),
                    "feature": r.get("feature", ""),
                })

        formal_union = sets[model]["3.1"] | sets[model]["3.2"] | sets[model]["3.3"]
        sets[model]["formal_union"] = formal_union
        sets[model]["strong_core_3_3"] = sets[model]["3.3"]
        sets[model]["supported_by_2plus"] = {
            s for s in formal_union
            if sum(s in sets[model][step] for step in ["3.1", "3.2", "3.3"]) >= 2
        }

    ref_df = pd.DataFrame(rows)
    return sets, ref_df


def load_sequence_categories(phase3: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    df = phase3.get("3.5", pd.DataFrame()).copy()
    if df.empty:
        return df
    keep = [
        "Species_key", "Species_clean", "categories", "is_early", "is_late", "is_progression",
        "is_ordered_increasing", "is_ordered_decreasing", "best_trend_direction", "best_trend_q",
        "q_crc_vs_control", "q_adenoma_vs_control", "q_crc_vs_adenoma",
    ]
    keep = [c for c in keep if c in df.columns]
    return df[keep].drop_duplicates("Species_key")


def validate_model(
    shap_df: pd.DataFrame,
    model: str,
    ref_sets: Dict[str, Dict[str, Set[str]]],
    sequence_df: pd.DataFrame,
    top_k: int,
    output_dir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = shap_df.copy()

    for step, label in [
        ("3.1", "In_Phase3_1_BlockedWilcoxon_rCLR"),
        ("3.2", "In_Phase3_2_CMH_PA"),
        ("3.3", "In_Phase3_3_MetaAnalysis_rCLR"),
    ]:
        df[label] = df["Species_key"].isin(ref_sets[model][step])

    df["In_Any_Formal_Phase3"] = df["Species_key"].isin(ref_sets[model]["formal_union"])
    df["In_2plus_Formal_Phase3"] = df["Species_key"].isin(ref_sets[model]["supported_by_2plus"])
    df["In_MetaAnalysis_Core"] = df["Species_key"].isin(ref_sets[model]["strong_core_3_3"])
    df["N_Formal_Phase3_Supports"] = (
        df["In_Phase3_1_BlockedWilcoxon_rCLR"].astype(int)
        + df["In_Phase3_2_CMH_PA"].astype(int)
        + df["In_Phase3_3_MetaAnalysis_rCLR"].astype(int)
    )

    def support_level(row: pd.Series) -> str:
        if bool(row["In_MetaAnalysis_Core"]):
            return "Strong: reproduced by random-effects meta-analysis"
        if int(row["N_Formal_Phase3_Supports"]) >= 2:
            return "Strong: supported by at least two formal Phase 3 tests"
        if bool(row["In_Any_Formal_Phase3"]):
            return "Moderate: supported by one formal Phase 3 test"
        return "ML-only: no formal Phase 3 overlap"

    df["Phase3_Support_Level"] = df.apply(support_level, axis=1)

    if not sequence_df.empty:
        df = df.merge(sequence_df, on="Species_key", how="left", suffixes=("", "_phase35"))

    df.to_csv(output_dir / f"{model}_all_SHAP_species_vs_Phase3.csv", index=False)
    df.head(top_k).to_csv(output_dir / f"{model}_top{top_k}_SHAP_vs_Phase3.csv", index=False)

    summary_rows = []
    for n in [10, 20, 50, top_k]:
        n = min(n, len(df))
        if n <= 0:
            continue
        top = df.head(n)
        summary_rows.append({
            "DisplayModel": model,
            "Top_N": n,
            "N_SHAP_species": len(top),
            "N_overlap_any_formal_Phase3": int(top["In_Any_Formal_Phase3"].sum()),
            "Percent_overlap_any_formal_Phase3": 100 * float(top["In_Any_Formal_Phase3"].mean()),
            "N_overlap_Phase3_1_BlockedWilcoxon": int(top["In_Phase3_1_BlockedWilcoxon_rCLR"].sum()),
            "N_overlap_Phase3_2_CMH": int(top["In_Phase3_2_CMH_PA"].sum()),
            "N_overlap_Phase3_3_MetaAnalysis": int(top["In_Phase3_3_MetaAnalysis_rCLR"].sum()),
            "N_supported_by_2plus_formal_tests": int(top["In_2plus_Formal_Phase3"].sum()),
            "RedFlag_No_Overlap": bool(top["In_Any_Formal_Phase3"].sum() == 0),
        })
    summary = pd.DataFrame(summary_rows).drop_duplicates(["DisplayModel", "Top_N"])
    return df, summary


def make_overlap_plot(summary: pd.DataFrame, output_path: Path, top_n: int) -> None:
    df = summary[summary["Top_N"] == top_n].copy()
    if df.empty:
        return
    x = np.arange(len(df))
    width = 0.22
    plt.figure(figsize=(9, 5))
    plt.bar(x - width, df["N_overlap_Phase3_1_BlockedWilcoxon"], width, label="3.1 Blocked Wilcoxon")
    plt.bar(x, df["N_overlap_Phase3_2_CMH"], width, label="3.2 CMH")
    plt.bar(x + width, df["N_overlap_Phase3_3_MetaAnalysis"], width, label="3.3 Meta-analysis")
    plt.xticks(x, df["DisplayModel"])
    plt.ylabel(f"Overlap count among top {top_n} SHAP species")
    plt.title(f"SHAP vs Phase 3 overlap, top {top_n}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def df_to_md(df: pd.DataFrame, max_rows: Optional[int] = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    try:
        return df.to_markdown(index=False)
    except Exception:
        # Fallback if tabulate is unavailable.
        return df.to_string(index=False)


def write_report(
    output_dir: Path,
    overlap_summary: pd.DataFrame,
    model_tables: Dict[str, pd.DataFrame],
    ref_df: pd.DataFrame,
    top_k: int,
    phase3_loaded: Dict[str, pd.DataFrame],
) -> None:
    lines: List[str] = []
    lines.append("# Step 5.2 — Cross-validation of SHAP Features with Phase 3\n")
    lines.append("Models included: **All**, **OnlyMale**, **OnlyFemale**. Adenoma / Models D,E,F were intentionally excluded.\n")
    lines.append("## Inputs used\n")
    lines.append("Formal Phase 3 validation sources:\n")
    lines.append("- 3.1 blocked Wilcoxon / Van Elteren significant rCLR species\n")
    lines.append("- 3.2 Cochran-Mantel-Haenszel significant PA species\n")
    lines.append("- 3.3 random-effects meta-analysis significant rCLR species\n")
    if "3.5" in phase3_loaded:
        lines.append("- 3.5 adenoma-carcinoma biomarker categories were added as optional biological annotation\n")
    lines.append("\n")

    lines.append("## Overlap summary\n")
    lines.append(df_to_md(overlap_summary))
    lines.append("\n")

    for model, df in model_tables.items():
        lines.append(f"## {model}\n")
        top = df.head(top_k).copy()
        cols = [
            "Species_clean", "SpeciesRank", "Species_MeanAbsSHAP", "Dominant_Signal",
            "In_Phase3_1_BlockedWilcoxon_rCLR", "In_Phase3_2_CMH_PA",
            "In_Phase3_3_MetaAnalysis_rCLR", "N_Formal_Phase3_Supports",
            "Phase3_Support_Level",
        ]
        extra = ["categories", "is_early", "is_late", "is_progression"]
        cols = [c for c in cols + extra if c in top.columns]
        lines.append(f"Top {top_k} SHAP species with Phase 3 support:\n")
        lines.append(df_to_md(top[cols]))
        lines.append("\n")

        n_any = int(top["In_Any_Formal_Phase3"].sum())
        n_meta = int(top["In_Phase3_3_MetaAnalysis_rCLR"].sum())
        if n_any == 0:
            lines.append("> Red flag: no top SHAP species overlapped with the formal Phase 3 species lists.\n\n")
        elif n_meta > 0:
            lines.append(f"> Interpretation: {n_any}/{len(top)} top SHAP species overlap with formal Phase 3 results, including {n_meta} supported by meta-analysis.\n\n")
        else:
            lines.append(f"> Interpretation: {n_any}/{len(top)} top SHAP species overlap with formal Phase 3 results, but none with the meta-analysis core list.\n\n")

    if not ref_df.empty:
        lines.append("## Phase 3 reference set sizes\n")
        ref_counts = (
            ref_df.groupby(["DisplayModel", "Phase3Step"], as_index=False)
            .agg(N_species=("Species_key", "nunique"))
        )
        lines.append(df_to_md(ref_counts))
        lines.append("\n")

    (output_dir / "step5_2_summary_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 5.2 SHAP vs Phase 3 cross-validation")
    p.add_argument("--shap-dir", required=True, help="Path to OUTPUTS/phase_5/5.1_shap")
    p.add_argument("--phase3-dir", required=True, help="Path to OUTPUTS/phase_3")
    p.add_argument("--output-dir", default="OUTPUTS/phase_5/5.2_phase3_validation")
    p.add_argument("--top-k", type=int, default=20, help="Top SHAP species to emphasize in per-model reports")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    shap_dir = Path(args.shap_dir)
    phase3_dir = Path(args.phase3_dir)
    output_dir = Path(args.output_dir)
    figures_dir = output_dir / "figures"
    safe_mkdir(output_dir)
    safe_mkdir(figures_dir)

    phase3 = load_phase3_tables(phase3_dir)
    ref_sets, ref_df = build_reference_sets(phase3)
    seq_df = load_sequence_categories(phase3)
    if not ref_df.empty:
        ref_df.to_csv(output_dir / "phase3_reference_sets.csv", index=False)

    model_tables: Dict[str, pd.DataFrame] = {}
    summaries: List[pd.DataFrame] = []

    for model in ["All", "OnlyMale", "OnlyFemale"]:
        shap_df = load_shap_species(shap_dir, model)
        validated, summary = validate_model(
            shap_df=shap_df,
            model=model,
            ref_sets=ref_sets,
            sequence_df=seq_df,
            top_k=args.top_k,
            output_dir=output_dir,
        )
        model_tables[model] = validated
        summaries.append(summary)

    overlap_summary = pd.concat(summaries, ignore_index=True)
    overlap_summary.to_csv(output_dir / "step5_2_overlap_summary.csv", index=False)
    make_overlap_plot(overlap_summary, figures_dir / f"step5_2_overlap_top{args.top_k}.png", top_n=args.top_k)

    write_report(
        output_dir=output_dir,
        overlap_summary=overlap_summary,
        model_tables=model_tables,
        ref_df=ref_df,
        top_k=args.top_k,
        phase3_loaded=phase3,
    )

    print("\nStep 5.2 complete.")
    print(f"Outputs written to: {output_dir}")
    print("Main report:", output_dir / "step5_2_summary_report.md")
    print("Overlap summary:", output_dir / "step5_2_overlap_summary.csv")


if __name__ == "__main__":
    main()
