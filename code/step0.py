"""
STEP 0 — DATA DISCOVERY & STRUCTURE IDENTIFICATION

This script performs an initial exploration of the project directory to identify
usable datasets and determine their roles in the analysis pipeline.

Main functionalities:
- Scans multiple directories (curatedMetagenomicData, harmonized, processed_data, etc.)
- Automatically detects tabular files (CSV/TSV)
- Infers file types (metadata vs abundance matrix) using heuristics
- Identifies key columns:
  - sample IDs (join keys)
  - condition / disease labels (e.g., CRC, Control)
  - sex / gender
  - study / cohort (for LOSO later)
- Prints summaries of each dataset (shape, columns, data types)
- Suggests the most likely candidates for:
  - feature matrix (microbiome abundances)
  - metadata table (labels + covariates)
- Saves a full dataset inventory to JSON for reproducibility

Blueprint steps covered:
- Data understanding and dataset identification
- Selection of feature matrix and metadata table
- Preparation for merging and downstream analysis
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import pandas as pd


# ---------------------------
# Config
# ---------------------------

PROJECT_ROOT = Path(".").resolve()
SEARCH_DIRS = [
    "curatedMetagenomicData",
    "harmonized",
    "processed_data",
    "species_level",
    "supplementary_tables",
]

TEXT_EXTENSIONS = {".csv", ".tsv", ".txt"}
MAX_ROWS_TO_READ = 200
MAX_COLS_TO_SHOW = 25
OUTPUT_JSON = "dataset_inventory.json"


# ---------------------------
# Helpers
# ---------------------------

@dataclass
class FileSummary:
    path: str
    delimiter: str
    n_rows: Optional[int]
    n_cols: Optional[int]
    columns_sample: list[str]
    index_name: Optional[str]
    first_col: Optional[str]
    likely_type: str
    sample_id_like_cols: list[str]
    sex_like_cols: list[str]
    target_like_cols: list[str]
    study_like_cols: list[str]
    notes: list[str]


def sniff_delimiter(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            head = "".join([f.readline() for _ in range(5)])
    except Exception:
        return ","

    if head.count("\t") > head.count(","):
        return "\t"
    return ","


def safe_read_table(path: Path, nrows: int = MAX_ROWS_TO_READ) -> Optional[pd.DataFrame]:
    delimiter = sniff_delimiter(path)
    try:
        df = pd.read_csv(path, sep=delimiter, nrows=nrows, low_memory=False)
        return df
    except Exception:
        return None


def normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def find_matching_columns(columns: list[str], patterns: list[str]) -> list[str]:
    matches = []
    normalized = {col: normalize_name(col) for col in columns}
    for col, norm in normalized.items():
        if any(p in norm for p in patterns):
            matches.append(col)
    return matches


def infer_file_type(df: pd.DataFrame, path: Path) -> tuple[str, list[str]]:
    notes = []
    cols = [str(c) for c in df.columns]

    sample_id_like = find_matching_columns(
        cols,
        ["sample", "sample_id", "run", "accession", "id"]
    )
    sex_like = find_matching_columns(
        cols,
        ["sex", "gender", "male", "female"]
    )
    target_like = find_matching_columns(
        cols,
        ["condition", "disease", "status", "class", "label", "crc", "adenoma", "control"]
    )
    study_like = find_matching_columns(
        cols,
        ["study", "cohort", "dataset", "project", "batch", "source"]
    )

    n_rows, n_cols = df.shape
    numeric_fraction = (df.dtypes.apply(lambda x: pd.api.types.is_numeric_dtype(x)).sum() / max(1, n_cols))

    likely_type = "unknown"

    # Heuristic: metadata table
    if (sex_like or target_like or study_like) and numeric_fraction < 0.7:
        likely_type = "metadata"
        notes.append("Contains phenotype/study-like columns.")
    # Heuristic: abundance table
    elif numeric_fraction > 0.7 and n_cols > 20:
        likely_type = "abundance_matrix"
        notes.append("Mostly numeric with many columns; looks like a feature matrix.")
    # Heuristic: supplementary
    elif "supp" in path.name.lower() or "table" in path.name.lower():
        likely_type = "supplementary"
        notes.append("Filename suggests supplementary material.")

    # More hints
    if any(k in path.name.lower() for k in ["metadata", "meta"]):
        likely_type = "metadata"
        notes.append("Filename suggests metadata.")
    if any(k in path.name.lower() for k in ["abundance", "feature", "species", "genus", "relative"]):
        if likely_type == "unknown":
            likely_type = "abundance_matrix"
        notes.append("Filename suggests taxa abundance/features.")

    return likely_type, notes


def summarize_file(path: Path) -> Optional[FileSummary]:
    if path.suffix.lower() not in TEXT_EXTENSIONS:
        return None

    df = safe_read_table(path)
    if df is None:
        return None

    cols = [str(c) for c in df.columns]
    likely_type, notes = infer_file_type(df, path)

    return FileSummary(
        path=str(path.relative_to(PROJECT_ROOT)),
        delimiter=sniff_delimiter(path),
        n_rows=int(df.shape[0]),
        n_cols=int(df.shape[1]),
        columns_sample=cols[:MAX_COLS_TO_SHOW],
        index_name=str(df.index.name) if df.index.name is not None else None,
        first_col=cols[0] if cols else None,
        likely_type=likely_type,
        sample_id_like_cols=find_matching_columns(cols, ["sample", "sample_id", "run", "accession", "id"]),
        sex_like_cols=find_matching_columns(cols, ["sex", "gender", "male", "female"]),
        target_like_cols=find_matching_columns(cols, ["condition", "disease", "status", "class", "label", "crc", "adenoma", "control"]),
        study_like_cols=find_matching_columns(cols, ["study", "cohort", "dataset", "project", "batch", "source"]),
        notes=notes,
    )


def scan_project() -> list[FileSummary]:
    summaries: list[FileSummary] = []

    for dirname in SEARCH_DIRS:
        root = PROJECT_ROOT / dirname
        if not root.exists():
            continue

        for path in root.rglob("*"):
            if path.is_file() and path.suffix.lower() in TEXT_EXTENSIONS:
                summary = summarize_file(path)
                if summary is not None:
                    summaries.append(summary)

    return summaries


def print_summary_table(summaries: list[FileSummary]) -> None:
    if not summaries:
        print("No readable text tables found.")
        return

    print("\n=== DATASET INVENTORY ===\n")
    for s in summaries:
        print(f"[{s.likely_type.upper()}] {s.path}")
        print(f"  shape (preview read): {s.n_rows} rows x {s.n_cols} cols")
        print(f"  delimiter: {repr(s.delimiter)}")
        print(f"  sample columns: {s.columns_sample[:10]}")
        if s.sample_id_like_cols:
            print(f"  sample-id-like cols: {s.sample_id_like_cols}")
        if s.sex_like_cols:
            print(f"  sex-like cols: {s.sex_like_cols}")
        if s.target_like_cols:
            print(f"  target-like cols: {s.target_like_cols}")
        if s.study_like_cols:
            print(f"  study-like cols: {s.study_like_cols}")
        if s.notes:
            print(f"  notes: {' | '.join(s.notes)}")
        print()


def choose_candidates(summaries: list[FileSummary]) -> tuple[list[FileSummary], list[FileSummary]]:
    metadata_candidates = [s for s in summaries if s.likely_type == "metadata"]
    abundance_candidates = [s for s in summaries if s.likely_type == "abundance_matrix"]

    metadata_candidates.sort(
        key=lambda s: (
            len(s.sex_like_cols) + len(s.target_like_cols) + len(s.study_like_cols),
            s.n_cols or 0,
            s.n_rows or 0,
        ),
        reverse=True,
    )

    abundance_candidates.sort(
        key=lambda s: (
            s.n_cols or 0,
            s.n_rows or 0,
            "harmonized" in s.path.lower(),
            "species" in s.path.lower(),
            "abundance" in s.path.lower(),
        ),
        reverse=True,
    )

    return metadata_candidates, abundance_candidates


def suggest_next_step(metadata_candidates: list[FileSummary], abundance_candidates: list[FileSummary]) -> None:
    print("\n=== BEST GUESSES ===\n")

    if abundance_candidates:
        top_ab = abundance_candidates[0]
        print(f"Main abundance matrix candidate: {top_ab.path}")
        print(f"  shape: {top_ab.n_rows} x {top_ab.n_cols}")
    else:
        print("No abundance matrix candidate found.")

    if metadata_candidates:
        top_md = metadata_candidates[0]
        print(f"Main metadata file candidate: {top_md.path}")
        print(f"  shape: {top_md.n_rows} x {top_md.n_cols}")
        print(f"  sex column candidates: {top_md.sex_like_cols}")
        print(f"  target column candidates: {top_md.target_like_cols}")
        print(f"  study column candidates: {top_md.study_like_cols}")
    else:
        print("No metadata candidate found.")

    if abundance_candidates and metadata_candidates:
        print("\nLikely analysis setup:")
        print("- abundance table = numeric taxa features")
        print("- metadata table = labels / covariates / LOSO groups")
        print("- next check = confirm the join key (sample ID)")
        print("- then decide genus vs species level")
    print()


def inspect_top_files(metadata_candidates: list[FileSummary], abundance_candidates: list[FileSummary]) -> None:
    top_files = []
    if abundance_candidates:
        top_files.append(abundance_candidates[0].path)
    if metadata_candidates:
        top_files.append(metadata_candidates[0].path)

    if not top_files:
        return

    print("\n=== HEAD OF TOP CANDIDATES ===\n")
    for rel_path in top_files:
        path = PROJECT_ROOT / rel_path
        df = safe_read_table(path, nrows=5)
        if df is None:
            continue

        print(f"--- {rel_path} ---")
        print(df.head())
        print()


def save_inventory(summaries: list[FileSummary]) -> None:
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump([asdict(s) for s in summaries], f, indent=2, ensure_ascii=False)
    print(f"Saved inventory to {OUTPUT_JSON}")


def main() -> None:
    print(f"Project root: {PROJECT_ROOT}\n")
    summaries = scan_project()
    print_summary_table(summaries)

    metadata_candidates, abundance_candidates = choose_candidates(summaries)
    suggest_next_step(metadata_candidates, abundance_candidates)
    inspect_top_files(metadata_candidates, abundance_candidates)
    save_inventory(summaries)


if __name__ == "__main__":
    main()
