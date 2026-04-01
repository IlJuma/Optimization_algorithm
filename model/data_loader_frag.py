"""
Fragment data loader.

Parses fragment FASTA files produced by the simulation pipeline into structured
FragmentRecord objects and optional pandas DataFrames. Includes full ground-truth
metadata (positions, copy index) for downstream evaluation, but can also be used
to provide reduced views for optimization algorithms.

Outputs
-------
fragments : List[FragmentRecord]
    Each record contains:
        - fragment_id : str
        - source : str
        - copy : int
        - frag_start : int
        - frag_end : int
        - frag_len : int
        - sequence : str

fragments_df : pd.DataFrame
    Columns:
        ["fragment_id", "source", "copy",
         "frag_start", "frag_end", "frag_len", "sequence"]

Typical variable names
----------------------
fragments
fragments_df
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import pandas as pd


# =========================
# DEFAULT INPUT PATHS
# =========================

DEFAULT_FRAGMENTS_FASTA = "data/fasta/fragments.fasta"


# =========================
# DATA CLASSES
# =========================

@dataclass
class FragmentRecord:
    fragment_id: str
    source: str
    copy: int
    frag_start: int
    frag_end: int
    frag_len: int
    sequence: str


# =========================
# LOW-LEVEL PARSERS
# =========================

def parse_key_value_header(header: str) -> Dict[str, str]:
    parts = header.split()
    data = {"record_id": parts[0]}
    for token in parts[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            data[key] = value
    return data


def read_fasta_records(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = None
    seq_parts: List[str] = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts).upper()))

    return records


# =========================
# FRAGMENT LOADING
# =========================

def load_fragments(fragments_fasta: str = DEFAULT_FRAGMENTS_FASTA) -> List[FragmentRecord]:
    raw_records = read_fasta_records(fragments_fasta)
    fragments: List[FragmentRecord] = []

    for header, seq in raw_records:
        meta = parse_key_value_header(header)

        fragments.append(
            FragmentRecord(
                fragment_id=meta["record_id"],
                source=meta["source"],
                copy=int(meta["copy"]),
                frag_start=int(meta["frag_start"]),
                frag_end=int(meta["frag_end"]),
                frag_len=int(meta["frag_len"]),
                sequence=seq,
            )
        )

    return fragments


def fragments_to_dataframe(fragments: List[FragmentRecord]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "fragment_id": f.fragment_id,
                "source": f.source,
                "copy": f.copy,
                "frag_start": f.frag_start,
                "frag_end": f.frag_end,
                "frag_len": f.frag_len,
                "sequence": f.sequence,
            }
            for f in fragments
        ]
    )


def export_fragments_dataframe(
    fragments_fasta: str = DEFAULT_FRAGMENTS_FASTA,
) -> pd.DataFrame:
    fragments = load_fragments(fragments_fasta)
    return fragments_to_dataframe(fragments)


# =========================
# EXAMPLE USAGE
# =========================

if __name__ == "__main__":
    fragments = load_fragments()

    print(f"Loaded fragments: {len(fragments)}")

    fragments_df = fragments_to_dataframe(fragments)

    print()
    print("Fragments dataframe columns:")
    print(list(fragments_df.columns))