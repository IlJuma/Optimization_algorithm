"""
Read data loader.

Parses FASTQ files (R1/R2) into structured ReadRecord objects, reconstructing
pair relationships and exposing sequencing metadata such as quality scores and
ground-truth genomic positions. Provides convenient mappings (pairing, lengths,
positions) and DataFrame export for downstream analysis.

Outputs
-------
read_data : LoadedReadData
    Container with:

    reads : List[ReadRecord]
        Each record contains:
            - read_id : str
            - pair_id : str
            - mate : int (1 or 2)
            - sequence : str
            - quality : str
            - read_length : int
            - fragment_id : Optional[str]
            - source : Optional[str]
            - copy : Optional[int]
            - frag_start : Optional[int]
            - frag_end : Optional[int]
            - read_start : Optional[int]
            - read_end : Optional[int]

    read_id_to_pair_id : Dict[str, str]
    pair_id_to_read_ids : Dict[str, List[str]]
    read_lengths : Dict[str, int]
    quality_scores : Dict[str, str]
    ground_truth_positions : Dict[str, Tuple[int, int]]
    copy_indices : Dict[str, Optional[int]]

reads_df : pd.DataFrame
    Columns:
        ["read_id", "pair_id", "mate",
         "sequence", "quality", "read_length",
         "fragment_id", "source", "copy",
         "frag_start", "frag_end",
         "read_start", "read_end"]

Typical variable names
----------------------
read_data
reads
reads_df
pair_id_to_read_ids
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import os

import pandas as pd


# =========================
# DEFAULT INPUT PATHS
# =========================

DEFAULT_READS_R1_FASTQ = "data/fastq/reads_R1.fastq"
DEFAULT_READS_R2_FASTQ = "data/fastq/reads_R2.fastq"


# =========================
# DATA CLASSES
# =========================

@dataclass
class ReadRecord:
    read_id: str
    pair_id: str
    mate: int
    sequence: str
    quality: str
    read_length: int

    fragment_id: Optional[str] = None
    source: Optional[str] = None
    copy: Optional[int] = None
    frag_start: Optional[int] = None
    frag_end: Optional[int] = None
    read_start: Optional[int] = None
    read_end: Optional[int] = None


@dataclass
class LoadedReadData:
    reads: List[ReadRecord]
    read_id_to_pair_id: Dict[str, str]
    pair_id_to_read_ids: Dict[str, List[str]]
    read_lengths: Dict[str, int]
    quality_scores: Dict[str, str]
    ground_truth_positions: Dict[str, Tuple[Optional[int], Optional[int]]]
    copy_indices: Dict[str, Optional[int]]


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


def read_fastq_records(path: str) -> List[Tuple[str, str, str]]:
    records: List[Tuple[str, str, str]] = []

    with open(path) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break

            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header in {path}: {header}")
            if plus != "+":
                raise ValueError(f"Invalid FASTQ record in {path}: expected '+' line")

            records.append((header[1:], seq, qual))

    return records


# =========================
# READ LOADING
# =========================

def infer_pair_id_and_mate(record_id: str) -> Tuple[str, Optional[int]]:
    """
    Expected record_id shape from your simulator:
    simread:<index>:<mate>
    Example:
    simread:18:1
    """
    parts = record_id.split(":")
    if len(parts) >= 3:
        pair_id = ":".join(parts[:2])
        try:
            mate = int(parts[2])
        except ValueError:
            mate = None
        return pair_id, mate
    return record_id, None


def load_reads_from_fastq(path: str) -> List[ReadRecord]:
    raw_records = read_fastq_records(path)
    reads: List[ReadRecord] = []

    for header, seq, qual in raw_records:
        meta = parse_key_value_header(header)
        pair_id, mate_from_id = infer_pair_id_and_mate(meta["record_id"])

        mate = int(meta["mate"]) if "mate" in meta else mate_from_id
        copy = int(meta["copy"]) if "copy" in meta else None
        frag_start = int(meta["frag_start"]) if "frag_start" in meta else None
        frag_end = int(meta["frag_end"]) if "frag_end" in meta else None
        read_start = int(meta["read_start"]) if "read_start" in meta else None
        read_end = int(meta["read_end"]) if "read_end" in meta else None

        reads.append(
            ReadRecord(
                read_id=meta["record_id"],
                pair_id=pair_id,
                mate=mate if mate is not None else -1,
                sequence=seq,
                quality=qual,
                read_length=len(seq),
                fragment_id=meta.get("fragment_id"),
                source=meta.get("source"),
                copy=copy,
                frag_start=frag_start,
                frag_end=frag_end,
                read_start=read_start,
                read_end=read_end,
            )
        )

    return reads


def load_reads(
    reads_r1_fastq: str = DEFAULT_READS_R1_FASTQ,
    reads_r2_fastq: Optional[str] = DEFAULT_READS_R2_FASTQ,
) -> LoadedReadData:
    reads: List[ReadRecord] = []

    if reads_r1_fastq and os.path.exists(reads_r1_fastq):
        reads.extend(load_reads_from_fastq(reads_r1_fastq))

    if reads_r2_fastq and os.path.exists(reads_r2_fastq):
        reads.extend(load_reads_from_fastq(reads_r2_fastq))

    read_id_to_pair_id: Dict[str, str] = {}
    pair_id_to_read_ids: Dict[str, List[str]] = {}
    read_lengths: Dict[str, int] = {}
    quality_scores: Dict[str, str] = {}
    ground_truth_positions: Dict[str, Tuple[Optional[int], Optional[int]]] = {}
    copy_indices: Dict[str, Optional[int]] = {}

    for read in reads:
        read_id_to_pair_id[read.read_id] = read.pair_id
        pair_id_to_read_ids.setdefault(read.pair_id, []).append(read.read_id)
        read_lengths[read.read_id] = read.read_length
        quality_scores[read.read_id] = read.quality
        ground_truth_positions[read.read_id] = (read.read_start, read.read_end)
        copy_indices[read.read_id] = read.copy

    return LoadedReadData(
        reads=reads,
        read_id_to_pair_id=read_id_to_pair_id,
        pair_id_to_read_ids=pair_id_to_read_ids,
        read_lengths=read_lengths,
        quality_scores=quality_scores,
        ground_truth_positions=ground_truth_positions,
        copy_indices=copy_indices,
    )


def reads_to_dataframe(reads: List[ReadRecord]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "read_id": r.read_id,
                "pair_id": r.pair_id,
                "mate": r.mate,
                "sequence": r.sequence,
                "quality": r.quality,
                "read_length": r.read_length,
                "fragment_id": r.fragment_id,
                "source": r.source,
                "copy": r.copy,
                "frag_start": r.frag_start,
                "frag_end": r.frag_end,
                "read_start": r.read_start,
                "read_end": r.read_end,
            }
            for r in reads
        ]
    )


def export_reads_dataframe(
    reads_r1_fastq: str = DEFAULT_READS_R1_FASTQ,
    reads_r2_fastq: Optional[str] = DEFAULT_READS_R2_FASTQ,
) -> pd.DataFrame:
    read_data = load_reads(reads_r1_fastq, reads_r2_fastq)
    return reads_to_dataframe(read_data.reads)


# =========================
# EXAMPLE USAGE
# =========================

if __name__ == "__main__":
    read_data = load_reads()

    print(f"Loaded reads: {len(read_data.reads)}")
    print(f"Read pairs: {len(read_data.pair_id_to_read_ids)}")

    reads_df = reads_to_dataframe(read_data.reads)

    print()
    print("Reads dataframe columns:")
    print(list(reads_df.columns))