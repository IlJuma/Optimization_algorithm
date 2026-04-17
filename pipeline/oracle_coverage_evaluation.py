"""
Oracle coverage evaluation.

Uses ground-truth metadata from fragments and reads to reconstruct their true
genomic placement and evaluate coverage along the reference genome. Computes
coverage profiles, detects gaps and contigs, and generates reports and plots
for both fragment-level (physical) and read-level (sequenced) coverage.

Fragment orientation note
-------------------------
Stored fragment orientation does not change reference placement for coverage.
Coverage is reconstructed from truth intervals:
    frag_start, frag_end
and for reads:
    read_start, read_end

So fragment orientation metadata can be ignored for interval-based coverage
evaluation.

Inputs
------
fragments : List[FragmentRecord]
reads : List[ReadRecord]
reference sequence (FASTA)

Internal representations
------------------------
fragment_intervals : List[Tuple[int, int, str]]
    (frag_start, frag_end, fragment_id)

read_intervals : List[Tuple[int, int, str]]
    (read_start, read_end, read_id)

coverage arrays
---------------
fragment_coverage : List[int]
read_coverage : List[int]

Outputs
-------
Coverage TSVs:
    - fragments_coverage.tsv
    - reads_coverage.tsv
    Columns:
        position, coverage

Gap intervals:
    List[Tuple[int, int]]  # (start, end)

Contigs:
    List[Tuple[int, int]]  # (start, end)

Plots:
    - fragments_coverage_plot.png
    - reads_coverage_plot.png

Report:
    reconstruction_report.txt

Typical variable names
----------------------
fragment_intervals
read_intervals
fragment_coverage
read_coverage
fragment_gaps
read_gaps
fragment_contigs
read_contigs
"""

import os
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt

from model.data_loader_frag import load_fragments, FragmentRecord
from model.data_loader_read import load_reads, ReadRecord

# =========================
# PARAMETERS (EDIT HERE)
# =========================

# Input layout
REFERENCE_FASTA = "data/fasta/simulated_chromosome.fasta"
FRAGMENTS_FASTA = "data/fasta/fragments.fasta"
READS_R1_FASTQ = "data/fastq/reads_R1.fastq"
READS_R2_FASTQ = "data/fastq/reads_R2.fastq"

# Always evaluate both reads and fragments
USE_R1 = True
USE_R2 = True

# Plotting
FRAGMENTS_COVERAGE_PLOT = "reports/fragments_coverage_plot.png"
READS_COVERAGE_PLOT = "reports/reads_coverage_plot.png"
PLOT_WIDTH = 14
PLOT_HEIGHT = 4
PLOT_DPI = 150

# Optional smoothing/binning for the plot only
PLOT_TARGET_BINS = 250  # aim for about this many plotted bins

# Outputs
OUTPUT_REPORT = "reports/reconstruction_report.txt"

FRAGMENTS_OUTPUT_GAPS = "reports/fragments_gaps.tsv"
FRAGMENTS_OUTPUT_CONTIGS = "reports/fragments_contigs.tsv"
FRAGMENTS_OUTPUT_COVERAGE = "reports/fragments_coverage.tsv"

READS_OUTPUT_GAPS = "reports/reads_gaps.tsv"
READS_OUTPUT_CONTIGS = "reports/reads_contigs.tsv"
READS_OUTPUT_COVERAGE = "reports/reads_coverage.tsv"

# =========================
# I/O HELPERS
# =========================

def ensure_reports_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def read_single_fasta(path: str) -> Tuple[str, str]:
    header = None
    seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:]
            else:
                seq_parts.append(line)

    if header is None:
        raise ValueError(f"No FASTA header found in {path}")

    return header, "".join(seq_parts).upper()


# =========================
# COVERAGE / INTERVAL LOGIC
# =========================

def increment_coverage(coverage: List[int], start: int, end: int) -> None:
    start = max(0, start)
    end = min(len(coverage), end)

    for i in range(start, end):
        coverage[i] += 1


def intervals_from_fragments(fragments: List[FragmentRecord]) -> List[Tuple[int, int, str]]:
    intervals = []
    for frag in fragments:
        intervals.append((frag.frag_start, frag.frag_end, frag.fragment_id))
    return intervals


def intervals_from_reads(
    reads: List[ReadRecord],
    use_r1: bool,
    use_r2: bool,
) -> List[Tuple[int, int, str]]:
    intervals = []

    for read in reads:
        if read.mate == 1 and not use_r1:
            continue
        if read.mate == 2 and not use_r2:
            continue
        if read.read_start is None or read.read_end is None:
            continue

        intervals.append((read.read_start, read.read_end, read.read_id))

    return intervals


def build_coverage(genome_length: int, intervals: List[Tuple[int, int, str]]) -> List[int]:
    coverage = [0] * genome_length
    for start, end, _ in intervals:
        increment_coverage(coverage, start, end)
    return coverage


def find_zero_coverage_gaps(coverage: List[int]) -> List[Tuple[int, int]]:
    gaps = []
    in_gap = False
    gap_start = None

    for i, cov in enumerate(coverage):
        if cov == 0 and not in_gap:
            in_gap = True
            gap_start = i
        elif cov > 0 and in_gap:
            gaps.append((gap_start, i))
            in_gap = False
            gap_start = None

    if in_gap:
        gaps.append((gap_start, len(coverage)))

    return gaps


def find_covered_contigs(coverage: List[int]) -> List[Tuple[int, int]]:
    contigs = []
    in_contig = False
    contig_start = None

    for i, cov in enumerate(coverage):
        if cov > 0 and not in_contig:
            in_contig = True
            contig_start = i
        elif cov == 0 and in_contig:
            contigs.append((contig_start, i))
            in_contig = False
            contig_start = None

    if in_contig:
        contigs.append((contig_start, len(coverage)))

    return contigs


def coverage_stats(coverage: List[int]) -> Dict[str, float]:
    genome_length = len(coverage)
    covered_positions = sum(1 for x in coverage if x > 0)
    uncovered_positions = genome_length - covered_positions
    total_coverage = sum(coverage)

    mean_coverage = total_coverage / genome_length if genome_length > 0 else 0.0
    max_coverage = max(coverage) if coverage else 0
    min_coverage = min(coverage) if coverage else 0
    covered_fraction = covered_positions / genome_length if genome_length > 0 else 0.0

    return {
        "genome_length": genome_length,
        "covered_positions": covered_positions,
        "uncovered_positions": uncovered_positions,
        "covered_fraction": covered_fraction,
        "mean_coverage": mean_coverage,
        "min_coverage": min_coverage,
        "max_coverage": max_coverage,
    }


# =========================
# OUTPUT
# =========================

def write_coverage_tsv(path: str, coverage: List[int]) -> None:
    ensure_reports_dir(path)
    with open(path, "w") as f:
        f.write("position\tcoverage\n")
        for pos, cov in enumerate(coverage):
            f.write(f"{pos}\t{cov}\n")


def write_intervals_tsv(path: str, intervals: List[Tuple[int, int]], kind: str) -> None:
    ensure_reports_dir(path)
    with open(path, "w") as f:
        f.write(f"{kind}_id\tstart\tend\tlength\n")
        for idx, (start, end) in enumerate(intervals, start=1):
            f.write(f"{kind}_{idx}\t{start}\t{end}\t{end - start}\n")


def write_section(
    f,
    label: str,
    interval_count: int,
    stats: Dict[str, float],
    gaps: List[Tuple[int, int]],
    contigs: List[Tuple[int, int]],
) -> None:
    f.write(f"{label}\n")
    f.write(f"{'-' * len(label)}\n")
    f.write(f"Placed intervals used for coverage: {interval_count}\n")
    f.write(f"Covered positions: {stats['covered_positions']}\n")
    f.write(f"Uncovered positions: {stats['uncovered_positions']}\n")
    f.write(f"Covered fraction: {stats['covered_fraction']:.6f}\n")
    f.write(f"Mean coverage: {stats['mean_coverage']:.4f}\n")
    f.write(f"Min coverage: {stats['min_coverage']}\n")
    f.write(f"Max coverage: {stats['max_coverage']}\n")
    f.write(f"Number of covered contigs: {len(contigs)}\n")
    f.write(f"Number of gaps: {len(gaps)}\n")

    if gaps:
        first_gap = gaps[0]
        last_gap = gaps[-1]
        f.write(f"First gap: {first_gap[0]}-{first_gap[1]} ({first_gap[1] - first_gap[0]} bp)\n")
        f.write(f"Last gap: {last_gap[0]}-{last_gap[1]} ({last_gap[1] - last_gap[0]} bp)\n")
    else:
        f.write("No gaps detected.\n")

    if contigs:
        largest = max(contigs, key=lambda x: x[1] - x[0])
        f.write(f"Largest covered contig: {largest[0]}-{largest[1]} ({largest[1] - largest[0]} bp)\n")

    f.write("\n")


def write_report(
    path: str,
    reference_header: str,
    reference_length: int,
    fragment_interval_count: int,
    fragment_stats: Dict[str, float],
    fragment_gaps: List[Tuple[int, int]],
    fragment_contigs: List[Tuple[int, int]],
    read_interval_count: int,
    read_stats: Dict[str, float],
    read_gaps: List[Tuple[int, int]],
    read_contigs: List[Tuple[int, int]],
) -> None:
    ensure_reports_dir(path)
    with open(path, "w") as f:
        f.write("Oracle reconstruction / coverage evaluation report\n")
        f.write("=================================================\n\n")
        f.write(f"Reference: {reference_header}\n")
        f.write(f"Reference length: {reference_length}\n\n")

        write_section(
            f,
            "Fragment-based physical coverage",
            fragment_interval_count,
            fragment_stats,
            fragment_gaps,
            fragment_contigs,
        )

        write_section(
            f,
            "Read-based sequenced coverage",
            read_interval_count,
            read_stats,
            read_gaps,
            read_contigs,
        )


def choose_plot_bin_size(genome_length: int, target_bins: int) -> int:
    if target_bins <= 0:
        raise ValueError("PLOT_TARGET_BINS must be > 0")
    return max(1, genome_length // target_bins)


def make_binned_coverage(
    coverage: List[int],
    bin_size: int,
) -> Tuple[List[float], List[float], List[bool]]:
    if bin_size <= 0:
        raise ValueError("PLOT_BIN_SIZE must be > 0")

    x = []
    y = []
    has_gap = []

    for start in range(0, len(coverage), bin_size):
        end = min(start + bin_size, len(coverage))
        window = coverage[start:end]
        x.append((start + end - 1) / 2.0)
        y.append(sum(window) / len(window))
        has_gap.append(any(v == 0 for v in window))

    return x, y, has_gap


def write_coverage_plot(path: str, coverage: List[int], target_bins: int, title: str) -> None:
    ensure_reports_dir(path)

    bin_size = choose_plot_bin_size(len(coverage), target_bins)
    x, y, has_gap = make_binned_coverage(coverage, bin_size)

    colors = ["red" if gap else "blue" for gap in has_gap]

    plt.figure(figsize=(PLOT_WIDTH, PLOT_HEIGHT))

    for i in range(len(x) - 1):
        if has_gap[i] and has_gap[i + 1]:
            color = "red"
        else:
            color = "blue"

        plt.plot(
            [x[i], x[i + 1]],
            [y[i], y[i + 1]],
            color=color,
        )

    if len(x) == 1:
        plt.plot(x, y, color=colors[0])

    plt.xlabel("Genome position (bp)")
    plt.ylabel(f"Mean coverage per {bin_size} bp bin")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=PLOT_DPI)
    plt.close()


# =========================
# MAIN
# =========================

def main() -> None:
    reference_header, reference_seq = read_single_fasta(REFERENCE_FASTA)
    genome_length = len(reference_seq)

    fragments = load_fragments(FRAGMENTS_FASTA)
    fragment_intervals = intervals_from_fragments(fragments)

    read_data = load_reads(
        reads_r1_fastq=READS_R1_FASTQ if USE_R1 else None,
        reads_r2_fastq=READS_R2_FASTQ if USE_R2 else None,
    )
    read_intervals = intervals_from_reads(read_data.reads, USE_R1, USE_R2)

    fragment_coverage = build_coverage(genome_length, fragment_intervals)
    read_coverage = build_coverage(genome_length, read_intervals)

    fragment_gaps = find_zero_coverage_gaps(fragment_coverage)
    fragment_contigs = find_covered_contigs(fragment_coverage)
    fragment_stats = coverage_stats(fragment_coverage)

    read_gaps = find_zero_coverage_gaps(read_coverage)
    read_contigs = find_covered_contigs(read_coverage)
    read_stats = coverage_stats(read_coverage)

    write_coverage_tsv(FRAGMENTS_OUTPUT_COVERAGE, fragment_coverage)
    write_intervals_tsv(FRAGMENTS_OUTPUT_GAPS, fragment_gaps, "gap")
    write_intervals_tsv(FRAGMENTS_OUTPUT_CONTIGS, fragment_contigs, "contig")

    write_coverage_tsv(READS_OUTPUT_COVERAGE, read_coverage)
    write_intervals_tsv(READS_OUTPUT_GAPS, read_gaps, "gap")
    write_intervals_tsv(READS_OUTPUT_CONTIGS, read_contigs, "contig")

    write_report(
        OUTPUT_REPORT,
        reference_header=reference_header,
        reference_length=genome_length,
        fragment_interval_count=len(fragment_intervals),
        fragment_stats=fragment_stats,
        fragment_gaps=fragment_gaps,
        fragment_contigs=fragment_contigs,
        read_interval_count=len(read_intervals),
        read_stats=read_stats,
        read_gaps=read_gaps,
        read_contigs=read_contigs,
    )

    write_coverage_plot(
        FRAGMENTS_COVERAGE_PLOT,
        fragment_coverage,
        PLOT_TARGET_BINS,
        "Fragment-based physical coverage along genome",
    )
    write_coverage_plot(
        READS_COVERAGE_PLOT,
        read_coverage,
        PLOT_TARGET_BINS,
        "Read-based sequenced coverage along genome",
    )

    print(f"Reference: {reference_header}")
    print(f"Genome length: {genome_length} bp")
    print()
    print("Fragments")
    print(f"  Intervals placed: {len(fragment_intervals)}")
    print(f"  Covered positions: {fragment_stats['covered_positions']}")
    print(f"  Uncovered positions: {fragment_stats['uncovered_positions']}")
    print(f"  Covered contigs: {len(fragment_contigs)}")
    print(f"  Gaps: {len(fragment_gaps)}")
    print(f"  Coverage TSV: {FRAGMENTS_OUTPUT_COVERAGE}")
    print(f"  Gaps TSV: {FRAGMENTS_OUTPUT_GAPS}")
    print(f"  Contigs TSV: {FRAGMENTS_OUTPUT_CONTIGS}")
    print(f"  Coverage plot: {FRAGMENTS_COVERAGE_PLOT}")
    print()
    print("Reads")
    print(f"  Intervals placed: {len(read_intervals)}")
    print(f"  Covered positions: {read_stats['covered_positions']}")
    print(f"  Uncovered positions: {read_stats['uncovered_positions']}")
    print(f"  Covered contigs: {len(read_contigs)}")
    print(f"  Gaps: {len(read_gaps)}")
    print(f"  Coverage TSV: {READS_OUTPUT_COVERAGE}")
    print(f"  Gaps TSV: {READS_OUTPUT_GAPS}")
    print(f"  Contigs TSV: {READS_OUTPUT_CONTIGS}")
    print(f"  Coverage plot: {READS_COVERAGE_PLOT}")
    print()
    print(f"Combined report: {OUTPUT_REPORT}")


if __name__ == "__main__":
    main()