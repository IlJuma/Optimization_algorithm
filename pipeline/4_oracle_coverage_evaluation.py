import os
from typing import Dict, List, Tuple, Iterable
import matplotlib.pyplot as plt

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


def parse_key_value_header(header: str) -> Dict[str, str]:
    parts = header.split()
    data = {"record_id": parts[0]}
    for token in parts[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            data[key] = value
    return data


def read_fragments_fasta(path: str) -> List[Dict]:
    records = []
    header = None
    seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    meta = parse_key_value_header(header)
                    meta["sequence"] = "".join(seq_parts).upper()
                    meta["frag_start"] = int(meta["frag_start"])
                    meta["frag_end"] = int(meta["frag_end"])
                    meta["frag_len"] = int(meta["frag_len"])
                    meta["copy"] = int(meta["copy"])
                    records.append(meta)

                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        meta = parse_key_value_header(header)
        meta["sequence"] = "".join(seq_parts).upper()
        meta["frag_start"] = int(meta["frag_start"])
        meta["frag_end"] = int(meta["frag_end"])
        meta["frag_len"] = int(meta["frag_len"])
        meta["copy"] = int(meta["copy"])
        records.append(meta)

    return records


def read_fastq_records(path: str) -> Iterable[Tuple[str, str, str]]:
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

            yield header[1:], seq, qual


def read_reads_fastq(path: str) -> List[Dict]:
    records = []

    for header, seq, qual in read_fastq_records(path):
        meta = parse_key_value_header(header)
        meta["sequence"] = seq
        meta["quality"] = qual
        meta["copy"] = int(meta["copy"])
        meta["frag_start"] = int(meta["frag_start"])
        meta["frag_end"] = int(meta["frag_end"])
        meta["read_start"] = int(meta["read_start"])
        meta["read_end"] = int(meta["read_end"])
        meta["mate"] = int(meta["mate"])
        records.append(meta)

    return records

# =========================
# COVERAGE / INTERVAL LOGIC
# =========================

def increment_coverage(coverage: List[int], start: int, end: int) -> None:
    start = max(0, start)
    end = min(len(coverage), end)

    for i in range(start, end):
        coverage[i] += 1


def intervals_from_fragments(fragments: List[Dict]) -> List[Tuple[int, int, str]]:
    intervals = []
    for frag in fragments:
        intervals.append((frag["frag_start"], frag["frag_end"], frag["record_id"]))
    return intervals


def intervals_from_reads(
    reads_r1: List[Dict],
    reads_r2: List[Dict],
    use_r1: bool,
    use_r2: bool,
) -> List[Tuple[int, int, str]]:
    intervals = []

    if use_r1:
        for read in reads_r1:
            intervals.append((read["read_start"], read["read_end"], read["record_id"]))

    if use_r2:
        for read in reads_r2:
            intervals.append((read["read_start"], read["read_end"], read["record_id"]))

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


def make_binned_coverage(coverage: List[int], bin_size: int) -> Tuple[List[float], List[float]]:
    if bin_size <= 0:
        raise ValueError("PLOT_BIN_SIZE must be > 0")

    x = []
    y = []

    for start in range(0, len(coverage), bin_size):
        end = min(start + bin_size, len(coverage))
        window = coverage[start:end]
        x.append((start + end - 1) / 2.0)
        y.append(sum(window) / len(window))

    return x, y


def write_coverage_plot(path: str, coverage: List[int], target_bins: int, title: str) -> None:
    ensure_reports_dir(path)

    bin_size = choose_plot_bin_size(len(coverage), target_bins)
    x, y = make_binned_coverage(coverage, bin_size)

    plt.figure(figsize=(PLOT_WIDTH, PLOT_HEIGHT))
    plt.plot(x, y)
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

    fragments = read_fragments_fasta(FRAGMENTS_FASTA)
    fragment_intervals = intervals_from_fragments(fragments)

    reads_r1 = read_reads_fastq(READS_R1_FASTQ) if USE_R1 else []
    reads_r2 = read_reads_fastq(READS_R2_FASTQ) if USE_R2 else []
    read_intervals = intervals_from_reads(reads_r1, reads_r2, USE_R1, USE_R2)

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