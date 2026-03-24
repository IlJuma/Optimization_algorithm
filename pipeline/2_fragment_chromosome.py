import random
from typing import List, Tuple, Dict

# =========================
# PARAMETERS
# =========================

INPUT_FASTA = "data/fasta/simulated_chromosome.fasta"         # Output from script #1
OUTPUT_FASTA = "data/fasta/fragments.fasta"
OUTPUT_REPORT = "reports/fragments_report.txt"

SEED = 123

# Number of independent genome molecules to fragment
N_GENOME_COPIES = 20

# Random cut-site model
MEAN_CUT_SPACING = 350

# Size selection after fragmentation
MIN_INSERT_SIZE = 250
MAX_INSERT_SIZE = 450

# Recovery / sampling after size selection
RECOVERY_FRACTION = 0.25

# Optional additional sampling cap
MAX_FRAGMENTS_TO_KEEP = None

FASTA_LINE_WIDTH = 80


# =========================
# FASTA I/O
# =========================

def read_fasta(path: str) -> Tuple[str, str]:
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
        header = "sequence"

    seq = "".join(seq_parts).upper()
    return header, seq


def write_fasta(fragments: List[Dict], path: str) -> None:
    with open(path, "w") as f:
        for frag in fragments:
            header = (
                f">{frag['fragment_id']} "
                f"source={frag['source_header']} "
                f"copy={frag['copy_index']} "
                f"frag_start={frag['start']} "
                f"frag_end={frag['end']} "
                f"frag_len={frag['length']}"
            )
            f.write(header + "\n")

            seq = frag["sequence"]
            for j in range(0, len(seq), FASTA_LINE_WIDTH):
                f.write(seq[j:j + FASTA_LINE_WIDTH] + "\n")


# =========================
# FRAGMENTATION
# =========================

def generate_cut_positions(genome_length: int, mean_cut_spacing: float, rng: random.Random) -> List[int]:
    if genome_length <= 0:
        return [0]

    if mean_cut_spacing <= 0:
        raise ValueError("MEAN_CUT_SPACING must be > 0")

    cuts = [0]
    pos = 0.0

    while True:
        step = rng.expovariate(1.0 / mean_cut_spacing)
        pos += step
        cut = int(pos)

        if cut >= genome_length:
            break
        if cut > cuts[-1]:
            cuts.append(cut)

    if cuts[-1] != genome_length:
        cuts.append(genome_length)

    return cuts


def fragment_one_copy(
    seq: str,
    copy_index: int,
    source_header: str,
    mean_cut_spacing: float,
    rng: random.Random,
) -> List[Dict]:
    cuts = generate_cut_positions(len(seq), mean_cut_spacing, rng)
    fragments = []

    for i in range(len(cuts) - 1):
        start = cuts[i]
        end = cuts[i + 1]
        if end <= start:
            continue

        fragments.append({
            "source_header": source_header,
            "copy_index": copy_index,
            "start": start,
            "end": end,
            "length": end - start,
            "sequence": seq[start:end],
        })

    return fragments


def size_select_fragments(
    fragments: List[Dict],
    min_insert_size: int,
    max_insert_size: int,
) -> List[Dict]:
    return [
        frag for frag in fragments
        if min_insert_size <= frag["length"] <= max_insert_size
    ]


def recover_fragments(
    fragments: List[Dict],
    recovery_fraction: float,
    rng: random.Random,
    max_fragments_to_keep: int = None,
) -> List[Dict]:
    if not (0.0 <= recovery_fraction <= 1.0):
        raise ValueError("RECOVERY_FRACTION must be between 0 and 1")

    kept = [frag for frag in fragments if rng.random() < recovery_fraction]

    if max_fragments_to_keep is not None and len(kept) > max_fragments_to_keep:
        kept = rng.sample(kept, max_fragments_to_keep)

    return kept


def assign_fragment_ids(fragments: List[Dict], source_header: str) -> List[Dict]:
    for i, frag in enumerate(fragments, start=1):
        frag["fragment_id"] = f"{source_header}_frag_{i}"
    return fragments


def simulate_library(
    genome_header: str,
    genome_seq: str,
    n_genome_copies: int,
    mean_cut_spacing: float,
    min_insert_size: int,
    max_insert_size: int,
    recovery_fraction: float,
    rng: random.Random,
    max_fragments_to_keep: int = None,
) -> Tuple[List[Dict], dict]:
    all_generated = []
    all_selected = []

    for copy_index in range(1, n_genome_copies + 1):
        generated = fragment_one_copy(
            seq=genome_seq,
            copy_index=copy_index,
            source_header=genome_header,
            mean_cut_spacing=mean_cut_spacing,
            rng=rng,
        )
        selected = size_select_fragments(generated, min_insert_size, max_insert_size)

        all_generated.extend(generated)
        all_selected.extend(selected)

    recovered = recover_fragments(
        fragments=all_selected,
        recovery_fraction=recovery_fraction,
        rng=rng,
        max_fragments_to_keep=max_fragments_to_keep,
    )

    recovered = assign_fragment_ids(recovered, genome_header)

    stats = {
        "genome_length": len(genome_seq),
        "n_genome_copies": n_genome_copies,
        "generated_fragments": len(all_generated),
        "size_selected_fragments": len(all_selected),
        "recovered_fragments": len(recovered),
    }

    return recovered, stats


# =========================
# REPORTING
# =========================

def mean(values: List[int]) -> float:
    return sum(values) / len(values) if values else 0.0


def length_summary(fragments: List[Dict]) -> dict:
    lengths = [f["length"] for f in fragments]
    if not lengths:
        return {"count": 0, "min": 0, "max": 0, "mean": 0.0}

    return {
        "count": len(lengths),
        "min": min(lengths),
        "max": max(lengths),
        "mean": mean(lengths),
    }


def write_report(path: str, stats: dict, fragments: List[Dict]) -> None:
    summary = length_summary(fragments)

    approx_bases = sum(f["length"] for f in fragments)
    approximate_mean_coverage = (
        approx_bases / stats["genome_length"]
        if stats["genome_length"] > 0 else 0.0
    )

    with open(path, "w") as f:
        f.write("Fragmentation simulation report\n")
        f.write("==============================\n\n")
        f.write(f"Genome length: {stats['genome_length']}\n")
        f.write(f"Genome copies fragmented: {stats['n_genome_copies']}\n")
        f.write(f"Fragments generated before size selection: {stats['generated_fragments']}\n")
        f.write(f"Fragments passing size selection: {stats['size_selected_fragments']}\n")
        f.write(f"Fragments recovered after random sampling: {stats['recovered_fragments']}\n\n")
        f.write("Recovered fragment length summary\n")
        f.write("--------------------------------\n")
        f.write(f"Count: {summary['count']}\n")
        f.write(f"Min: {summary['min']}\n")
        f.write(f"Max: {summary['max']}\n")
        f.write(f"Mean: {summary['mean']:.2f}\n\n")
        f.write(f"Total recovered bases: {approx_bases}\n")
        f.write(f"Approximate mean physical coverage: {approximate_mean_coverage:.2f}x\n")


# =========================
# MAIN
# =========================

def main() -> None:
    rng = random.Random(SEED)

    genome_header, genome_seq = read_fasta(INPUT_FASTA)

    fragments, stats = simulate_library(
        genome_header=genome_header,
        genome_seq=genome_seq,
        n_genome_copies=N_GENOME_COPIES,
        mean_cut_spacing=MEAN_CUT_SPACING,
        min_insert_size=MIN_INSERT_SIZE,
        max_insert_size=MAX_INSERT_SIZE,
        recovery_fraction=RECOVERY_FRACTION,
        rng=rng,
        max_fragments_to_keep=MAX_FRAGMENTS_TO_KEEP,
    )

    write_fasta(fragments, OUTPUT_FASTA)
    write_report(OUTPUT_REPORT, stats, fragments)

    print(f"Input FASTA: {INPUT_FASTA}")
    print(f"Genome length: {len(genome_seq)} bp")
    print(f"Fragments generated: {stats['generated_fragments']}")
    print(f"Fragments size-selected: {stats['size_selected_fragments']}")
    print(f"Fragments recovered: {stats['recovered_fragments']}")
    print(f"Output FASTA: {OUTPUT_FASTA}")
    print(f"Report: {OUTPUT_REPORT}")


if __name__ == "__main__":
    main()