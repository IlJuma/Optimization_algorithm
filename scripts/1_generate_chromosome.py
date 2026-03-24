import random
import os
from typing import List, Optional

# =========================
# PARAMETERS (EDIT HERE)
# =========================

LENGTH = 1_000_000                  # Total sequence length
SEED = 42                           # Random seed for reproducibility

# Global / local GC control
MEAN_GC_CONTENT = 0.50              # Mean GC across the whole sequence
GC_WINDOW_SIZE = 10_000             # GC is assigned per window
GC_STDDEV = 0.08                    # Variation in GC between windows
MIN_GC = 0.20                       # Clamp local GC to this minimum
MAX_GC = 0.80                       # Clamp local GC to this maximum

# Background sequence features
MAX_HOMOPOLYMER = 10                # Max allowed homopolymer length (None to disable)

# Output
OUTPUT_FASTA = "fasta/simulated.fasta"
FASTA_HEADER = "simulated_sequence"
OUTPUT_REPORT = "reports/simulated_report.txt"


# =========================
# HELPER FUNCTIONS
# =========================

def clamp(value: float, minimum: float, maximum: float) -> float:
    return max(minimum, min(maximum, value))


def random_base(gc_content: float, rng: random.Random) -> str:
    at = (1.0 - gc_content) / 2.0
    gc = gc_content / 2.0
    return rng.choices(["A", "C", "G", "T"], weights=[at, gc, gc, at], k=1)[0]


# =========================
# GC WINDOWS
# =========================

def generate_gc_profile(
    length: int,
    window_size: int,
    mean_gc: float,
    stddev: float,
    min_gc: float,
    max_gc: float,
    rng: random.Random,
) -> List[float]:
    gc_profile = []
    pos = 0

    while pos < length:
        local_gc = clamp(rng.gauss(mean_gc, stddev), min_gc, max_gc)
        this_window = min(window_size, length - pos)
        gc_profile.extend([local_gc] * this_window)
        pos += this_window

    return gc_profile


# =========================
# BACKGROUND GENERATION
# =========================

def generate_background_sequence(
    length: int,
    gc_profile: List[float],
    rng: random.Random,
    max_homopolymer: Optional[int],
) -> List[str]:
    seq = []

    for i in range(length):
        base = random_base(gc_profile[i], rng)

        if max_homopolymer and len(seq) >= max_homopolymer:
            run = seq[-max_homopolymer:]
            if all(b == run[0] for b in run) and base == run[0]:
                choices = ["A", "C", "G", "T"]
                choices.remove(run[0])
                base = rng.choice(choices)

        seq.append(base)

    return seq


# =========================
# OUTPUT / REPORTING
# =========================

def compute_gc_fraction(seq: str) -> float:
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq) if seq else 0.0


def write_fasta(seq: str, path: str, header: str, width: int = 80) -> None:
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def write_report(
    report_path: str,
    gc_windows: List[float],
    final_seq: str,
) -> None:
    with open(report_path, "w") as f:
        f.write("Simulation report\n")
        f.write("=================\n\n")
        f.write(f"Final length: {len(final_seq)}\n")
        f.write(f"Final GC fraction: {compute_gc_fraction(final_seq):.4f}\n")
        f.write(f"Number of GC windows: {len(gc_windows)}\n")
        f.write(f"Mean window GC: {sum(gc_windows) / len(gc_windows):.4f}\n")


# =========================
# MAIN
# =========================

def main():
    rng = random.Random(SEED)

    # Create output folders
    os.makedirs("fasta", exist_ok=True)
    os.makedirs("reports", exist_ok=True)

    gc_profile = generate_gc_profile(
        length=LENGTH,
        window_size=GC_WINDOW_SIZE,
        mean_gc=MEAN_GC_CONTENT,
        stddev=GC_STDDEV,
        min_gc=MIN_GC,
        max_gc=MAX_GC,
        rng=rng,
    )

    gc_windows = [
        gc_profile[i]
        for i in range(0, len(gc_profile), GC_WINDOW_SIZE)
    ]

    seq = generate_background_sequence(
        length=LENGTH,
        gc_profile=gc_profile,
        rng=rng,
        max_homopolymer=MAX_HOMOPOLYMER,
    )

    seq_str = "".join(seq)

    write_fasta(seq_str, OUTPUT_FASTA, FASTA_HEADER)
    write_report(OUTPUT_REPORT, gc_windows, seq_str)

    print(f"Sequence generated: {len(seq_str)} bp")
    print(f"Final GC fraction: {compute_gc_fraction(seq_str):.4f}")
    print(f"FASTA written to: {OUTPUT_FASTA}")
    print(f"Report written to: {OUTPUT_REPORT}")


if __name__ == "__main__":
    main()