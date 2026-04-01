import random
from typing import List, Tuple, Dict

# =========================
# PARAMETERS
# =========================

INPUT_FASTA = "data/fasta/fragments.fasta"          # Output from fragment_chromosome.py
OUTPUT_FASTQ_R1 = "data/fastq/reads_R1.fastq"
OUTPUT_FASTQ_R2 = "data/fastq/reads_R2.fastq"
OUTPUT_REPORT = "reports/sequencing_report.txt"

SEED = 456

PAIRED_END = True
READ_LENGTH = 150

FRAGMENT_SAMPLING_FRACTION = 1.0
MAX_FRAGMENTS_TO_SEQUENCE = None

START_QUALITY = 36
END_QUALITY = 28
QUALITY_JITTER_STDDEV = 2.0
MIN_QUALITY = 2
MAX_QUALITY = 40

READ_NAME_PREFIX = "simread"


# =========================
# FASTA I/O
# =========================

def parse_header_metadata(header: str) -> Dict[str, str]:
    parts = header.split()
    data = {"record_id": parts[0]}

    for token in parts[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            data[key] = value

    return data


def read_fasta_records(path: str) -> List[Dict]:
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
                    meta = parse_header_metadata(header)
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
        meta = parse_header_metadata(header)
        meta["sequence"] = "".join(seq_parts).upper()
        meta["frag_start"] = int(meta["frag_start"])
        meta["frag_end"] = int(meta["frag_end"])
        meta["frag_len"] = int(meta["frag_len"])
        meta["copy"] = int(meta["copy"])
        records.append(meta)

    return records


# =========================
# BASIC DNA UTILITIES
# =========================

def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def random_base_except(base: str, rng: random.Random) -> str:
    choices = ["A", "C", "G", "T"]
    choices.remove(base)
    return rng.choice(choices)


# =========================
# QUALITY PROFILE
# =========================

def clamp(value: int, minimum: int, maximum: int) -> int:
    return max(minimum, min(maximum, value))


def generate_quality_profile(
    read_length: int,
    start_quality: int,
    end_quality: int,
    jitter_stddev: float,
    min_quality: int,
    max_quality: int,
    rng: random.Random,
) -> List[int]:
    if read_length <= 0:
        return []

    if read_length == 1:
        means = [start_quality]
    else:
        means = [
            start_quality + (end_quality - start_quality) * i / (read_length - 1)
            for i in range(read_length)
        ]

    qs = []
    for q in means:
        if jitter_stddev > 0:
            q = int(round(rng.gauss(q, jitter_stddev)))
        else:
            q = int(round(q))
        qs.append(clamp(q, min_quality, max_quality))

    return qs


def phred_to_ascii(q: int) -> str:
    return chr(q + 33)


# =========================
# ERROR SIMULATION
# =========================

def phred_error_probability(q: int) -> float:
    return 10 ** (-q / 10.0)


def apply_sequencing_errors(
    template_seq: str,
    quality_scores: List[int],
    rng: random.Random,
) -> Tuple[str, str, int]:
    read_bases = []
    qual_chars = []
    error_count = 0

    for base, q in zip(template_seq, quality_scores):
        p_error = phred_error_probability(q)

        if rng.random() < p_error:
            base = random_base_except(base, rng)
            error_count += 1

        read_bases.append(base)
        qual_chars.append(phred_to_ascii(q))

    return "".join(read_bases), "".join(qual_chars), error_count


# =========================
# READ SELECTION
# =========================

def choose_fragments_to_sequence(
    fragments: List[Dict],
    sampling_fraction: float,
    max_fragments: int,
    rng: random.Random,
) -> List[Dict]:
    if not (0.0 <= sampling_fraction <= 1.0):
        raise ValueError("FRAGMENT_SAMPLING_FRACTION must be between 0 and 1")

    chosen = [rec for rec in fragments if rng.random() < sampling_fraction]

    if max_fragments is not None and len(chosen) > max_fragments:
        chosen = rng.sample(chosen, max_fragments)

    return chosen


# =========================
# READ GENERATION
# =========================

def make_read_header(
    read_index: int,
    mate: int,
    fragment: Dict,
    read_start: int,
    read_end: int,
) -> str:
    return (
        f"@{READ_NAME_PREFIX}:{read_index}:{mate} "
        f"fragment_id={fragment['record_id']} "
        f"source={fragment['source']} "
        f"copy={fragment['copy']} "
        f"frag_start={fragment['frag_start']} "
        f"frag_end={fragment['frag_end']} "
        f"read_start={read_start} "
        f"read_end={read_end} "
        f"mate={mate}"
    )


def make_single_end_read(
    fragment: Dict,
    read_length: int,
    rng: random.Random,
    read_index: int,
) -> Tuple[str, str, str, int]:
    template = fragment["sequence"][:read_length]
    read_start = fragment["frag_start"]
    read_end = fragment["frag_start"] + len(template)

    qualities = generate_quality_profile(
        read_length=len(template),
        start_quality=START_QUALITY,
        end_quality=END_QUALITY,
        jitter_stddev=QUALITY_JITTER_STDDEV,
        min_quality=MIN_QUALITY,
        max_quality=MAX_QUALITY,
        rng=rng,
    )

    read_seq, qual_str, errors = apply_sequencing_errors(template, qualities, rng)
    header = make_read_header(read_index, 1, fragment, read_start, read_end)

    return header, read_seq, qual_str, errors


def make_paired_end_reads(
    fragment: Dict,
    read_length: int,
    rng: random.Random,
    read_index: int,
) -> Tuple[Tuple[str, str, str, int], Tuple[str, str, str, int]]:
    seq = fragment["sequence"]
    frag_start = fragment["frag_start"]
    frag_end = fragment["frag_end"]

    r1_template = seq[:read_length]
    r2_template = reverse_complement(seq[-read_length:])

    r1_start = frag_start
    r1_end = frag_start + len(r1_template)

    r2_start = frag_end - len(r2_template)
    r2_end = frag_end

    q1 = generate_quality_profile(
        read_length=len(r1_template),
        start_quality=START_QUALITY,
        end_quality=END_QUALITY,
        jitter_stddev=QUALITY_JITTER_STDDEV,
        min_quality=MIN_QUALITY,
        max_quality=MAX_QUALITY,
        rng=rng,
    )
    q2 = generate_quality_profile(
        read_length=len(r2_template),
        start_quality=START_QUALITY,
        end_quality=END_QUALITY,
        jitter_stddev=QUALITY_JITTER_STDDEV,
        min_quality=MIN_QUALITY,
        max_quality=MAX_QUALITY,
        rng=rng,
    )

    r1_seq, r1_qual, e1 = apply_sequencing_errors(r1_template, q1, rng)
    r2_seq, r2_qual, e2 = apply_sequencing_errors(r2_template, q2, rng)

    h1 = make_read_header(read_index, 1, fragment, r1_start, r1_end)
    h2 = make_read_header(read_index, 2, fragment, r2_start, r2_end)

    return (h1, r1_seq, r1_qual, e1), (h2, r2_seq, r2_qual, e2)


# =========================
# FASTQ OUTPUT
# =========================

def write_fastq(records: List[Tuple[str, str, str]], path: str) -> None:
    with open(path, "w") as f:
        for header, seq, qual in records:
            f.write(header + "\n")
            f.write(seq + "\n")
            f.write("+\n")
            f.write(qual + "\n")


# =========================
# REPORTING
# =========================

def write_report(
    path: str,
    total_fragments_input: int,
    total_fragments_selected: int,
    total_reads_r1: int,
    total_reads_r2: int,
    total_bases_r1: int,
    total_bases_r2: int,
    total_errors_r1: int,
    total_errors_r2: int,
) -> None:
    error_rate_r1 = total_errors_r1 / total_bases_r1 if total_bases_r1 else 0.0
    error_rate_r2 = total_errors_r2 / total_bases_r2 if total_bases_r2 else 0.0

    with open(path, "w") as f:
        f.write("Sequencing simulation report\n")
        f.write("===========================\n\n")
        f.write(f"Input fragments available: {total_fragments_input}\n")
        f.write(f"Fragments selected for sequencing: {total_fragments_selected}\n\n")
        f.write(f"R1 reads written: {total_reads_r1}\n")
        f.write(f"R1 total bases: {total_bases_r1}\n")
        f.write(f"R1 total introduced errors: {total_errors_r1}\n")
        f.write(f"R1 observed base error rate: {error_rate_r1:.6f}\n\n")
        f.write(f"R2 reads written: {total_reads_r2}\n")
        f.write(f"R2 total bases: {total_bases_r2}\n")
        f.write(f"R2 total introduced errors: {total_errors_r2}\n")
        f.write(f"R2 observed base error rate: {error_rate_r2:.6f}\n")


# =========================
# MAIN
# =========================

def main() -> None:
    rng = random.Random(SEED)

    fragments = read_fasta_records(INPUT_FASTA)
    selected_fragments = choose_fragments_to_sequence(
        fragments=fragments,
        sampling_fraction=FRAGMENT_SAMPLING_FRACTION,
        max_fragments=MAX_FRAGMENTS_TO_SEQUENCE,
        rng=rng,
    )

    r1_records = []
    r2_records = []

    total_errors_r1 = 0
    total_errors_r2 = 0
    total_bases_r1 = 0
    total_bases_r2 = 0

    read_index = 1

    for fragment in selected_fragments:
        if PAIRED_END:
            r1, r2 = make_paired_end_reads(
                fragment=fragment,
                read_length=READ_LENGTH,
                rng=rng,
                read_index=read_index,
            )

            h1, s1, q1, e1 = r1
            h2, s2, q2, e2 = r2

            r1_records.append((h1, s1, q1))
            r2_records.append((h2, s2, q2))

            total_errors_r1 += e1
            total_errors_r2 += e2
            total_bases_r1 += len(s1)
            total_bases_r2 += len(s2)

        else:
            h1, s1, q1, e1 = make_single_end_read(
                fragment=fragment,
                read_length=READ_LENGTH,
                rng=rng,
                read_index=read_index,
            )

            r1_records.append((h1, s1, q1))
            total_errors_r1 += e1
            total_bases_r1 += len(s1)

        read_index += 1

    write_fastq(r1_records, OUTPUT_FASTQ_R1)
    if PAIRED_END:
        write_fastq(r2_records, OUTPUT_FASTQ_R2)

    write_report(
        path=OUTPUT_REPORT,
        total_fragments_input=len(fragments),
        total_fragments_selected=len(selected_fragments),
        total_reads_r1=len(r1_records),
        total_reads_r2=len(r2_records),
        total_bases_r1=total_bases_r1,
        total_bases_r2=total_bases_r2,
        total_errors_r1=total_errors_r1,
        total_errors_r2=total_errors_r2,
    )

    print(f"Input fragments: {len(fragments)}")
    print(f"Fragments selected for sequencing: {len(selected_fragments)}")
    print(f"R1 reads written: {len(r1_records)}")
    if PAIRED_END:
        print(f"R2 reads written: {len(r2_records)}")
    print(f"R1 FASTQ: {OUTPUT_FASTQ_R1}")
    if PAIRED_END:
        print(f"R2 FASTQ: {OUTPUT_FASTQ_R2}")
    print(f"Report: {OUTPUT_REPORT}")


if __name__ == "__main__":
    main()