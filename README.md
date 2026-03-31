# Optimization Algorithm – Genome Simulation Pipeline

This project simulates a simplified DNA sequencing workflow:

1. **Genome generation**

   * Random DNA sequence (A/C/G/T only)
   * GC content variation across windows

2. **Fragmentation**

   * Random cut-site model (Poisson process)
   * Size selection (Illumina-like insert sizes)
   * Random recovery of fragments

3. **Sequencing simulation**

   * Paired-end reads
   * Phred-based error model
   * Provenance tracking (true genomic coordinates)

4. **Oracle coverage evaluation**

   * Ground-truth coordinate reconstruction from fragments and reads
   * Per-base coverage calculation across the reference genome
   * Gap and contig detection (uncovered vs covered regions)
   * Summary statistics (coverage depth, breadth, extremes)
   * Export of coverage tracks, intervals (TSV), and diagnostic plots

## Structure

```
pipeline/   # simulation scripts

- `pipeline/` – generation and validation of simulated genomic data
- `model/` – data loading, scoring, constraints, and problem representation
- `algorithms/` – random search, simulated annealing, and genetic algorithm
- `experiments/` – optimization runs and result analysis
- `data/fasta/` - generated sequences (ignored)
- `data/fastq/` - generated reads (ignored)
- `reports/` - simulation reports (ignored)
```

## Usage

Run the scripts in order; the first script initializes all required directories.

```
python pipeline/1_generate_chromosome.py
python pipeline/2_fragment_chromosome.py
python pipeline/3_simulate_illumina.py
python pipeline/4_oracle_coverage_evaluation.py
```

## Notes

* Output folders (`data/fasta/`, `data/fastq/`, `reports/`) are created automatically
* All outputs are excluded from version control
* Designed for testing genome reconstruction and assembly accuracy
