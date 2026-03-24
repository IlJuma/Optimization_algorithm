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

## Structure

```
scripts/   # simulation scripts
fasta/     # generated sequences (ignored)
reports/   # simulation reports (ignored)
```

## Usage

Run in order:

```
python scripts/1_generate_genome.py
python scripts/2_fragment_chromosome.py
python scripts/3_simulate_illumina.py
```

## Notes

* Output folders (`fasta/`, `reports/`) are created automatically
* All outputs are excluded from version control
* Designed for testing genome reconstruction and assembly accuracy
