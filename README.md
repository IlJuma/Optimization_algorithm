# Optimization of DNA Fragment Assembly

## 1. Project Description

This project formulates **DNA fragment assembly as a combinatorial optimization problem** and evaluates bio-inspired algorithms to reconstruct a genome from unordered fragments.

DNA is represented as a sequence over the alphabet {A, C, G, T}. We simulate a genome and generate multiple copies, which are randomly fragmented. Due to the fragmentation and recovery process, the resulting dataset has the following properties:

- fragments partially overlap  
- genome coverage is incomplete (gaps may exist)  
- redundancy occurs due to multiple copies  
- fragment orientation is unknown (forward or reverse complement)

The goal is to reconstruct the original chromosome from an **unordered set of fragments with unknown orientation**.

---

### Mathematical formulation

Let  
V = {F₁, F₂, …, Fₙ}  
be the set of fragments.

A solution is a permutation with orientation:  
π = (π₁, π₂, …, πₙ),   πₖ ∈ {Fᶠᵢ, Fʳᵢ}

Each fragment is used exactly once, in either forward or reverse orientation.

Fragments are connected using suffix-prefix overlap:  
overlap(i, j)

A minimum overlap threshold is enforced:  
overlap(i, j) ≥ o_min   (default: 20)

If not satisfied → **break (new contig)**.

---

### Cost function

For two fragments i, j, let:  
ℓᵢⱼ = min(|i|, |j|)  
P = γ · ℓᵢⱼ  

The pairwise cost is:

c(i, j) =  
    ℓᵢⱼ − overlap(i, j)          if overlap ≥ o_min  
    P − β · overlap(i, j)        if 0 < overlap < o_min  
    P                            if overlap = 0  

with:  
- γ = 2 (break penalty)  
- β = 1 (weak overlap scaling)  

The total cost:  
f(π) = Σ c(πₖ, πₖ₊₁)

Goal:  
minimize f(π)

Search space:  
|S| = n! · 2ⁿ

---

### Optimization methods

We compare:

- Random Search (baseline)  
- Simulated Annealing (local search with probabilistic acceptance)  
- Genetic Algorithm (population-based global search)  

All methods are evaluated under a fixed budget of **15,000 evaluations**.

---

## 2. Project Structure

```text
Optimisation_and_Bio-inspired_algorithms/
├── algorithms/
├── experiments/
├── model/
├── pipeline/
├── data/
│   ├── fasta/    (generated, ignored)
│   ├── fastq/    (optional, ignored)
├── reports/      (generated, ignored)
├── README.md
└── requirements.txt
```

---

## 3. Modules Overview

### pipeline/
Simulation of genome and fragment generation.

- generate_chromosome.py — generate synthetic genome  
- fragment_chromosome.py — fragment genome and randomize orientation  
- oracle_coverage_evaluation.py — evaluate ground-truth coverage  
- simulate_illumina.py (optional, not used in final model)  

---

### model/
Core problem definition.

- config.py — parameters (overlap threshold, penalties, experiment setup)  
- data_loader_frag.py — load fragments and generate orientations  
- problem.py — cost function and evaluation logic  

---

### algorithms/
Optimization methods.

- random_search.py — baseline random sampling  
- simulated_annealing.py — temperature-based local search  
- Genetic_Algorithm.py — evolutionary optimization  
- oracle_solution.py — ground-truth ordering (benchmark)  

---

### experiments/
Running and analyzing experiments.

- run_experiments.py — main experiment pipeline  
- analyze_results.py — comparison plots  
- analyze_random_search.py — RS diagnostics  
- analyze_simulated_annealing.py — SA diagnostics  
- Genetic_Algorithm_optimization.py — GA hyperparameter tuning  
- run_grid_search.py — parameter sweeps  
- evaluate_oracle.py — evaluate oracle solution  

---

## 4. Running the Project

### Setup environment

```bash
python3 -m venv .venv
source .venv/bin/activate   # Linux/macOS/WSL

pip install -r requirements.txt
```

### Step 1 — Generate fragments

```bash
python3 -m pipeline.generate_chromosome
python3 -m pipeline.fragment_chromosome
```

### Step 2 — Run optimization experiments

```bash
python3 -m experiments.run_experiments
```

### Step 3 — Analyze results

```bash
python3 -m experiments.analyze_results
```

### Optional analysis

```bash
python3 -m experiments.analyze_random_search
python3 -m experiments.analyze_simulated_annealing
python3 -m experiments.Genetic_Algorithm_optimization
``` 

---

## Notes

- Output data is stored in data/ and reports/ and is excluded from version control
- Commands for running the project are written for Linux/WSL environments using `python3`.
