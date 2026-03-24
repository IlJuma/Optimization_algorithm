Here is your **clean, copy-pasteable Markdown version** with properly formatted LaTeX:

---

# 🧬 Project Summary: Fragment Ordering from Paired-End Sequencing Reads

## 1. Problem Statement

Given a set of sequencing reads

$$
V = \{1,\dots,n\}
$$

originating from **paired-end sequencing of multiple copies of a genome**, find an ordering (permutation)

$$
\pi \in S_n
$$

that reconstructs the original genome as accurately as possible.

### Key properties of the data

- reads originate from **fragmented copies of the same genome**
- regions are **covered multiple times (coverage redundancy)**
- reads come in **paired-end relationships**
- paired reads originate from the **same DNA fragment** but:
  - may be separated by a gap
  - may be bridged by other reads
- fragment lengths are size-selected to:

$$
L_{\text{frag}} \in [250, 450]\ \text{bp}
$$

- read lengths are:

$$
L_{\text{read}} \le 150\ \text{bp}
$$

---

## Problem type

This is a **combinatorial optimization problem** analogous to the Traveling Salesman Problem (TSP):

* nodes = reads
* edge weights = pairwise alignment cost
* solution = permutation
* objective = minimize total adjacency cost

---

## 2. Feasible Set (Constraints)

### Core constraint

$$
\pi \in S_n
$$

Each read appears **exactly once**.

---

## Biological constraints

For each adjacent pair $(i,j)$:

---

### (1) Minimum overlap

$$
\text{overlap}(i,j) \ge o_{\min}
$$

Ensures only plausible sequence adjacencies.

---

### (2) Phred-based mismatch constraint

Let:

- $M(i,j)$ = observed mismatches
- $p_k = 10^{-Q_k/10}$ = base error probability
- $E(i,j) = \sum_{k \in \text{overlap}} p_k$ = expected mismatches over the overlapping region

Constraint:

$$
M(i,j) \le \alpha \, E(i,j), \qquad \alpha \approx 2\text{--}5
$$

Ensures mismatches are consistent with sequencing error.

---

### (3) Paired-end distance constraint (IMPORTANT)

Each read $i$ has a paired read $\text{pair}(i)$.

Paired reads:

- originate from the **same DNA fragment**
- are separated by an **unknown gap**
- may have **other reads between them**
- must satisfy a **distance constraint derived from fragment length and read length**

#### Gap definition

Let the gap between paired reads be:

$$
\text{gap}(i,\text{pair}(i)) = L_{\text{frag}} - L_{\text{R1}} - L_{\text{R2}}
$$

Given:

$$
L_{\text{frag}} \in [250,450], \qquad L_{\text{R1}} = L_{\text{R2}} = 150
$$

we obtain:

$$
-50 \le \text{gap} \le 150
$$

This constraint is directly derived from the fragmentation and sequencing parameters used in the simulation.

---

#### Constraint formulation

Let $\hat d(i,j)$ denote the estimated genomic distance between reads.

We enforce:

$$
-50 \le \hat d(i,\text{pair}(i)) \le 150
$$

#### Simplified version (recommended)

In practice, the upper bound is most important:

$$
\hat d(i,\text{pair}(i)) \le 150
$$

---

#### Distance estimation

Approximate:

$$
\hat d(i,j) \approx \text{sum of overlaps and gaps along the path between } i \text{ and } j
$$

For direct adjacency:

$$
\hat d(i,j) \approx |i| + |j| - \text{overlap}(i,j)
$$

---

#### Interpretation

- gap < 0 → reads overlap  
- gap = 0 → reads meet  
- gap > 0 → unsequenced region  
- large distances → biologically impossible → rejected

---

### Final feasible set

$$
\mathcal{X} =
\{ \pi \in S_n :
\text{overlap}(i,j) \ge o_{\min},\;
M(i,j) \le \alpha E(i,j),\;
-50 \le \hat d(i,\text{pair}(i)) \le 150
\}
$$

---

## 3. Cost Function (Permutation Cost)

### Pairwise cost

$$
c_{ij} = L - \text{overlap}(i,j)
$$

where

$$
L = \min(|i|, |j|)
$$

### Total permutation cost

$$
f(\pi) = \sum_{k=1}^{n-1} c_{\pi_k,\pi_{k+1}}
$$

---

### Interpretation

* large overlap → low cost
* small overlap → high cost

---

### Optional extension

$$
c_{ij} = L - \text{overlap}(i,j) + \lambda \cdot \text{mismatches}(i,j)
$$

⚠️ Only valid if mismatch constraints are enforced.

---

## 4. Optimization Goal

The objective is to minimize the total adjacency cost over all feasible permutations.

$$
\min_{\pi \in \mathcal{X}} f(\pi)
= \sum_{k=1}^{n-1} c_{\pi_k,\pi_{k+1}}
$$

where:

- $\pi$: a permutation (ordering) of all reads
- $\mathcal{X}$: the feasible set of permutations satisfying the biological constraints
- $f(\pi)$: the total cost of a permutation
- $n$: the number of reads
- $k$: the index of a position in the permutation, from $1$ to $n-1$
- $\pi_k$: the read at position $k$
- $\pi_{k+1}$: the read immediately after $\pi_k$
- $c_{ij}$: the pairwise cost of placing read $j$ after read $i$

---

### Interpretation

The objective sums the costs of all **adjacent read pairs** in the ordering.  
A good reconstruction is one where consecutive reads have **large overlaps and low mismatch**, resulting in a low total cost.

---

## 5. Computational Budget Constraint

$$
C \le B
$$

where:

- $C$ = computational cost (e.g. evaluations)
- $B$ = fixed budget

---

### Interpretation

> Study genome reconstruction quality under limited computational resources.

---

## 6. Optimization Methods

### Random Search

* random permutations
* keep best

---

### Simulated Annealing

- local modifications (swap, reverse)
- probabilistic acceptance:

$$
P = e^{-\Delta/T}
$$

---

### Genetic Algorithm

* population of permutations
* selection, crossover, mutation

---

## 7. Evaluation Metrics

### (A) Objective-based

- best cost:

$$
\min f(\pi)
$$

- average cost
- convergence curves

---

### (B) Biological accuracy

Using known ground truth:

- adjacency accuracy:

$$
\frac{\text{correct adjacent pairs}}{n-1}
$$

- position correlation

---

### (C) Robustness

* variance across runs
* sensitivity to:

  * sequencing error
  * coverage level
  * parameter choices

---

### (D) Efficiency

* performance vs computational budget
* runtime / evaluations

---

## 8. Experimental Design

* subset size: 50–200 reads
* fixed computational budget $B$
* multiple runs (20–50)

Compare:

* best result
* mean ± variance
* convergence behavior

---

## 9. Expected Insights

* Random Search: baseline, inefficient
* Simulated Annealing: strong local optimization
* Genetic Algorithm: robust global exploration

---

## 10. Modeling Considerations

### Avoid

* unrealistic full overlaps
* ignoring sequencing errors
* violating paired-end structure

### Additional consideration: coverage redundancy

Due to fragmentation of multiple genome copies, most genomic regions are covered by multiple reads. This redundancy:

- introduces **multiple valid overlaps**
- increases **ambiguity in ordering**
- creates a **non-unique optimal solution space**

The optimization therefore seeks a **globally consistent ordering**, rather than a single exact reconstruction.

---

### Ensure

* biologically valid overlaps
* error-aware matching
* consistent paired-read placement
* handling of redundant coverage

---

## 11. Alternative Modeling Options

### (A) Alternative cost functions

1. Alignment score:

$$
c_{ij} = -(\text{matches} - \mu \cdot \text{mismatches})
$$

2. Normalized mismatch:

$$
c_{ij} = L - \text{overlap}(i,j) + \lambda \cdot \frac{\text{mismatches}(i,j)}{\text{overlap}(i,j)}
$$

---

### (B) Graph formulation

* weighted directed graph
* path optimization

---

### (C) Multi-objective optimization

$$
\min (f(\pi), C)
$$

---

### (D) Advanced quality modeling

* Phred-weighted mismatches
* likelihood-based scoring

---

### (E) Soft paired-end constraint

$$
c ;+=; \gamma \cdot \max(0, \hat d(i, \text{pair}(i)) - X)
$$

---

## 12. Final Positioning

> This project formulates genome fragment ordering as a constrained combinatorial optimization problem and evaluates stochastic optimization methods under limited computational resources. The model incorporates biologically informed constraints based on sequence overlap, sequencing error probabilities, paired-end geometry (including gap constraints derived from fragment and read lengths), and coverage redundancy, ensuring realistic solutions while enabling meaningful comparison of optimization strategies. The paired-end distance constraint is derived directly from the simulation parameters: with fragment sizes between 250–450 bp and read lengths of 150 bp, the maximum gap between paired reads is 150 bp, while shorter fragments may produce overlapping reads.