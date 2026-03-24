Here’s a **combined, clean plan** for both the **team code architecture** and the **development scaling strategy**.

---

# Team architecture + development size strategy

You have a good 4-person split, but the most important rule is:

> **All algorithms must use the exact same problem representation, scoring logic, constraints, and result format.**

Otherwise, the comparison will not be fair and integration will become messy.

The cleanest overall structure is:

1. **shared data / model layer**
2. **objective + constraints layer**
3. **algorithm layer**
4. **experiment / evaluation layer**

With 4 people, a natural division is:

* **Person A:** shared model, scoring, constraints, preprocessing
* **Person B:** random search
* **Person C:** simulated annealing
* **Person D:** genetic algorithm + experiment runner/integration

---

## 1. Shared code architecture

### Input / preprocessing

Use your existing scripts for chromosome generation, fragmentation, and sequencing. Add a shared loader layer that converts outputs into a common in-memory representation.

Suggested module:

```python
data_loader.py
```

This should produce:

* list of reads
* mapping from read ID to paired read ID
* read lengths
* quality scores
* optional ground-truth positions from the simulator
* optional copy indices

---

### Problem model

This is the central shared component that everyone uses.

Suggested module:

```python
problem.py
```

This should define:

* what a solution is
* how a solution is scored
* what constraints are enforced
* helper methods shared by all algorithms

A solution should be represented consistently, for example:

```python
solution = [read_id_1, read_id_2, ..., read_id_n]
```

or more efficiently:

```python
solution = [0, 5, 2, 1, ...]   # indices into the read list
```

Recommended practical choice:

* use **indices internally**
* convert to read IDs only for logs / reporting

---

### Algorithms

Each method lives in its own file:

```python
random_search.py
simulated_annealing.py
genetic_algorithm.py
```

Each should expose the same interface:

```python
def optimize(problem, config, rng):
    ...
    return result
```

Each `result` should contain:

* best solution
* best score
* score history
* evaluation count
* runtime
* optional diagnostics

---

### Experiment runner

Use one shared script to run all methods under identical conditions.

Suggested module:

```python
run_experiments.py
```

This should:

* load the simulated data
* create the `AssemblyProblem`
* call each algorithm with the same budget
* save results
* generate tables / plots

This is what guarantees a fair comparison.

---

## 2. Recommended team split

### Person A — shared model, constraints, scoring, preprocessing

This is the most central role.

Responsibilities:

* define read representation
* define solution representation
* implement pairwise overlap and mismatch logic
* implement global scoring
* implement constraint checking
* precompute reusable structures

This person should provide functions like:

* `evaluate(solution)`
* `pair_cost(i, j)`
* `is_edge_feasible(i, j)`
* `pair_distance_penalty(solution)`
* `random_solution(rng)`

This is the foundation that all other team members depend on.

---

### Person B — random search

Responsibilities:

* implement random search baseline
* use the shared evaluation interface
* return results in the common result format

This person can also help with:

* baseline sanity tests
* debugging result formatting
* checking whether the score improves at all versus random

---

### Person C — simulated annealing

Responsibilities:

* implement neighbor generation
* implement annealing schedule
* integrate with shared evaluator

Important choices:

* swap move
* insertion move
* segment reversal
* cooling schedule
* acceptance probability

This person may also provide a shared helper:

* `propose_neighbor(solution, rng)`

---

### Person D — genetic algorithm + integration

Responsibilities:

* implement permutation-based GA
* selection, crossover, mutation
* population evolution
* integrate GA into experiment runner
* potentially own result aggregation and plotting

Because GA is the most complex algorithm, this is a good standalone role.

---

## 3. Global interfaces everyone must agree on first

Before anyone codes, agree on these shared structures.

### Read representation

Example:

```python
Read = {
    "id": str,
    "seq": str,
    "qual": list[int],
    "pair_id": int | None,
    "true_pos": int | None,
    "copy_index": int | None,
}
```

Minimum fields:

* id
* sequence
* qualities
* pair ID

Optional but useful:

* true position
* originating genome copy

---

### Solution representation

Pick one and stick to it.

Recommended:

```python
solution = list[int]
```

where each integer is the index of a read in `problem.reads`.

---

### Problem object

All algorithms should use the same object.

Example:

```python
class AssemblyProblem:
    reads: list[Read]
    overlap_matrix: list[list[int]]
    mismatch_matrix: list[list[int]]
    feasible_edge_matrix: list[list[bool]]
    pair_map: dict[int, int]

    def evaluate(self, solution: list[int]) -> float: ...
    def pair_cost(self, i: int, j: int) -> float: ...
    def is_edge_feasible(self, i: int, j: int) -> bool: ...
    def pair_distance_penalty(self, solution: list[int]) -> float: ...
    def random_solution(self, rng) -> list[int]: ...
```

This is the single most important architectural decision.

---

### Common result format

All algorithms should return the same structure.

Example:

```python
Result = {
    "method": str,
    "best_solution": list[int],
    "best_score": float,
    "history": list[float],
    "evaluations": int,
    "runtime_sec": float,
}
```

That makes experiment code simple and fair.

---

## 4. Global parameters to define together

Put them in one shared file:

```python
config.py
```

### Data / model parameters

```python
MIN_OVERLAP = 20
ALPHA_MISMATCH = 3.0
MIN_PAIR_GAP = -50
MAX_PAIR_GAP = 150
USE_MISMATCH_PENALTY = False
LAMBDA_MISMATCH = 1.0
BIG_M = 1e9
```

### Instance / simulation size

Since you decided to control difficulty via **genome length**, this section should focus on that instead of read subsets.

```python
GENOME_LENGTH = 10_000
N_GENOME_COPIES = 20
RANDOM_SEED = 123
```

### Evaluation budget

```python
MAX_EVALUATIONS = 10_000
MAX_TIME_SEC = 60
```

Every algorithm should use the same budget.

### Random search

```python
RS_NUM_SAMPLES = 10_000
```

### Simulated annealing

```python
SA_INITIAL_TEMP = 10.0
SA_COOLING_RATE = 0.995
SA_MIN_TEMP = 1e-3
SA_NEIGHBORHOOD = "swap"
```

### Genetic algorithm

```python
GA_POP_SIZE = 100
GA_NUM_GENERATIONS = 200
GA_MUTATION_RATE = 0.1
GA_CROSSOVER_RATE = 0.8
GA_ELITISM = 2
```

---

## 5. What should be precomputed centrally

To avoid duplicated work, Person A should precompute:

### Overlap matrix

For every pair `(i, j)`:

* overlap length

### Mismatch matrix

For every pair `(i, j)`:

* mismatches in overlap

### Feasible edge matrix

For every pair `(i, j)`:

* whether overlap + mismatch constraints are satisfied

### Pair mapping

For each read:

* paired read index

This allows the algorithms to focus on search, not bioinformatics plumbing.

---

## 6. Recommended evaluation design

### Important design choice: hard constraints vs soft penalties

Recommended:

* use **hard constraints** for local edge feasibility:

  * minimum overlap
  * mismatch threshold
* use **soft penalties** for paired-end distance

This is easier to implement and gives a smoother search landscape.

So evaluation can be decomposed as:

```python
total_cost = adjacency_cost + pair_distance_penalty
```

and infeasible local edges can receive a huge penalty:

```python
BIG_M = 1e9
```

### Useful score decomposition

For debugging and analysis:

```python
score = (
    overlap_cost
    + mismatch_penalty
    + pair_distance_penalty
)
```

Optional detailed return:

```python
{
    "total": ...,
    "overlap_cost": ...,
    "pair_penalty": ...,
    "num_invalid_edges": ...,
}
```

This helps all team members interpret failures.

---

## 7. Development scaling strategy

Instead of subsetting reads, you should control problem size through **genome length**. That is a better choice because it preserves:

* coverage structure
* paired-end relationships
* redundancy patterns
* overall biological realism

### Why 1,000,000 bp is too large for development

With your current parameters:

* genome length = 1,000,000
* copies = 100
* mean cut spacing ≈ 350 bp

Per copy:

$$
\frac{1{,}000{,}000}{350} \approx 2857 \text{ fragments}
$$

Across 100 copies:

$$
2857 \times 100 \approx 285{,}700 \text{ fragments}
$$

After size selection and 25% recovery:

$$
\sim 70{,}000 \text{ fragments}
$$

That means:

* one permutation already has about 70,000 elements
* one evaluation requires about 70,000 adjacency checks
* simulated annealing and GA become very slow
* even random search becomes impractical

So **1,000,000 bp is too large for development**.

---

### Recommended progressive scaling

#### Stage 1 — debugging

```python
GENOME_LENGTH = 5_000
N_GENOME_COPIES = 20
```

Roughly a few hundred fragments total.
Use this for:

* testing interfaces
* verifying scores
* debugging constraints
* checking algorithm output format

#### Stage 2 — algorithm tuning

```python
GENOME_LENGTH = 20_000
N_GENOME_COPIES = 20
```

Roughly 1,500–2,000 fragments.
Use this for:

* tuning SA parameters
* tuning GA population / mutation
* measuring convergence behavior

#### Stage 3 — final experiments

```python
GENOME_LENGTH = 50_000
N_GENOME_COPIES = 50
```

Roughly several thousand fragments.
This is large enough to look meaningful, but still manageable if you precompute and sparsify well.

---

### Strong recommendation

Do **not** use read subsets if you can avoid it.
Use smaller simulated genomes instead.

This preserves the structure of the problem and makes the comparison more scientifically defensible.

A good report sentence is:

> For computational feasibility, experiments were conducted on simulated genomes of reduced length. This preserves key structural properties such as coverage redundancy and paired-end relationships while allowing meaningful comparison of optimization methods within a fixed computational budget.

---

## 8. Extra speed-up trick

Even at moderate size, do not allow a full dense (n^2) graph if you can help it.

Pre-filter candidate edges:

```python
overlap(i, j) >= MIN_OVERLAP
```

This gives you a **sparse graph** instead of evaluating every possible pair. That can make a very large difference for all three algorithms.

---

## 9. Practical folder structure

```text
project/
│
├── data/
├── reports/
├── scripts/
│   ├── 1_generate_chromosome.py
│   ├── 2_fragment_chromosome.py
│   ├── 3_simulate_illumina.py
│
├── model/
│   ├── config.py
│   ├── data_loader.py
│   ├── problem.py
│   ├── scoring.py
│   ├── constraints.py
│
├── algorithms/
│   ├── random_search.py
│   ├── simulated_annealing.py
│   ├── genetic_algorithm.py
│
├── experiments/
│   ├── run_experiments.py
│   ├── analyze_results.py
│
└── tests/
```

---

## 10. Best workflow

### Step 1

Agree on:

* solution representation
* cost function
* constraints
* result format
* development genome size

### Step 2

Person A finishes:

* `problem.py`
* `scoring.py`
* `constraints.py`
* `data_loader.py`

### Step 3

Algorithm people implement against a tiny mock instance first.

### Step 4

Integrate all methods into `run_experiments.py`.

### Step 5

Run all methods on the same genome lengths and the same evaluation budgets.

---

## Final recommendation

Do **not** let each person implement their own scoring logic.
That causes inconsistent comparisons and report problems.

Centralize:

* scoring
* constraints
* data representation
* result format

Then split the algorithms across team members.

### Best final split

* **Person A:** shared model, scoring, constraints, overlap preprocessing
* **Person B:** random search + baselines + sanity checks
* **Person C:** simulated annealing + neighborhood operators
* **Person D:** genetic algorithm + experiment runner + aggregation