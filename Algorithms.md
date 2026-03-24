## 6. Optimization Methods

All methods operate on the same representation:

* a solution is a permutation
  $\pi = (\pi_1, \dots, \pi_n)$
* each element $\pi_k$ is a **single read** (R1 or R2)
* paired reads are linked via `pair(i)` but may be separated in the permutation

The objective function is:

$$
f(\pi) = \sum_{k=1}^{n-1} c_{\pi_k,\pi_{k+1}} + \text{pair-penalty}(\pi)
$$

where the pair-penalty enforces the constraint:

$$
-50 \le \hat d(i, \text{pair}(i)) \le 150
$$

---

## 6.1 Random Search

### Idea

Random search explores the solution space by sampling random permutations of reads and keeping the best one found.

It ignores problem structure and serves as a **baseline method**.

---

### Algorithm

1. Initialize:

   * sample a random permutation $\pi$
   * compute cost $f(\pi)$
   * store as best solution

2. Repeat until computational budget is exhausted:

   * sample a new random permutation $\pi'$
   * compute cost $f(\pi')$
   * if $f(\pi') < f(\pi)$:

     * update best solution

3. Return best solution found

---

### Implementation notes

* permutations must include **all reads exactly once**
* paired reads are **not enforced structurally**, but evaluated via penalty
* cost evaluation includes:

  * adjacency overlap cost
  * mismatch constraint
  * paired-read distance penalty

---

### Use of cost function

* evaluates full permutations
* uses cost only for comparison
* no local improvement

---

### Characteristics

* simple baseline
* no parameters
* inefficient (random exploration)
* high variance

---

## 6.2 Simulated Annealing

### Idea

Simulated annealing is a **local search method** that iteratively improves a permutation, while occasionally accepting worse solutions to escape local minima.

---

### Algorithm

1. Initialize:

   * sample a random permutation $\pi$
   * compute cost $f(\pi)$
   * set initial temperature $T = T_0$

2. Repeat until stopping condition:

   * generate a neighboring permutation $\pi'$
     (e.g. swap, insert, reverse)

   * compute cost difference:
     $$
     \Delta = f(\pi') - f(\pi)
     $$

   * if $\Delta < 0$:

     * accept $\pi'$

   * else:

     * accept $\pi'$ with probability:
       $$
       P = \exp(-\Delta / T)
       $$

   * decrease temperature:
     $T \leftarrow \alpha T$

3. Return best solution found

---

### Neighborhood moves

Typical moves on permutations:

* swap two reads
* insert read at a new position
* reverse subsequence

These create **small local changes**.

---

### Handling paired reads

* paired reads may move independently
* moves may temporarily violate pair distance
* penalty in $f(\pi)$ discourages invalid configurations

---

### Use of cost function

* evaluates local modifications
* uses cost difference $\Delta$ to guide acceptance
* balances exploration vs exploitation

---

### Characteristics

* efficient use of computational budget
* escapes local minima via probabilistic acceptance
* sensitive to:

  * initial temperature $T_0$
  * cooling rate $\alpha$

---

### Example

```
Start: [R3, R1, R5, R2, R4]   cost = 75

Swap → [R1, R3, R5, R2, R4]   cost = 68   (accept)
Swap → [R1, R2, R5, R3, R4]   cost = 60   (accept)
Swap → [R1, R2, R3, R5, R4]   cost = 50   (accept)
Swap → [R1, R2, R3, R4, R5]   cost = 40   (optimal)
```

---

## 6.3 Genetic Algorithm

### Idea

Genetic algorithms evolve a **population of permutations** using selection, crossover, and mutation.

They explore multiple regions of the search space simultaneously.

---

### Algorithm

1. Initialize population:

   * generate $N$ random permutations

2. Repeat for each generation:

   * evaluate each individual using $f(\pi)$

   * select parents (prefer lower cost)

   * apply crossover to produce offspring

   * apply mutation

   * form new population (optionally keep best individuals)

3. Return best solution found

---

### Representation

* each individual = permutation of reads
* fitness defined as:

  * lower cost = higher fitness
    (e.g. fitness = $-f(\pi)$)

---

### Selection

Common strategies:

* tournament selection
* rank-based selection

---

### Crossover (permutation-safe)

* order crossover (OX)
* partially matched crossover (PMX)

These preserve valid permutations.

---

### Mutation

* swap two reads
* reverse subsequence
* insert move

Maintains diversity.

---

### Handling paired reads

* pairs are **not fixed together**
* GA may separate them
* constraint enforced via penalty:

  * large distances → high cost → selection removes them

---

### Use of cost function

* evaluates all individuals
* drives selection toward better permutations
* used as fitness

---

### Characteristics

* global exploration
* robust across runs
* slower per iteration than SA
* requires parameter tuning:

  * population size
  * mutation rate
  * crossover rate

---

### Example

```
Initial population:
P1: [R1, R3, R2, R5, R4]   cost 65
P2: [R3, R1, R5, R2, R4]   cost 75
P3: [R1, R2, R4, R3, R5]   cost 55

Selection → P1, P3

Crossover:
[R1, R2, R3, R5, R4]

Mutation:
[R1, R2, R3, R4, R5]

→ cost = 40
```

---

## 6.4 Comparison of Methods

| Method              | Strengths                       | Weaknesses                |
| ------------------- | ------------------------------- | ------------------------- |
| Random Search       | Simple baseline                 | Very inefficient          |
| Simulated Annealing | Efficient, escapes local minima | Parameter-sensitive       |
| Genetic Algorithm   | Global exploration, robust      | Higher computational cost |

---

## 6.5 Role of the Cost Function

All methods rely on:

$$
f(\pi) = \sum_{k=1}^{n-1} c_{\pi_k,\pi_{k+1}} + \text{pair-penalty}(\pi)
$$

* Random Search: evaluates and compares permutations
* Simulated Annealing: uses **cost differences**
* Genetic Algorithm: uses cost as **fitness**