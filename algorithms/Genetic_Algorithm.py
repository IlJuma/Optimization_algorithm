# -*- coding: utf-8 -*-
"""
Genetic Algorithm for DNA Fragment Assembly
"""

import time
import numpy as np
import matplotlib.pyplot as plt

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem
from model.config import (
    SEED,
    MAX_EVALUATIONS,
    MAX_TIME_SEC,
    GA_POP_SIZE,
    GA_NUM_GENERATIONS,
    GA_MUTATION_RATE,
    GA_CROSSOVER_RATE,
    GA_ELITISM,
)

print("Import done")


# =========================
# GA PARAMETERS
# =========================

POP_SIZE = GA_POP_SIZE
CROSSOVER_RATE = GA_CROSSOVER_RATE
MUTATION_RATE = GA_MUTATION_RATE
NUM_GENERATIONS = GA_NUM_GENERATIONS
ELITISM = GA_ELITISM


# =========================
# HELPER FUNCTIONS
# =========================

def base_fragment_id(oriented_fragment_id):
    if oriented_fragment_id.endswith("_F") or oriented_fragment_id.endswith("_R"):
        return oriented_fragment_id[:-2]
    raise ValueError(f"Invalid oriented fragment id: {oriented_fragment_id}")


def random_orientation(rng):
    return rng.choice(["F", "R"])


def build_random_oriented_solution(problem, rng):
    """
    Build a valid random solution:
    - each base fragment is used exactly once
    - orientation is chosen randomly for each fragment
    """
    base_fragment_ids = problem.base_fragment_ids.copy()
    rng.shuffle(base_fragment_ids)

    solution = []
    for fragment_id in base_fragment_ids:
        solution.append(f"{fragment_id}_{random_orientation(rng)}")

    return solution


def calculate_fitness(solution, problem):
    return problem.evaluate(solution)


def crossover(parent1, parent2, rng):
    """
    Order-preserving crossover on base fragments, with orientation inherited
    from the parent contributing the chosen occurrence.
    """
    if rng.random() >= CROSSOVER_RATE:
        return parent1.copy(), parent2.copy()

    n = len(parent1)
    point = rng.integers(1, n)

    prefix1 = parent1[:point]
    prefix2 = parent2[:point]

    used1 = {base_fragment_id(x) for x in prefix1}
    used2 = {base_fragment_id(x) for x in prefix2}

    child1 = prefix1 + [x for x in parent2 if base_fragment_id(x) not in used1]
    child2 = prefix2 + [x for x in parent1 if base_fragment_id(x) not in used2]

    return child1, child2


def mutate(solution, rng):
    """
    Mutation uses two operators:
    - swap two fragments
    - flip orientation of a fragment
    """
    mutated = solution.copy()
    n = len(mutated)

    for i in range(n):
        if rng.random() < MUTATION_RATE:
            if rng.random() < 0.5:
                j = rng.integers(0, n - 1)
                if j >= i:
                    j += 1
                mutated[i], mutated[j] = mutated[j], mutated[i]
            else:
                if mutated[i].endswith("_F"):
                    mutated[i] = mutated[i][:-2] + "_R"
                else:
                    mutated[i] = mutated[i][:-2] + "_F"

    return mutated


# =========================
# MAIN GA FUNCTION
# =========================

print("Calculating GA function")


def optimize(problem=None, config=None, rng=None):
    if problem is None:
        fragments = load_fragments("data/fasta/fragments.fasta")
        problem = AssemblyProblem(fragments)

    if rng is None:
        rng = np.random.default_rng(SEED)

    pop_size = getattr(config, "GA_POP_SIZE", POP_SIZE) if config is not None else POP_SIZE
    crossover_rate = getattr(config, "GA_CROSSOVER_RATE", CROSSOVER_RATE) if config is not None else CROSSOVER_RATE
    mutation_rate = getattr(config, "GA_MUTATION_RATE", MUTATION_RATE) if config is not None else MUTATION_RATE
    num_generations = getattr(config, "GA_NUM_GENERATIONS", NUM_GENERATIONS) if config is not None else NUM_GENERATIONS
    elitism = getattr(config, "GA_ELITISM", ELITISM) if config is not None else ELITISM
    max_evaluations = getattr(config, "MAX_EVALUATIONS", MAX_EVALUATIONS) if config is not None else MAX_EVALUATIONS
    max_time_sec = getattr(config, "MAX_TIME_SEC", MAX_TIME_SEC) if config is not None else MAX_TIME_SEC

    start_time = time.time()
    evaluations = 0

    population = [build_random_oriented_solution(problem, rng) for _ in range(pop_size)]

    best_solution = None
    best_score = float("inf")

    history = []
    current_cost_history = []
    contigs_history = []
    overlap_history = []

    for gen in range(num_generations):
        if evaluations >= max_evaluations:
            break
        if time.time() - start_time > max_time_sec:
            break

        # Evaluate population
        fitness_scores = []
        for sol in population:
            score = calculate_fitness(sol, problem)
            fitness_scores.append(score)
            evaluations += 1

            if evaluations >= max_evaluations:
                break
            if time.time() - start_time > max_time_sec:
                break

        fitness_scores = np.array(fitness_scores)

        if len(fitness_scores) == 0:
            break

        # Trim population if time/eval limit interrupted evaluation
        population = population[:len(fitness_scores)]

        gen_best_score = float(fitness_scores.min())
        history.append(gen_best_score)
        current_cost_history.append(gen_best_score)

        best_idx = int(np.argmin(fitness_scores))
        best_sol = population[best_idx]

        contigs_history.append(problem.count_contigs(best_sol))
        overlap_history.append(problem.total_overlap(best_sol))

        if gen_best_score < best_score:
            best_score = gen_best_score
            best_solution = best_sol.copy()

        # Selection
        inv_fitness = 1 / (fitness_scores + 1e-6)
        probs = inv_fitness / inv_fitness.sum()
        selected_indices = rng.choice(
            range(len(population)),
            size=len(population),
            replace=True,
            p=probs,
        )

        # Crossover
        offspring = []
        for i in range(0, len(population) - 1, 2):
            p1 = population[selected_indices[i]]
            p2 = population[selected_indices[i + 1]]

            old_rate = CROSSOVER_RATE
            globals()["CROSSOVER_RATE"] = crossover_rate
            c1, c2 = crossover(p1, p2, rng)
            globals()["CROSSOVER_RATE"] = old_rate

            offspring.extend([c1, c2])

        # Mutation
        old_mutation_rate = MUTATION_RATE
        globals()["MUTATION_RATE"] = mutation_rate
        offspring = [mutate(child, rng) for child in offspring]
        globals()["MUTATION_RATE"] = old_mutation_rate

        # Survivor selection
        combined = population + offspring
        combined_scores = [calculate_fitness(sol, problem) for sol in combined]
        evaluations += len(combined)

        ranked = sorted(zip(combined, combined_scores), key=lambda x: x[1])
        population = [sol for sol, _ in ranked[:pop_size]]

        # Elitism
        if elitism > 0:
            elites = [sol for sol, _ in ranked[:elitism]]
            population[-elitism:] = elites

        # Logging
        if (gen + 1) % 20 == 0 or gen == 0:
            print(f"Generation {gen + 1}/{num_generations} | Best cost: {best_score:.2f}")

        if evaluations >= max_evaluations:
            break
        if time.time() - start_time > max_time_sec:
            break

    runtime = time.time() - start_time

    final_breaks = problem.count_breaks(best_solution) if best_solution else 0
    final_contigs = problem.count_contigs(best_solution) if best_solution else 0
    final_total_overlap = problem.total_overlap(best_solution) if best_solution else 0

    return {
        "method": "Genetic Algorithm",
        "best_solution": best_solution,
        "best_score": best_score,
        "history": history,
        "current_cost_history": current_cost_history,
        "contigs_history": contigs_history,
        "overlap_history": overlap_history,
        "evaluations": evaluations,
        "runtime_sec": runtime,
        "population_size": pop_size,
        "num_generations": num_generations,
        "best_breaks": final_breaks,
        "best_contigs": final_contigs,
        "best_total_overlap": final_total_overlap,
        "breaks": final_breaks,
        "contigs": final_contigs,
        "total_overlap": final_total_overlap,
    }


def genetic_algorithm(problem=None, config=None, rng=None):
    return optimize(problem=problem, config=config, rng=rng)


# =========================
# RUN SCRIPT
# =========================

if __name__ == "__main__":
    print("Loading fragments...")
    fragments = load_fragments("data/fasta/fragments.fasta")
    problem = AssemblyProblem(fragments)

    print("Running Genetic Algorithm...")
    result = optimize(problem=problem)

    print("\n=== BEST SOLUTION FOUND ===")
    print(f"Best score (total cost): {result['best_score']:.2f}")
    print(f"Number of fragments: {problem.n_fragments}")
    print(f"Fragment IDs order: {result['best_solution']}")

    # 1. Cost convergence
    plt.figure()
    plt.plot(result["history"])
    plt.xlabel("Generation")
    plt.ylabel("Best cost")
    plt.title("GA Convergence")
    plt.grid(True)

    # 2. Contigs vs iteration
    plt.figure()
    plt.plot(result["contigs_history"])
    plt.xlabel("Generation")
    plt.ylabel("Number of contigs")
    plt.title("Contigs vs Iteration")
    plt.grid(True)

    # 3. Overlap vs iteration
    plt.figure()
    plt.plot(result["overlap_history"])
    plt.xlabel("Generation")
    plt.ylabel("Total overlap")
    plt.title("Overlap vs Iteration")
    plt.grid(True)

    plt.show()