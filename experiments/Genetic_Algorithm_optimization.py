# -*- coding: utf-8 -*-
"""
Genetic Algorithm for DNA Fragment Assembly + Hyperparameter Search
Author: Inés Calvo Esteva
"""

import numpy as np
import matplotlib.pyplot as plt
import json

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem
from model.config import (
    SEED,
    GA_POP_SIZES,
    GA_NUM_GENERATIONS,
    GA_MUTATION_RATES,
    GA_CROSSOVER_RATES,
    GA_ELITISM_VALUES
)

# =========================
# FITNESS
# =========================
def calculate_fitness(solution, problem):
    return problem.evaluate(solution)


# =========================
# GA CORE (PARAMETRIZADO)
# =========================
def genetic_algorithm(problem, pop_size, mutation_rate, crossover_rate, elitism):

    rng = np.random.default_rng(SEED)

    population = [
        rng.permutation(problem.fragment_ids).tolist()
        for _ in range(pop_size)
    ]

    best_solution = None
    best_score = float("inf")

    history = []
    contigs_history = []
    overlap_history = []

    for gen in range(GA_NUM_GENERATIONS):

        fitness_scores = np.array([
            calculate_fitness(sol, problem) for sol in population
        ])

        history.append(fitness_scores.min())

        best_idx = np.argmin(fitness_scores)
        best_sol = population[best_idx]

        contigs_history.append(problem.count_contigs(best_sol))
        overlap_history.append(problem.total_overlap(best_sol))

        # ------------------------
        # Selection (roulette)
        # ------------------------
        inv_fitness = 1 / (fitness_scores + 1e-6)
        probs = inv_fitness / inv_fitness.sum()

        selected_idx = rng.choice(
            range(pop_size),
            size=pop_size,
            replace=True,
            p=probs
        )

        # ------------------------
        # Crossover
        # ------------------------
        offspring = []
        for i in range(0, pop_size, 2):
            p1 = population[selected_idx[i]]
            p2 = population[selected_idx[i + 1]]

            if rng.random() < crossover_rate:
                point = rng.integers(1, len(p1))
                c1 = p1[:point] + [x for x in p2 if x not in p1[:point]]
                c2 = p2[:point] + [x for x in p1 if x not in p2[:point]]
            else:
                c1, c2 = p1.copy(), p2.copy()

            offspring.extend([c1, c2])

        # ------------------------
        # Mutation
        # ------------------------
        mutated = []
        for sol in offspring:
            sol = sol.copy()
            for i in range(len(sol)):
                if rng.random() < mutation_rate:
                    j = rng.integers(0, len(sol))
                    sol[i], sol[j] = sol[j], sol[i]
            mutated.append(sol)

        # ------------------------
        # Survivor selection
        # ------------------------
        combined = population + mutated
        combined.sort(key=lambda s: calculate_fitness(s, problem))

        population = combined[:pop_size]

        # ------------------------
        # Elitism
        # ------------------------
        if elitism > 0:
            elites = combined[:elitism]
            population[-elitism:] = elites

        # ------------------------
        # Best update
        # ------------------------
        current_best = calculate_fitness(population[0], problem)
        if current_best < best_score:
            best_score = current_best
            best_solution = population[0]

        if gen % 20 == 0 or gen == 0:
            print(f"Gen {gen} | Best: {best_score:.2f}")

    return {
        "best_solution": best_solution,
        "best_score": best_score,
        "history": history,
        "contigs_history": contigs_history,
        "overlap_history": overlap_history,
    }


# =========================
# MAIN: HYPERPARAMETER SEARCH
# =========================
if __name__ == "__main__":

    print("Loading fragments...")
    fragments = load_fragments("data/fasta/fragments.fasta")
    problem = AssemblyProblem(fragments)

    best_global = float("inf")
    best_config = None
    best_result = None

    all_results = []

    print("\nRunning hyperparameter search...\n")

    for pop_size in GA_POP_SIZES:
        for mutation in GA_MUTATION_RATES:
            for crossover in GA_CROSSOVER_RATES:
                for elitism in GA_ELITISM_VALUES:

                    print(f"Testing -> pop={pop_size}, mut={mutation}, cross={crossover}, elit={elitism}")

                    result = genetic_algorithm(
                        problem,
                        pop_size,
                        mutation,
                        crossover,
                        elitism
                    )

                    score = result["best_score"]

                    print(f" -> Score: {score:.2f}\n")

                    all_results.append({
                        "pop_size": pop_size,
                        "mutation": mutation,
                        "crossover": crossover,
                        "elitism": elitism,
                        "score": float(score)
                    })

                    if score < best_global:
                        best_global = score
                        best_config = {
                            "pop_size": pop_size,
                            "mutation": mutation,
                            "crossover": crossover,
                            "elitism": elitism
                        }
                        best_result = result

    # =========================
    # SAVE RESULTS
    # =========================
    with open("best_config.json", "w") as f:
        json.dump({
            "best_config": best_config,
            "best_score": float(best_global)
        }, f, indent=4)

    with open("all_results.json", "w") as f:
        json.dump(all_results, f, indent=4)

    print("\n=== BEST CONFIG FOUND ===")
    print(best_config)
    print("Best score:", best_global)

    # =========================
    # PLOTS (BEST RUN)
    # =========================
    result = best_result

    plt.figure()
    plt.plot(result["history"], color="red")
    plt.title("Best Cost Convergence")
    plt.xlabel("Generation")
    plt.ylabel("Cost")
    plt.grid()

    plt.figure()
    plt.plot(result["contigs_history"], color="blue")
    plt.title("Contigs Evolution")
    plt.grid()

    plt.figure()
    plt.plot(result["overlap_history"], color="green")
    plt.title("Overlap Evolution")
    plt.grid()

    plt.show()