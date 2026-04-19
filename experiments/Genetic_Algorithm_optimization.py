# -*- coding: utf-8 -*-
"""
Genetic Algorithm for DNA Fragment Assembly + Hyperparameter Search
Author: Inés Calvo Esteva
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem
from model.config import (
    SEED,
    GA_POP_SIZES,
    GA_NUM_GENERATIONS,
    GA_MUTATION_RATES,
    GA_CROSSOVER_RATES,
    GA_ELITISM_VALUES,
    MAX_EVALUATIONS,
    MAX_TIME_SEC,
)
from algorithms.genetic_algorithm import optimize


# =========================
# MAIN: HYPERPARAMETER SEARCH
# =========================
if __name__ == "__main__":

    print("Loading fragments...")
    fragments = load_fragments("data/fasta/fragments.fasta")
    problem = AssemblyProblem(fragments)

    output_dir = "reports/genetic_algorithm"
    os.makedirs(output_dir, exist_ok=True)

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

                    config = SimpleNamespace(
                        SEED=SEED,
                        GA_POP_SIZE=pop_size,
                        GA_NUM_GENERATIONS=GA_NUM_GENERATIONS,
                        GA_MUTATION_RATE=mutation,
                        GA_CROSSOVER_RATE=crossover,
                        GA_ELITISM=elitism,
                        MAX_EVALUATIONS=MAX_EVALUATIONS,
                        MAX_TIME_SEC=MAX_TIME_SEC,
                    )

                    shared_rng = np.random.default_rng(SEED)

                    result = optimize(
                        problem=problem,
                        config=config,
                        shared_rng=shared_rng,
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
    with open(os.path.join(output_dir, "best_config.json"), "w") as f:
        json.dump({
            "best_config": best_config,
            "best_score": float(best_global)
        }, f, indent=4)

    with open(os.path.join(output_dir, "all_results.json"), "w") as f:
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
    plt.xlabel("Evaluation")
    plt.ylabel("Cost")
    plt.grid()
    plt.savefig(os.path.join(output_dir, "best_cost_convergence.png"), bbox_inches="tight")
    plt.close()

    plt.figure()
    plt.plot(result["contigs_history"], color="blue")
    plt.title("Contigs Evolution")
    plt.xlabel("Evaluation")
    plt.ylabel("Contigs")
    plt.grid()
    plt.savefig(os.path.join(output_dir, "contigs_evolution.png"), bbox_inches="tight")
    plt.close()

    plt.figure()
    plt.plot(result["overlap_history"], color="green")
    plt.title("Overlap Evolution")
    plt.xlabel("Evaluation")
    plt.ylabel("Overlap")
    plt.grid()
    plt.savefig(os.path.join(output_dir, "overlap_evolution.png"), bbox_inches="tight")
    plt.close()

    print(f"\nResults saved in: {output_dir}")