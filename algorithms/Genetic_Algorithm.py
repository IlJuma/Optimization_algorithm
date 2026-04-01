# -*- coding: utf-8 -*-
"""
Genetic Algorithm for DNA Fragment Assembly
"""
# Load libraries
import numpy as np
import matplotlib.pyplot as plt

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem
from model.config import (
    SEED,
    GA_POP_SIZE,
    GA_NUM_GENERATIONS,
    GA_MUTATION_RATE,
    GA_CROSSOVER_RATE,
    GA_ELITISM
)

print("Import done")


# GA PARAMETERS

POP_SIZE = GA_POP_SIZE
CROSSOVER_RATE = GA_CROSSOVER_RATE
MUTATION_RATE = GA_MUTATION_RATE
NUM_GENERATIONS = GA_NUM_GENERATIONS
ELITISM = GA_ELITISM


# HELPER FUNCTIONS

def calculate_fitness(solution, problem):
    return problem.evaluate(solution)


def crossover(parent1, parent2, rng):
    if rng.random() < CROSSOVER_RATE:
        n = len(parent1)
        point = rng.integers(1, n)
        child1 = parent1[:point] + [f for f in parent2 if f not in parent1[:point]]
        child2 = parent2[:point] + [f for f in parent1 if f not in parent2[:point]]
        return child1, child2
    else:
        return parent1.copy(), parent2.copy()


def mutate(solution, rng):
    n = len(solution)
    mutated = solution.copy()
    for i in range(n):
        if rng.random() < MUTATION_RATE:
            j = rng.integers(0, n-1)
            if j >= i:
                j += 1
            mutated[i], mutated[j] = mutated[j], mutated[i]
    return mutated



# MAIN GA FUNCTION

print("Calculating GA function")

def genetic_algorithm(problem):
    rng = np.random.default_rng(SEED)

    population = [rng.permutation(problem.fragment_ids).tolist() for _ in range(POP_SIZE)]

    best_solution = None
    best_score = float('inf')

    history = []
    contigs_history = []
    overlap_history = []

    for gen in range(NUM_GENERATIONS):
        # Evaluate population
        fitness_scores = np.array([calculate_fitness(sol, problem) for sol in population])
        history.append(fitness_scores.min())

        # Best solution of this generation
        best_idx = np.argmin(fitness_scores)
        best_sol = population[best_idx]

        # Metrics
        contigs_history.append(problem.count_contigs(best_sol))
        overlap_history.append(problem.total_overlap(best_sol))

        # Selection
        inv_fitness = 1 / (fitness_scores + 1e-6)
        probs = inv_fitness / inv_fitness.sum()
        selected_indices = rng.choice(
            range(POP_SIZE),
            size=POP_SIZE,
            replace=True,
            p=probs
        )

        # Crossover
        offspring = []
        for i in range(0, POP_SIZE, 2):
            p1 = population[selected_indices[i]]
            p2 = population[selected_indices[i+1]]
            c1, c2 = crossover(p1, p2, rng)
            offspring.extend([c1, c2])

        # Mutation
        offspring = [mutate(child, rng) for child in offspring]

        # Survivor selection
        combined = population + offspring
        combined.sort(key=lambda sol: calculate_fitness(sol, problem))
        population = combined[:POP_SIZE]

        # Elitism
        if ELITISM > 0:
            elites = combined[:ELITISM]
            population[-ELITISM:] = elites

        # Update best global solution
        current_best_score = calculate_fitness(population[0], problem)
        if current_best_score < best_score:
            best_score = current_best_score
            best_solution = population[0]

        # Logging
        if (gen + 1) % 20 == 0 or gen == 0:
            print(f"Generation {gen+1}/{NUM_GENERATIONS} | Best cost: {best_score:.2f}")

    return {
        "best_solution": best_solution,
        "best_score": best_score,
        "history": history,
        "contigs_history": contigs_history,
        "overlap_history": overlap_history,
        "population_size": POP_SIZE,
        "num_generations": NUM_GENERATIONS
    }



# RUN SCRIPT

if __name__ == "__main__":
    print("Loading fragments...")
    fragments = load_fragments("data/fasta/fragments.fasta")
    problem = AssemblyProblem(fragments)

    print("Running Genetic Algorithm...")
    result = genetic_algorithm(problem)

    print("\n=== BEST SOLUTION FOUND ===")
    print(f"Best score (total cost): {result['best_score']:.2f}")
    print(f"Number of fragments: {problem.n}")
    print(f"Fragment IDs order: {result['best_solution']}")

    
    # Plotting of the results
    

    # 1. Cost convergence
    plt.figure()
    plt.plot(result['history'])
    plt.xlabel("Generation")
    plt.ylabel("Best cost")
    plt.title("GA Convergence")
    plt.grid(True)

    # 2. Contigs vs iteration
    plt.figure()
    plt.plot(result['contigs_history'])
    plt.xlabel("Generation")
    plt.ylabel("Number of contigs")
    plt.title("Contigs vs Iteration")
    plt.grid(True)

    # 3. Overlap vs iteration
    plt.figure()
    plt.plot(result['overlap_history'])
    plt.xlabel("Generation")
    plt.ylabel("Total overlap")
    plt.title("Overlap vs Iteration")
    plt.grid(True)

    plt.show()
