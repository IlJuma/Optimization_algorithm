import time
import numpy as np

def calculate_fitness(solution, problem):
    return problem.evaluate(solution)

def crossover(parent1, parent2, rng, crossover_rate):
    if rng.random() < crossover_rate:
        n = len(parent1)
        point = rng.integers(1, n)
        child1 = parent1[:point] + [f for f in parent2 if f not in parent1[:point]]
        child2 = parent2[:point] + [f for f in parent1 if f not in parent2[:point]]
        return child1, child2
    else:
        return parent1.copy(), parent2.copy()

def mutate(solution, rng, mutation_rate):
    n = len(solution)
    mutated = solution.copy()
    for i in range(n):
        if rng.random() < mutation_rate:
            j = rng.integers(0, n-1)
            if j >= i:
                j += 1
            mutated[i], mutated[j] = mutated[j], mutated[i]
    return mutated

def optimize(problem, config, shared_rng):
    """
    Standardized wrapper for the Experiment Runner.
    """
    start_time = time.time()
    
    # Extract config parameters (with fallbacks just in case)
    pop_size = getattr(config, 'GA_POP_SIZE', 50)
    crossover_rate = getattr(config, 'GA_CROSSOVER_RATE', 0.8)
    mutation_rate = getattr(config, 'GA_MUTATION_RATE', 0.05)
    elitism = getattr(config, 'GA_ELITISM', 2)
    max_evals = getattr(config, 'MAX_EVALUATIONS', 10000)
    
    # Create a numpy random generator from the shared random seed to ensure reproducibility
    numpy_seed = shared_rng.randint(0, 2**32 - 1)
    rng = np.random.default_rng(numpy_seed)

    # Initialize Population
    population = [rng.permutation(problem.fragment_ids).tolist() for _ in range(pop_size)]

    best_solution = None
    best_score = float('inf')

    history = []
    current_cost_history = []
    contigs_history = []
    overlap_history = []
    
    evaluations = 0

    while evaluations < max_evals:
        # Evaluate population
        fitness_scores = []
        for sol in population:
            if evaluations >= max_evals:
                break
            score = calculate_fitness(sol, problem)
            fitness_scores.append(score)
            current_cost_history.append(score) # Track every single guess for the Violin Plot
            evaluations += 1
            
        fitness_array = np.array(fitness_scores)

        # Update best global solution
        current_best_idx = np.argmin(fitness_array)
        if fitness_array[current_best_idx] < best_score:
            best_score = fitness_array[current_best_idx]
            best_solution = population[current_best_idx].copy()

        # Record metrics for the *best solution found so far*
        history.append(best_score)
        contigs_history.append(problem.count_contigs(best_solution))
        overlap_history.append(problem.total_overlap(best_solution))

        if evaluations >= max_evals:
            break

        # Selection (Roulette Wheel)
        inv_fitness = 1 / (fitness_array - fitness_array.min() + 1e-6)
        probs = inv_fitness / inv_fitness.sum()
        selected_indices = rng.choice(range(len(population)), size=pop_size, replace=True, p=probs)

        # Crossover
        offspring = []
        for i in range(0, pop_size - 1, 2):
            p1 = population[selected_indices[i]]
            p2 = population[selected_indices[i+1]]
            c1, c2 = crossover(p1, p2, rng, crossover_rate)
            offspring.extend([c1, c2])

        # Handle odd population sizes
        if len(offspring) < pop_size:
            offspring.append(population[selected_indices[-1]].copy())

        # Mutation
        offspring = [mutate(child, rng, mutation_rate) for child in offspring]

        # Survivor selection & Elitism
        combined = population + offspring
        combined.sort(key=lambda sol: calculate_fitness(sol, problem))
        
        if elitism > 0:
            population = combined[:elitism] + offspring[:-elitism]
        else:
            population = offspring[:pop_size]

    runtime = time.time() - start_time
    final_breaks = problem.count_breaks(best_solution)

    return {
        "method": "Genetic Algorithm",
        "best_solution": best_solution,
        "best_score": float(best_score),                                  # Cast to normal float for JSON serialization
        "history": [float(x) for x in history],                           
        "current_cost_history": [float(x) for x in current_cost_history], 
        "contigs_history": [int(x) for x in contigs_history],             
        "overlap_history": [int(x) for x in overlap_history],             
        "evaluations": int(evaluations),
        "runtime_sec": float(runtime),
        "breaks": int(final_breaks)
    }