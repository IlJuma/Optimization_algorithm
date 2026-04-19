import time
import numpy as np
import random
from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem

def base_fragment_id(oriented_fragment_id):
    """Strips the _F or _R suffix to get the core fragment ID."""
    if oriented_fragment_id.endswith("_F") or oriented_fragment_id.endswith("_R"):
        return oriented_fragment_id[:-2]
    raise ValueError(f"Invalid oriented fragment id: {oriented_fragment_id}")

def build_random_oriented_solution(problem, rng):
    base_fragment_ids = problem.base_fragment_ids.copy()
    rng.shuffle(base_fragment_ids)
    return [f"{fid}_{rng.choice(['F', 'R'])}" for fid in base_fragment_ids]

def calculate_fitness(solution, problem):
    return problem.evaluate(solution)

def crossover(parent1, parent2, rng, crossover_rate):
    """Oriented Order Crossover (OX1): Preserves both sequence order AND orientation."""
    if rng.random() > crossover_rate:
        return parent1.copy(), parent2.copy()
        
    n = len(parent1)
    
    def ox1(p1, p2):
        start, end = np.sort(rng.choice(n, 2, replace=False))
        child = [None] * n
        
        # Copy the swath from parent 1 (keeps parent 1's orientation)
        child[start:end+1] = p1[start:end+1]
        
        # Track which base fragments have been used
        used_bases = {base_fragment_id(x) for x in child[start:end+1]}
        
        # Filter parent 2 for unused base fragments (keeps parent 2's orientation)
        p2_filtered = [x for x in p2 if base_fragment_id(x) not in used_bases]
        
        # Fill remaining slots
        p2_idx = 0
        for i in range(n):
            if child[i] is None:
                child[i] = p2_filtered[p2_idx]
                p2_idx += 1
        return child

    return ox1(parent1, parent2), ox1(parent2, parent1)

def mutate(solution, rng, mutation_rate):
    """Oriented Mutation: Sequence Reversal AND Orientation Flipping."""
    mutated = solution.copy()
    
    if rng.random() < mutation_rate:
        if rng.random() < 0.5:
            # Move 1: Biological Inversion (Reverse chunk AND flip orientations)
            n = len(solution)
            i, j = rng.choice(n, 2, replace=False)
            if i > j:
                i, j = j, i
            
            chunk = mutated[i:j+1][::-1]
            flipped_chunk = []
            for frag in chunk:
                if frag.endswith("_F"):
                    flipped_chunk.append(frag[:-2] + "_R")
                else:
                    flipped_chunk.append(frag[:-2] + "_F")
            mutated[i:j+1] = flipped_chunk
        else:
            # Move 2: Point Mutation (Flip one single fragment's orientation)
            idx = rng.integers(0, len(mutated))
            frag = mutated[idx]
            if frag.endswith("_F"):
                mutated[idx] = frag[:-2] + "_R"
            else:
                mutated[idx] = frag[:-2] + "_F"
                
    return mutated

def optimize(problem=None, config=None, shared_rng=None):
    """Standardized wrapper for the Experiment Runner."""
    if problem is None:
        fragments = load_fragments()
        problem = AssemblyProblem(fragments)

    # Safely convert standard Python random to numpy random
    if isinstance(shared_rng, random.Random):
        numpy_seed = shared_rng.randint(0, 2**32 - 1)
        rng = np.random.default_rng(numpy_seed)
    elif shared_rng is None:
        rng = np.random.default_rng(42)
    else:
        rng = shared_rng

    start_time = time.time()
    
    # Extract config parameters
    pop_size = getattr(config, 'GA_POP_SIZE', 50)
    crossover_rate = getattr(config, 'GA_CROSSOVER_RATE', 0.8)
    mutation_rate = getattr(config, 'GA_MUTATION_RATE', 0.2) 
    elitism = getattr(config, 'GA_ELITISM', 2)
    max_evals = getattr(config, 'MAX_EVALUATIONS', 10000)
    max_time_sec = getattr(config, 'MAX_TIME_SEC', 60)

    # Initialize Population with proper F/R orientations
    population = [build_random_oriented_solution(problem, rng) for _ in range(pop_size)]
    fitness_scores = [] # Track scores to prevent duplicate evaluations

    best_solution = None
    best_score = float('inf')

    history = []
    current_cost_history = []
    contigs_history = []
    overlap_history = []
    
    evaluations = 0

    # 1. Initial Evaluation of the starting population
    for sol in population:
        if evaluations >= max_evals or time.time() - start_time > max_time_sec:
            break
        score = calculate_fitness(sol, problem)
        fitness_scores.append(score)
        current_cost_history.append(score) 
        
        if score < best_score:
            best_score = score
            best_solution = sol.copy()

        history.append(best_score)
        contigs_history.append(problem.count_contigs(best_solution))
        overlap_history.append(problem.total_overlap(best_solution))
        evaluations += 1

    # 2. Main Generation Loop
    while evaluations < max_evals and (time.time() - start_time) <= max_time_sec:
            
        fitness_array = np.array(fitness_scores)

        # Selection (Tournament)
        selected_indices = []
        tournament_size = 3
        for _ in range(pop_size):
            fighters = rng.choice(range(len(population)), size=tournament_size, replace=False)
            winner = fighters[np.argmin(fitness_array[fighters])]
            selected_indices.append(winner)

        # Crossover
        offspring = []
        for i in range(0, pop_size - 1, 2):
            p1 = population[selected_indices[i]]
            p2 = population[selected_indices[i+1]]
            c1, c2 = crossover(p1, p2, rng, crossover_rate)
            offspring.extend([c1, c2])

        if len(offspring) < pop_size:
            offspring.append(population[selected_indices[-1]].copy())

        # Mutation
        offspring = [mutate(child, rng, mutation_rate) for child in offspring]

        # Evaluate ONLY the new offspring (Massive speed boost!)
        offspring_scores = []
        for sol in offspring:
            if evaluations >= max_evals or time.time() - start_time > max_time_sec:
                break
            score = calculate_fitness(sol, problem)
            offspring_scores.append(score)
            current_cost_history.append(score)
            
            if score < best_score:
                best_score = score
                best_solution = sol.copy()

            history.append(best_score)
            contigs_history.append(problem.count_contigs(best_solution))
            overlap_history.append(problem.total_overlap(best_solution))
            evaluations += 1

        if len(offspring_scores) == 0:
            break

        # Survivor selection & Elitism (Combine parents and ONLY the evaluated offspring)
        combined = population + offspring[:len(offspring_scores)]
        combined_scores = fitness_scores + offspring_scores
        
        ranked = sorted(zip(combined, combined_scores), key=lambda x: x[1])
        population = [sol for sol, _ in ranked[:pop_size]]
        fitness_scores = [score for _, score in ranked[:pop_size]]
        
        if elitism > 0:
            elites = [sol for sol, _ in ranked[:elitism]]
            elite_scores = [score for _, score in ranked[:elitism]]
            
            population[-elitism:] = elites
            fitness_scores[-elitism:] = elite_scores

    runtime = time.time() - start_time
    final_breaks = problem.count_breaks(best_solution)

    return {
        "method": "Genetic Algorithm",
        "best_solution": best_solution,
        "best_score": float(best_score),
        "history": [float(x) for x in history],                           
        "current_cost_history": [float(x) for x in current_cost_history], 
        "contigs_history": [int(x) for x in contigs_history],             
        "overlap_history": [int(x) for x in overlap_history],             
        "evaluations": int(evaluations),
        "runtime_sec": float(runtime),
        "breaks": int(final_breaks)
    }