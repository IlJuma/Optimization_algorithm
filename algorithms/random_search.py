import time

def optimize(problem, config, rng):
    """
    Random Search Optimization Algorithm.
    
    Explores the solution space by sampling random permutations of reads 
    and tracking the best one found. Serves as a baseline method.
    """
    start_time = time.time()
    
    best_solution = None
    best_score = float('inf')
    history = []
    evaluations = 0
    
    # extract budgets from config, falling back to defaults if missing
    max_evals = getattr(config, 'MAX_EVALUATIONS', 10000)
    max_time_sec = getattr(config, 'MAX_TIME_SEC', 60)
    
    while evaluations < max_evals:
        # check time limit
        if time.time() - start_time > max_time_sec:
            break
            
        # sample a new random permutation
        # Person A's problem model to generate a valid random permutation
        current_solution = problem.random_solution(rng)
        
        # evaluate the permutation
        score = problem.evaluate(current_solution)
        evaluations += 1
        
        # update the best solution if the current one is better
        if score < best_score:
            best_score = score
            best_solution = current_solution.copy()
            
        # track history (tracking every step allows for convergence plotting later)
        history.append(best_score)
        
    runtime = time.time() - start_time
    
    # Return results in the global format expected by analyze_results.py
    return {
        "method": "Random Search",
        "best_solution": best_solution,
        "best_score": best_score,
        "history": history,
        "evaluations": evaluations,
        "runtime_sec": runtime
    }