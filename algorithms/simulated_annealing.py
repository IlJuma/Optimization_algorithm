from __future__ import annotations

import numpy as np

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem

# =========================
# NEIGHBORHOOD GENERATION
# =========================

def propose_neighbor(solution, rng, move_type=None):
    """
    Generate a neighboring permutation using one of three moves:
    - 'swap': exchange two reads
    - 'insert': remove and reinsert at different position
    - 'reverse': reverse a subsequence
    """
    neighbor = solution.copy()
    n = len(neighbor)
    
    if move_type is None:
        move_type = rng.choice(['swap', 'insert', 'reverse']) # ne sceglie uno casuale
    
    if move_type == 'swap':
        # Swap two random positions
        i, j = rng.choice(n, 2, replace=False)
        neighbor[i], neighbor[j] = neighbor[j], neighbor[i]
    
    elif move_type == 'insert':
        # Remove element at random position and insert at another
        i = rng.integers(0, n)
        j = rng.integers(0, n)
        element = neighbor.pop(i)
        neighbor.insert(j, element)
    
    elif move_type == 'reverse':
        # Reverse a random subsequence
        i = rng.integers(0, n)
        j = rng.integers(i + 1, n + 1)
        neighbor[i:j] = neighbor[i:j][::-1]
    
    return neighbor


# =========================
# ACCEPTANCE LOGIC
# =========================

def acceptance_probability(delta, temperature):
    """
    Calculate acceptance probability for a move.
    P = exp(-Δ/T) where Δ = f(π') - f(π)
    """
    if delta <= 0:  # Always accept improvements
        return 1.0

    if temperature <= 0:
        return 0.0

    return float(np.exp(-delta / temperature))


def cooling_schedule(T0, alpha, iteration):
    """
    Exponential cooling: T(t) = T0 * α^t
    """
    return T0 * (alpha ** iteration)

# =========================
# MAIN SIMULATED ANNEALING
# =========================

def simulated_annealing(
    problem=None,
    fragments=None,
    T0=100.0,
    alpha=0.995,
    max_iterations=10000,
    pair_indices=None,
    seed=42,
    verbose=True,
    log_interval=500,
    early_stop_temperature=1e-6,
):
    """
    Simulated Annealing algorithm for fragment ordering.
    
    Parameters:
    -----------
    fragments : list of str
        DNA sequences to order
    T0 : float
        Initial temperature
    alpha : float
        Cooling rate (0 < alpha < 1)
    max_iterations : int
        Maximum number of iterations
    pair_indices : dict
        Mapping of paired reads (optional)
    seed : int
        Random seed for reproducibility
    verbose : bool
        Print progress
    
    Returns:
    --------
    dict with keys:
        - 'best_solution': best permutation found
        - 'best_cost': cost of best solution
        - 'cost_history': list of costs over iterations
        - 'temperature_history': list of temperatures
        - 'acceptance_rate_history': acceptance rates per iteration
    """
    if problem is None:
        if fragments is None:
            fragments = load_fragments()
        problem = AssemblyProblem(fragments)
    
    rng = np.random.default_rng(seed)
    
    # Initialize with random permutation
    current_solution = problem.fragment_ids.copy()
    rng.shuffle(current_solution)
    current_cost = problem.evaluate(current_solution)

    best_solution = current_solution.copy()
    best_cost = current_cost
    
    # Tracking
    cost_history = [current_cost]
    temperature_history = [T0]
    acceptance_rate_history = []

    if pair_indices is None:
        pair_indices = {}
    
    # Main loop
    for iteration in range(max_iterations):
        temperature = cooling_schedule(T0, alpha, iteration)
        
        # Generate neighbor
        neighbor = propose_neighbor(current_solution, rng)
        neighbor_cost = problem.evaluate(neighbor)
        
        # Accept/reject decision
        delta = neighbor_cost - current_cost
        
        if acceptance_probability(delta, temperature) > rng.random():
            # Accept neighbor
            current_solution = neighbor
            current_cost = neighbor_cost
            accepted = True
        else:
            accepted = False
        
        # Track best solution
        if current_cost < best_cost:
            best_solution = current_solution.copy()
            best_cost = current_cost
            if verbose and iteration % 1000 == 0:
                print(f"Iter {iteration}: new best cost = {best_cost:.1f}")
        
        # Record history
        cost_history.append(current_cost)
        temperature_history.append(temperature)
        acceptance_rate_history.append(1.0 if accepted else 0.0)

        if verbose and (
            iteration == 0
            or (iteration + 1) % log_interval == 0
            or iteration == max_iterations - 1
        ):
            recent_window = acceptance_rate_history[max(0, len(acceptance_rate_history) - 100):]
            recent_acceptance = float(np.mean(recent_window)) if recent_window else 0.0

            print(
                f"[iter {iteration + 1:6d}] "
                f"T={temperature:10.4f} "
                f"current={current_cost:12.2f} "
                f"best={best_cost:12.2f} "
                f"accept_rate={recent_acceptance:6.3f}"
            )
        
        # Early stopping if temperature too low
        if temperature < early_stop_temperature:
            if verbose:
                print(f"Early stopping at iteration {iteration} (T < {early_stop_temperature})")
            break
    
    # Smooth acceptance rate for visualization
    acceptance_rate_smoothed = []
    window = min(100, max(1, max_iterations // 50))
    for i in range(len(acceptance_rate_history)):
        start = max(0, i - window // 2)
        end = min(len(acceptance_rate_history), i + window // 2)
        avg_rate = float(np.mean(acceptance_rate_history[start:end]))
        acceptance_rate_smoothed.append(avg_rate)

    best_breaks = problem.count_breaks(best_solution)
    best_contigs = problem.count_contigs(best_solution)
    best_overlap = problem.total_overlap(best_solution)

    if verbose:
        print(f"\nSimulated Annealing completed!")
        print(f"Initial cost: {cost_history[0]:.1f}")
        print(f"Final best cost: {best_cost:.1f}")

        if cost_history[0] != 0:
            improvement = cost_history[0] - best_cost
            relative_improvement = 100 * improvement / abs(cost_history[0])
            print(f"Improvement: {improvement:.1f} ({relative_improvement:.1f}%)")
        else:
            print("Improvement: initial cost was 0")

    return {
        "best_solution": best_solution,
        "best_cost": best_cost,
        "cost_history": cost_history,
        "temperature_history": temperature_history,
        "acceptance_rate_history": acceptance_rate_smoothed,
        "current_solution": current_solution,
        "current_cost": current_cost,
        "best_breaks": best_breaks,
        "best_contigs": best_contigs,
        "best_total_overlap": best_overlap,
        "seed": seed,
        "T0": T0,
        "alpha": alpha,
        "max_iterations": max_iterations,
    }


def load_problem():
    fragments = load_fragments()
    return AssemblyProblem(fragments)

import time

def optimize(problem, config, rng):
    """
    Standardized wrapper to integrate Simulated Annealing into the experiment runner.
    Translates the shared config into SA parameters and formats the output.
    """
    start_time = time.time()
    
    # 1. Extract shared budget and SA-specific config
    max_evals = getattr(config, 'MAX_EVALUATIONS', 10000)
    t0 = getattr(config, 'SA_INITIAL_TEMP', 10.0)
    alpha = getattr(config, 'SA_COOLING_RATE', 0.995)
    min_temp = getattr(config, 'SA_MIN_TEMP', 1e-3)
    
    # Generate a numpy-safe seed from the shared random generator
    numpy_seed = rng.randint(0, 2**32 - 1)
    
    # 2. Call Person C's core algorithm
    raw_result = simulated_annealing(
        problem=problem,
        T0=t0,
        alpha=alpha,
        max_iterations=max_evals,
        seed=numpy_seed,
        early_stop_temperature=min_temp,
        verbose=False  # Turn off their internal prints so the runner stays clean
    )
    
    runtime = time.time() - start_time
    
    # 3. Format the history arrays
    # Person C's 'cost_history' tracks the noisy current cost. 
    # We need to calculate the "best score so far" history for the generic plotting script.
    best_history = []
    current_best = float('inf')
    for cost in raw_result['cost_history']:
        if cost < current_best:
            current_best = cost
        best_history.append(current_best)
        
    # 4. Return the globally agreed-upon Result dictionary
    return {
        "method": "Simulated Annealing",
        "best_solution": raw_result["best_solution"],
        "best_score": raw_result["best_cost"],
        "history": best_history,                             # Smooth best-cost staircase
        "current_cost_history": raw_result["cost_history"],  # Noisy exploration cloud
        "evaluations": len(raw_result["cost_history"]) - 1,  # Number of actual steps taken
        "runtime_sec": runtime,
        "breaks": raw_result["best_breaks"]
    }

if __name__ == "__main__":
    problem = load_problem()

    result = simulated_annealing(problem=problem)

    print("\n" + "=" * 70)
    print("BEST SOLUTION FOUND")
    print("=" * 70)
    print(f"Fragments: {problem.n}")
    print(f"Best cost: {result['best_cost']:.2f}")
    print(f"Breaks: {result['best_breaks']}")
    print(f"Contigs: {result['best_contigs']}")
    print(f"Total overlap: {result['best_total_overlap']}")
    print("Best solution (fragment IDs):")
    print(result["best_solution"])