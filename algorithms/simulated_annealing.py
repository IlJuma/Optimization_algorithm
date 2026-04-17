from __future__ import annotations

import time
import numpy as np

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem


# =========================
# ORIENTATION HELPERS
# =========================

def flip_orientation(oriented_fragment_id):
    if oriented_fragment_id.endswith("_F"):
        return oriented_fragment_id[:-2] + "_R"
    if oriented_fragment_id.endswith("_R"):
        return oriented_fragment_id[:-2] + "_F"
    raise ValueError(f"Invalid oriented fragment id: {oriented_fragment_id}")


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
        orientation = rng.choice(["F", "R"])
        solution.append(f"{fragment_id}_{orientation}")

    return solution


# =========================
# NEIGHBORHOOD GENERATION
# =========================

def propose_neighbor(solution, rng, move_type=None):
    """
    Generate a neighboring permutation using one of four moves:
    - 'swap': exchange two fragments
    - 'insert': remove and reinsert at different position
    - 'reverse': reverse a subsequence
    - 'flip': flip the orientation of one fragment
    """
    neighbor = solution.copy()
    n = len(neighbor)

    if n < 2:
        if n == 1 and (move_type == "flip" or move_type is None):
            neighbor[0] = flip_orientation(neighbor[0])
        return neighbor

    if move_type is None:
        move_type = rng.choice(["swap", "insert", "reverse", "flip"])

    if move_type == "swap":
        i, j = rng.choice(n, 2, replace=False)
        neighbor[i], neighbor[j] = neighbor[j], neighbor[i]

    elif move_type == "insert":
        i = rng.integers(0, n)
        j = rng.integers(0, n)
        element = neighbor.pop(i)
        neighbor.insert(j, element)

    elif move_type == "reverse":
        i = rng.integers(0, n)
        j = rng.integers(i + 1, n + 1)
        neighbor[i:j] = neighbor[i:j][::-1]

    elif move_type == "flip":
        i = rng.integers(0, n)
        neighbor[i] = flip_orientation(neighbor[i])

    else:
        raise ValueError(f"Unknown move_type: {move_type}")

    return neighbor


# =========================
# ACCEPTANCE LOGIC
# =========================

def acceptance_probability(delta, temperature):
    """
    Calculate acceptance probability for a move.
    P = exp(-Δ/T) where Δ = f(π') - f(π)
    """
    if delta <= 0:
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

def optimize(
    problem=None,
    fragments=None,
    config=None,
    rng=None,
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
    Simulated Annealing algorithm for oriented fragment ordering.
    """
    if problem is None:
        if fragments is None:
            fragments = load_fragments()
        problem = AssemblyProblem(fragments)

    if config is not None:
        T0 = getattr(config, "SA_INITIAL_TEMP", T0)
        alpha = getattr(config, "SA_COOLING_RATE", alpha)
        max_iterations = getattr(config, "MAX_EVALUATIONS", max_iterations)
        early_stop_temperature = getattr(config, "SA_MIN_TEMP", early_stop_temperature)
        seed = getattr(config, "RANDOM_SEED", seed)
        max_time_sec = getattr(config, "MAX_TIME_SEC", None)
    else:
        max_time_sec = None

    if rng is None:
        rng = np.random.default_rng(seed)

    start_time = time.time()

    # Initialize with random oriented permutation
    current_solution = build_random_oriented_solution(problem, rng)
    current_cost = problem.evaluate(current_solution)

    best_solution = current_solution.copy()
    best_cost = current_cost

    # Tracking
    cost_history = [current_cost]
    temperature_history = [T0]
    acceptance_rate_history = []
    current_cost_history = [current_cost]
    evaluations = 1

    if pair_indices is None:
        pair_indices = {}

    # Main loop
    for iteration in range(max_iterations):
        if max_time_sec is not None and (time.time() - start_time > max_time_sec):
            break

        temperature = cooling_schedule(T0, alpha, iteration)

        # Generate neighbor
        neighbor = propose_neighbor(current_solution, rng)
        neighbor_cost = problem.evaluate(neighbor)
        evaluations += 1

        # Accept/reject decision
        delta = neighbor_cost - current_cost

        if acceptance_probability(delta, temperature) > rng.random():
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
        current_cost_history.append(current_cost)
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
    runtime = time.time() - start_time

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
        "method": "Simulated Annealing",
        "best_solution": best_solution,
        "best_score": best_cost,
        "best_cost": best_cost,
        "history": cost_history,
        "cost_history": cost_history,
        "current_cost_history": current_cost_history,
        "temperature_history": temperature_history,
        "acceptance_rate_history": acceptance_rate_smoothed,
        "current_solution": current_solution,
        "current_cost": current_cost,
        "evaluations": evaluations,
        "runtime_sec": runtime,
        "best_breaks": best_breaks,
        "best_contigs": best_contigs,
        "best_total_overlap": best_overlap,
        "breaks": best_breaks,
        "contigs": best_contigs,
        "total_overlap": best_overlap,
        "seed": seed,
        "T0": T0,
        "alpha": alpha,
        "max_iterations": max_iterations,
    }


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
    return optimize(
        problem=problem,
        fragments=fragments,
        config=None,
        rng=None,
        T0=T0,
        alpha=alpha,
        max_iterations=max_iterations,
        pair_indices=pair_indices,
        seed=seed,
        verbose=verbose,
        log_interval=log_interval,
        early_stop_temperature=early_stop_temperature,
    )


def load_problem():
    fragments = load_fragments()
    return AssemblyProblem(fragments)


if __name__ == "__main__":
    problem = load_problem()

    result = optimize(problem=problem)

    print("\n" + "=" * 70)
    print("BEST SOLUTION FOUND")
    print("=" * 70)
    print(f"Fragments: {problem.n_fragments}")
    print(f"Best cost: {result['best_cost']:.2f}")
    print(f"Breaks: {result['best_breaks']}")
    print(f"Contigs: {result['best_contigs']}")
    print(f"Total overlap: {result['best_total_overlap']}")
    print("Best solution (fragment IDs):")
    print(result["best_solution"])