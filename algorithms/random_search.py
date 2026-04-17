import time
import random

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem


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


def optimize(problem=None, config=None, rng=None):
    """
    Random Search Optimization Algorithm.

    Explores the solution space by sampling random permutations of oriented
    fragments and tracking the best one found. Serves as a baseline method.
    """
    if problem is None:
        fragments = load_fragments()
        problem = AssemblyProblem(fragments)

    if rng is None:
        rng = random.Random(42)

    start_time = time.time()

    best_solution = None
    best_score = float("inf")
    history = []
    current_cost_history = []
    evaluations = 0

    # extract budgets from config, falling back to defaults if missing
    max_evals = getattr(config, "MAX_EVALUATIONS", 10000) if config is not None else 10000
    max_time_sec = getattr(config, "MAX_TIME_SEC", 60) if config is not None else 60

    while evaluations < max_evals:
        # check time limit
        if time.time() - start_time > max_time_sec:
            break

        # sample a new random valid oriented permutation
        current_solution = build_random_oriented_solution(problem, rng)

        # evaluate the permutation
        score = problem.evaluate(current_solution)
        evaluations += 1

        # update the best solution if the current one is better
        if score < best_score:
            best_score = score
            best_solution = current_solution.copy()

        # track history
        history.append(best_score)
        current_cost_history.append(score)

    runtime = time.time() - start_time

    final_breaks = problem.count_breaks(best_solution) if best_solution else 0
    final_contigs = problem.count_contigs(best_solution) if best_solution else 0
    final_total_overlap = problem.total_overlap(best_solution) if best_solution else 0

    return {
        "method": "Random Search",
        "best_solution": best_solution,
        "best_score": best_score,
        "history": history,
        "current_cost_history": current_cost_history,
        "evaluations": evaluations,
        "runtime_sec": runtime,
        "best_breaks": final_breaks,
        "best_contigs": final_contigs,
        "best_total_overlap": final_total_overlap,
        "breaks": final_breaks,
        "contigs": final_contigs,
        "total_overlap": final_total_overlap,
    }


def random_search(problem=None, config=None, rng=None):
    """
    Compatibility wrapper for notebook usage:
        from algorithms.random_search import random_search
    """
    return optimize(problem=problem, config=config, rng=rng)


if __name__ == "__main__":
    fragments = load_fragments()
    problem = AssemblyProblem(fragments)

    result = optimize(problem=problem)

    print("Random Search evaluation:")
    print(f"Fragments: {problem.n_fragments}")
    print(f"Best score: {result['best_score']:.2f}")
    print(f"Breaks: {result['best_breaks']}")
    print(f"Contigs: {result['best_contigs']}")
    print(f"Total overlap: {result['best_total_overlap']}")