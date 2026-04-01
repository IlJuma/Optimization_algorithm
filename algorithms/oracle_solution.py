"""
Oracle solution module.

Provides a ground-truth baseline for fragment ordering based on known
genomic coordinates from the simulation pipeline.

This acts as a reference solution against which optimization algorithms
(random search, simulated annealing, GA) can be compared.
"""

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem


# =========================
# CORE LOGIC
# =========================

def build_oracle_permutation(fragments):
    """
    Construct the oracle permutation using true genomic positions.

    Fragments are ordered by:
        1. frag_start (primary)
        2. frag_end (secondary)
        3. fragment_id (tie-breaker)

    Returns:
        List[str] : ordered fragment_ids
    """
    ordered = sorted(
        fragments,
        key=lambda f: (f.frag_start, f.frag_end, f.fragment_id)
    )
    return [f.fragment_id for f in ordered]


def evaluate_oracle(problem, fragments):
    """
    Evaluate the oracle solution using the AssemblyProblem scoring model.

    Returns:
        dict compatible with other algorithm outputs
    """
    solution = build_oracle_permutation(fragments)

    score = problem.evaluate(solution)
    breaks = problem.count_breaks(solution)
    contigs = problem.count_contigs(solution)
    total_overlap = problem.total_overlap(solution)

    return {
        "method": "Oracle Ground Truth",
        "best_solution": solution,
        "best_score": score,
        "history": [score],          # single evaluation
        "evaluations": 1,
        "runtime_sec": 0.0,
        "best_breaks": breaks,
        "best_contigs": contigs,
        "best_total_overlap": total_overlap,
    }


# =========================
# CONVENIENCE ENTRY POINT
# =========================

def oracle_solution(problem=None, fragments=None):
    """
    Convenience wrapper so run_experiments.py can call this in one line.

    Usage:
        result = oracle_solution(problem=problem)

    If fragments are not provided, they will be loaded automatically.
    """
    if fragments is None:
        fragments = load_fragments()

    if problem is None:
        problem = AssemblyProblem(fragments)

    return evaluate_oracle(problem, fragments)


# =========================
# CLI USAGE
# =========================

if __name__ == "__main__":
    fragments = load_fragments()
    problem = AssemblyProblem(fragments)

    result = oracle_solution(problem=problem, fragments=fragments)

    print("Oracle evaluation:")
    print(f"Fragments: {problem.n}")
    print(f"Score: {result['best_score']:.2f}")
    print(f"Breaks: {result['best_breaks']}")
    print(f"Contigs: {result['best_contigs']}")
    print(f"Total overlap: {result['best_total_overlap']}")