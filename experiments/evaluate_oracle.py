import json
import os
import sys

# Dynamically add the project root to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem
from model import config
from algorithms.oracle_solution import oracle_solution


def main():
    print("Evaluating the Oracle (True Genomic Sequence)...")

    fragments = load_fragments()
    problem = AssemblyProblem(
        fragments=fragments,
        min_overlap=config.MIN_OVERLAP,
        break_factor=config.BREAK_FACTOR,
        partial_overlap_factor=config.PARTIAL_OVERLAP_FACTOR,
        dense_threshold=config.DENSE_THRESHOLD,
        max_neighbors_per_fragment=config.MAX_NEIGHBORS_PER_FRAGMENT,
    )

    result = oracle_solution(problem=problem, fragments=fragments)

    oracle_data = {
        "method": result["method"],
        "best_solution": result["best_solution"],
        "best_score": result["best_score"],
        "history": result["history"],
        "evaluations": result["evaluations"],
        "runtime_sec": result["runtime_sec"],
        "best_breaks": result["best_breaks"],
        "best_contigs": result["best_contigs"],
        "best_total_overlap": result["best_total_overlap"],
        "breaks": result["breaks"],
        "contigs": result["contigs"],
        "total_overlap": result["total_overlap"],
    }

    os.makedirs("reports", exist_ok=True)
    out_path = "reports/oracle_baseline.json"
    with open(out_path, "w") as f:
        json.dump(oracle_data, f, indent=4)

    print("\n--- ORACLE METRICS ---")
    print(f"Cost:    {oracle_data['best_score']:.2f}")
    print(f"Contigs: {oracle_data['best_contigs']}")
    print(f"Overlap: {oracle_data['best_total_overlap']}")
    print(f"Breaks:  {oracle_data['best_breaks']}")
    print("\nOracle solution uses oriented fragment IDs:")
    print("  <fragment_id>_F = stored fragment sequence from fragments.fasta")
    print("  <fragment_id>_R = reverse complement of stored fragment sequence")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()