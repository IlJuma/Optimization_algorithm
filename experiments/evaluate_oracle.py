import json
import os
import sys

# Dynamically add the project root to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem
from model import config

def main():
    print("Evaluating the Oracle (True Genomic Sequence)...")
    
    fragments = load_fragments()
    problem = AssemblyProblem(fragments=fragments, min_overlap=config.MIN_OVERLAP)

    # The Oracle sequence is exactly the fragments sorted by their actual starting position
    sorted_fragments = sorted(fragments, key=lambda f: f.frag_start)
    oracle_solution = [f.fragment_id for f in sorted_fragments]

    # Evaluate the perfect sequence
    oracle_cost = problem.evaluate(oracle_solution)
    oracle_contigs = problem.count_contigs(oracle_solution)
    oracle_overlap = problem.total_overlap(oracle_solution)
    oracle_breaks = problem.count_breaks(oracle_solution)

    oracle_data = {
        "method": "Oracle (True Genome)",
        "best_score": oracle_cost,
        "contigs": oracle_contigs,
        "overlap": oracle_overlap,
        "breaks": oracle_breaks
    }

    os.makedirs("reports", exist_ok=True)
    with open("reports/oracle_baseline.json", "w") as f:
        json.dump(oracle_data, f, indent=4)

    print("\n--- ORACLE METRICS ---")
    print(f"Cost:    {oracle_cost:.2f}")
    print(f"Contigs: {oracle_contigs}")
    print(f"Overlap: {oracle_overlap}")
    print(f"Breaks:  {oracle_breaks}")
    print("\nSaved to reports/oracle_baseline.json")

if __name__ == "__main__":
    main()