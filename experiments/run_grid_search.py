import os
import sys
import time
import random
import itertools
import pandas as pd
from types import SimpleNamespace

# Dynamically add the project root to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem
from model import config as default_config
from algorithms import simulated_annealing
from algorithms import genetic_algorithm

def run_grid_search(problem, base_config, rng, algo_name, algo_func, param_grid):
    """Generates combinations, runs the algorithm, and returns a DataFrame of results."""
    
    # 1. Generate all possible combinations from the dictionary
    keys = param_grid.keys()
    values = param_grid.values()
    combinations = list(itertools.product(*values))
    
    print(f"\nStarting Grid Search for {algo_name}: {len(combinations)} combinations.")
    
    results = []
    
    for i, combo in enumerate(combinations):
        # Create a temporary config combining defaults with this specific combo
        combo_dict = dict(zip(keys, combo))
        temp_config_dict = {k: getattr(base_config, k) for k in dir(base_config) if not k.startswith('__')}
        temp_config_dict.update(combo_dict)
        
        # Convert dictionary back to an object so getattr() works inside the algorithms
        search_config = SimpleNamespace(**temp_config_dict)
        
        print(f"[{i+1}/{len(combinations)}] Testing: {combo_dict} ... ", end="")
        
        # Run the algorithm
        start_time = time.time()
        res = algo_func(problem=problem, config=search_config, rng=rng)
        runtime = time.time() - start_time
        
        print(f"Score: {res['best_score']:.1f}")
        
        # Save the data
        row = combo_dict.copy()
        row["best_score"] = res["best_score"]
        row["breaks"] = res.get("breaks", "N/A")
        row["evaluations"] = res["evaluations"]
        row["runtime_sec"] = round(runtime, 2)
        results.append(row)
        
    return pd.DataFrame(results)


def main():
    print("Loading fragments and building problem...")
    fragments = load_fragments()
    problem = AssemblyProblem(fragments=fragments, min_overlap=default_config.MIN_OVERLAP)
    rng = random.Random(42)

    os.makedirs("reports", exist_ok=True)

    #SA grid search
    sa_grid = {
        "SA_INITIAL_TEMP": [100.0, 500.0],
        "SA_COOLING_RATE": [0.999, 0.9995, 0.9999], # Push cooling to be extremely slow
        "SA_MIN_TEMP": [1e-4],                      # Drop the floor lower
        "MAX_EVALUATIONS": [15000]                  # Raise the ceiling
    }
    
    sa_df = run_grid_search(problem, default_config, rng, "Simulated Annealing", simulated_annealing.optimize, sa_grid)
    
    # Sort by best score and save
    sa_df = sa_df.sort_values(by="best_score")
    sa_out = "reports/grid_search_sa_phase2.csv"
    sa_df.to_csv(sa_out, index=False)
    print(f"SA Grid Search complete. Saved to {sa_out}")


    # GA grid search
    ga_grid = {
        "GA_POP_SIZE": [30, 50],             # Smaller populations = more generations
        "GA_MUTATION_RATE": [0.3, 0.4, 0.5], # Push the Inversion mutation to the limit
        "GA_CROSSOVER_RATE": [0.8, 0.9],     # Test slightly more aggressive crossover
        "GA_ELITISM": [1, 2],                # Keep elitism low to prevent stagnation
        "MAX_EVALUATIONS": [5000]
    }
    
    ga_df = run_grid_search(problem, default_config, rng, "Genetic Algorithm", genetic_algorithm.optimize, ga_grid)
    
    # Sort by best score and save
    ga_df = ga_df.sort_values(by="best_score")
    ga_out = "reports/grid_search_ga_phase2.csv"
    ga_df.to_csv(ga_out, index=False)
    print(f"GA Grid Search complete. Saved to {ga_out}")

if __name__ == "__main__":
    main()