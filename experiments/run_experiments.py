import json
import random
import time
import os

# Import the shared config and problem model
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from model import config
# from model.problem import AssemblyProblem  # Uncomment when ready

# Import the algorithms
from algorithms import random_search
# from algorithms import simulated_annealing # Uncomment when ready
# from algorithms import genetic_algorithm   # Uncomment when ready

def main():
    print("Initializing Genome Assembly Optimization Experiment...")
    rng = random.Random(config.RANDOM_SEED)
    
    # 1 Load the problem
    # problem = AssemblyProblem(reads_file="data/fastq/reads_R1.fastq", config=config)
    problem = None # Placeholder
    
    results = []
    
    # 2 Run Random Search
    print("\nRunning Baseline: Random Search...")
    rs_result = random_search.optimize(problem, config, rng)
    results.append(rs_result)
    
    # 3 Run Simulated Annealing
    # print("\nRunning Simulated Annealing...")
    # sa_result = simulated_annealing.optimize(problem, config, rng)
    # results.append(sa_result)
    
    # 4 Run Genetic Algorithm
    # print("\nRunning Genetic Algorithm...")
    # ga_result = genetic_algorithm.optimize(problem, config, rng)
    # results.append(ga_result)
    
    # Save results to JSON for analyze_results.py script
    os.makedirs("reports", exist_ok=True)
    out_file = "reports/experiment_results.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=4)
        
    print(f"\nExperiments done. Results saved to {out_file}")
    print("Run 'python experiments/analyze_results.py' to generate plots.")

if __name__ == "__main__":
    main()