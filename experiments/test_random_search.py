import random
import os
import sys

# Dynamically add the project root to the path so imports work smoothly
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the actual Random Search algorithm
from algorithms import random_search

# Import the real global config
from model import config

# 1. Import the REAL loaders and problem class
from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem

# ---------------------------------------------------------
# Run Algorithm
# ---------------------------------------------------------
if __name__ == "__main__":
    rng = random.Random(42)
    
    # Temporarily override the real config parameters for a quick test
    config.MAX_EVALUATIONS = 5000
    config.MAX_TIME_SEC = 5.0
    
    # 2. Load the real data and build the real problem
    print("Loading real fragments...")
    fragments = load_fragments() 
    
    print("Building real AssemblyProblem...")
    problem = AssemblyProblem(fragments=fragments, min_overlap=config.MIN_OVERLAP)
    
    print(f"Starting Random Search on {problem.n} real fragments...")
    result = random_search.optimize(problem, config, rng)
    
    print("\n--- Results ---")
    print(f"Method:        {result['method']}")
    print(f"Evaluations:   {result['evaluations']}")
    print(f"Runtime:       {result['runtime_sec']:.4f} seconds")
    print(f"Best Score:    {result['best_score']}")
    print(f"Breaks:        {result.get('breaks', 'N/A')}")
    # Note: We skip printing 'best_solution' here because printing thousands of IDs will flood your terminal!