import json
import random
import os
import sys
from types import SimpleNamespace

# Import the shared config and problem model
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from model import config
from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem


# Import the algorithms
from algorithms.oracle_solution import oracle_solution
from algorithms import random_search
from algorithms import simulated_annealing
from algorithms import genetic_algorithm


def _as_list(value):
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]

def main():
    print("Initializing Genome Assembly Optimization Experiment...")
    rng = random.Random(config.RANDOM_SEED)
    
    # Load the fragments
    print("Loading fragments...")
    fragments = load_fragments()

    # 1 Load the problem
    print("Building Assembly Problem...")
    problem = AssemblyProblem(
        fragments=fragments, 
        min_overlap=config.MIN_OVERLAP
    )
    print(f"Problem initialized with {problem.n} fragments.")

    results = []
    
    # 2 Run Consult Oracle    
    print("\nRunning Oracle Ground Truth...")
    oracle_result = oracle_solution(problem=problem, fragments=fragments)
    results.append(oracle_result)
    
    # 3 Run Random Search
    print("\nRunning Baseline: Random Search...")
    rs_result = random_search.optimize(problem, config, rng)
    results.append(rs_result)
    
    # 4 Run Simulated Annealing (hyperparameter sweep)
    sa_results = []
    if getattr(config, "SWEEP_OPT", True):
        print("\nRunning Simulated Annealing Hyperparameter Sweep...")
        sa_initial_temps = _as_list(getattr(config, "OPT_SA_INITIAL_TEMP", getattr(config, "SA_INITIAL_TEMP", 10.0)))
        sa_cooling_rates = _as_list(getattr(config, "OPT_SA_COOLING_RATE", getattr(config, "SA_COOLING_RATE", 0.995)))

        run_id = 0
        for t0 in sa_initial_temps: # optimization
            for alpha in sa_cooling_rates:
                run_id += 1
                print(f"  -> SA run {run_id}: T0={t0}, alpha={alpha}")

                sa_config = SimpleNamespace(
                    SA_INITIAL_TEMP=t0,
                    SA_COOLING_RATE=alpha,
                    SA_MIN_TEMP=getattr(config, "OPT_SA_MIN_TEMP", getattr(config, "SA_MIN_TEMP", 1e-3)),
                    MAX_EVALUATIONS=getattr(config, "MAX_EVALUATIONS", 10_000),
                    MAX_TIME_SEC=getattr(config, "MAX_TIME_SEC", 60.0),
                    RANDOM_SEED=getattr(config, "RANDOM_SEED", 123),
                )
        
                sa_result = simulated_annealing.optimize(problem=problem, config=sa_config, rng=None)
                sa_result["method"] = "Simulated Annealing"
                sa_result["method_label"] = f"SA (T0={t0}, alpha={alpha})"
                sa_result["sweep_run_id"] = run_id
                sa_result["sweep_T0"] = t0
                sa_result["sweep_alpha"] = alpha

                sa_results.append(sa_result)
                results.append(sa_result)
    else:
        print("\nRunning Simulated Annealing (single configuration)...")
        t0 = getattr(config, "SA_INITIAL_TEMP", 10.0)
        alpha = getattr(config, "SA_COOLING_RATE", 0.995)

        sa_config = SimpleNamespace(
            SA_INITIAL_TEMP=t0,
            SA_COOLING_RATE=alpha,
            SA_MIN_TEMP=getattr(config, "OPT_SA_MIN_TEMP", getattr(config, "SA_MIN_TEMP", 1e-3)),
            MAX_EVALUATIONS=getattr(config, "MAX_EVALUATIONS", 10_000),
            MAX_TIME_SEC=getattr(config, "MAX_TIME_SEC", 60.0),
            RANDOM_SEED=getattr(config, "RANDOM_SEED", 123),
        )

        sa_result = simulated_annealing.optimize(problem=problem, config=sa_config, rng=None)
        sa_result["method"] = "Simulated Annealing"
        sa_result["method_label"] = f"SA (single, T0={t0}, alpha={alpha})"
        sa_result["sweep_run_id"] = 1
        sa_result["sweep_T0"] = t0
        sa_result["sweep_alpha"] = alpha

        sa_results.append(sa_result)
        results.append(sa_result)

    #print(f"Completed SA sweep with {len(sa_results)} runs.") # not needed
    
    # 5 Run Genetic Algorithm
    print("\nRunning Genetic Algorithm...")
    ga_result = genetic_algorithm.optimize(problem, config, rng)
    results.append(ga_result)
    
    # Save results to JSON for analyze_results.py script
    os.makedirs("reports", exist_ok=True)
    out_file = "reports/experiment_results.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=4)

    sa_out_file = "reports/simulated_annealing/sa_sweep_results.json"
    os.makedirs("reports/simulated_annealing", exist_ok=True)
    with open(sa_out_file, "w") as f:
        json.dump(sa_results, f, indent=4)

if __name__ == "__main__":
    main()