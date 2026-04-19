import json
import random
import os
import sys
from types import SimpleNamespace
from itertools import product

# Import the shared config and problem model
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from model import config
from model.data_loader_frag import load_fragments
from model.problem import AssemblyProblem


# Import the algorithms
from algorithms.oracle_solution import oracle_solution
from algorithms import random_search
from algorithms import simulated_annealing # Uncomment when ready
from algorithms import Genetic_Algorithm


def _as_list(value):
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def _best_result(results):
    if not results:
        return None
    return min(results, key=lambda result: result.get("best_score", float("inf")))


def _build_ga_run_configs(base_config):
    if getattr(base_config, "GA_SWEEP_OPT", False):
        pop_sizes = _as_list(getattr(base_config, "OPT_GA_POP_SIZE", getattr(base_config, "GA_POP_SIZE", 150)))
        generation_counts = _as_list(getattr(base_config, "OPT_GA_NUM_GENERATIONS", getattr(base_config, "GA_NUM_GENERATIONS", 300)))
        mutation_rates = _as_list(getattr(base_config, "OPT_GA_MUTATION_RATE", getattr(base_config, "GA_MUTATION_RATE", 0.05)))
        crossover_rates = _as_list(getattr(base_config, "OPT_GA_CROSSOVER_RATE", getattr(base_config, "GA_CROSSOVER_RATE", 0.8)))
    else:
        pop_sizes = [getattr(base_config, "GA_POP_SIZE", 150)]
        generation_counts = [getattr(base_config, "GA_NUM_GENERATIONS", 300)]
        mutation_rates = [getattr(base_config, "GA_MUTATION_RATE", 0.05)]
        crossover_rates = [getattr(base_config, "GA_CROSSOVER_RATE", 0.8)]

    run_configs = []
    run_id = 0
    for pop_size, num_generations, mutation_rate, crossover_rate in product(
        pop_sizes,
        generation_counts,
        mutation_rates,
        crossover_rates,
    ):
        run_id += 1
        run_configs.append(
            (
                run_id,
                SimpleNamespace(
                    GA_POP_SIZE=pop_size,
                    GA_NUM_GENERATIONS=num_generations,
                    GA_MUTATION_RATE=mutation_rate,
                    GA_CROSSOVER_RATE=crossover_rate,
                    GA_ELITISM=getattr(base_config, "GA_ELITISM", 3),
                    GA_VERBOSE=getattr(base_config, "GA_VERBOSE", False),
                    MAX_EVALUATIONS=getattr(base_config, "MAX_EVALUATIONS", 10_000),
                    MAX_TIME_SEC=getattr(base_config, "MAX_TIME_SEC", 60.0),
                    RANDOM_SEED=getattr(base_config, "RANDOM_SEED", 123),
                ),
            )
        )
    return run_configs


def run_genetic_algorithm_suite(problem, base_config):
    ga_runs = []
    run_configs = _build_ga_run_configs(base_config)
    is_sweep = getattr(base_config, "GA_SWEEP_OPT", False)

    if is_sweep:
        print("\nRunning Genetic Algorithm Hyperparameter Sweep...")
    else:
        print("\nRunning Genetic Algorithm (single configuration)...")

    for run_id, ga_config in run_configs:
        pop_size = ga_config.GA_POP_SIZE
        num_generations = ga_config.GA_NUM_GENERATIONS
        mutation_rate = ga_config.GA_MUTATION_RATE
        crossover_rate = ga_config.GA_CROSSOVER_RATE

        print(
            f"  -> GA run {run_id}: "
            f"pop={pop_size}, gen={num_generations}, mut={mutation_rate}, cross={crossover_rate}"
        )

        ga_result = Genetic_Algorithm.optimize(problem=problem, config=ga_config, rng=None)
        ga_result["method"] = "Genetic Algorithm"
        ga_result["method_label"] = (
            f"GA (pop={pop_size}, gen={num_generations}, "
            f"mut={mutation_rate}, cross={crossover_rate})"
        )
        ga_result["sweep_run_id"] = run_id
        ga_result["sweep_population_size"] = pop_size
        ga_result["sweep_num_generations"] = num_generations
        ga_result["sweep_mutation_rate"] = mutation_rate
        ga_result["sweep_crossover_rate"] = crossover_rate

        ga_runs.append(ga_result)

    return _best_result(ga_runs), ga_runs

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

    # 5 Run Genetic Algorithm (sweep stored separately, best run kept for comparison)
    best_ga_result, ga_results = run_genetic_algorithm_suite(problem=problem, base_config=config)
    if best_ga_result is not None:
        results.append(best_ga_result)
    
    # Save results to JSON for analyze_results.py script
    os.makedirs("reports", exist_ok=True)
    out_file = "reports/experiment_results.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=4)

    sa_out_file = "reports/simulated_annealing/sa_sweep_results.json"
    os.makedirs("reports/simulated_annealing", exist_ok=True)
    with open(sa_out_file, "w") as f:
        json.dump(sa_results, f, indent=4)

    ga_out_file = "reports/genetic_algorithm/ga_sweep_results.json"
    os.makedirs("reports/genetic_algorithm", exist_ok=True)
    with open(ga_out_file, "w") as f:
        json.dump(ga_results, f, indent=4)

if __name__ == "__main__":
    main()
