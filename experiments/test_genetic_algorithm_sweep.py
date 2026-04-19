import os
import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from algorithms import Genetic_Algorithm
from experiments import run_experiments


class MockAssemblyProblem:
    def __init__(self, num_fragments):
        self.base_fragment_ids = [f"frag_{i}" for i in range(num_fragments)]

    def evaluate(self, solution):
        cost = 0.0
        for index, fragment in enumerate(solution):
            base_id = fragment.rsplit("_", 1)[0]
            expected_index = int(base_id.split("_")[1])
            cost += abs(index - expected_index)
        return cost

    def count_breaks(self, solution):
        return 0

    def count_contigs(self, solution):
        return 1 if solution else 0

    def total_overlap(self, solution):
        return len(solution)


class GeneticAlgorithmSweepTest(unittest.TestCase):
    def test_optimize_suppresses_generation_logs_when_ga_verbose_is_false(self):
        problem = MockAssemblyProblem(num_fragments=6)
        config = SimpleNamespace(
            GA_POP_SIZE=8,
            GA_NUM_GENERATIONS=3,
            GA_MUTATION_RATE=0.2,
            GA_CROSSOVER_RATE=0.7,
            GA_ELITISM=1,
            MAX_EVALUATIONS=100,
            MAX_TIME_SEC=5.0,
            GA_VERBOSE=False,
        )

        with patch("builtins.print") as mock_print:
            Genetic_Algorithm.optimize(
                problem=problem,
                config=config,
                rng=np.random.default_rng(123),
            )

        generation_logs = [
            call.args[0]
            for call in mock_print.call_args_list
            if call.args and isinstance(call.args[0], str) and call.args[0].startswith("Generation ")
        ]
        self.assertEqual(generation_logs, [])

    def test_run_genetic_algorithm_suite_keeps_only_best_run_for_main_results(self):
        config = SimpleNamespace(
            RANDOM_SEED=123,
            MAX_EVALUATIONS=100,
            MAX_TIME_SEC=5.0,
            GA_ELITISM=3,
            GA_VERBOSE=False,
            GA_SWEEP_OPT=True,
            GA_POP_SIZE=150,
            GA_NUM_GENERATIONS=300,
            GA_MUTATION_RATE=0.05,
            GA_CROSSOVER_RATE=0.8,
            OPT_GA_POP_SIZE=[100, 200],
            OPT_GA_NUM_GENERATIONS=[20],
            OPT_GA_MUTATION_RATE=[0.02],
            OPT_GA_CROSSOVER_RATE=[0.7, 0.9],
        )

        score_by_combo = {
            (100, 20, 0.02, 0.7): 90.0,
            (100, 20, 0.02, 0.9): 75.0,
            (200, 20, 0.02, 0.7): 60.0,
            (200, 20, 0.02, 0.9): 45.0,
        }

        def fake_optimize(problem=None, config=None, rng=None):
            combo = (
                config.GA_POP_SIZE,
                config.GA_NUM_GENERATIONS,
                config.GA_MUTATION_RATE,
                config.GA_CROSSOVER_RATE,
            )
            return {
                "method": "Genetic Algorithm",
                "best_solution": ["frag_0_F"],
                "best_score": score_by_combo[combo],
                "history": [score_by_combo[combo]],
                "current_cost_history": [score_by_combo[combo]],
                "contigs_history": [1],
                "overlap_history": [1],
                "evaluations": 1,
                "runtime_sec": 0.1,
                "breaks": 0,
                "contigs": 1,
                "total_overlap": 1,
            }

        with patch.object(run_experiments.Genetic_Algorithm, "optimize", side_effect=fake_optimize):
            best_run, sweep_runs = run_experiments.run_genetic_algorithm_suite(
                problem=object(),
                base_config=config,
            )

        self.assertEqual(len(sweep_runs), 4)
        self.assertEqual(best_run["best_score"], 45.0)
        self.assertEqual(best_run["sweep_run_id"], 4)
        self.assertEqual(best_run["sweep_population_size"], 200)
        self.assertEqual(best_run["sweep_num_generations"], 20)
        self.assertEqual(best_run["sweep_mutation_rate"], 0.02)
        self.assertEqual(best_run["sweep_crossover_rate"], 0.9)
        self.assertIn("GA (pop=200, gen=20, mut=0.02, cross=0.9)", best_run["method_label"])


if __name__ == "__main__":
    unittest.main()
