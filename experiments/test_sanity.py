import random
import time
import sys
import os
from test_random_search import MockAssemblyProblem # reuse your mock from earlier
sys.path.append(os.path.abspath("../algorithms"))
from random_search import optimize


class StrictConfig:
    MAX_EVALUATIONS = 100
    MAX_TIME_SEC = 0.001 # 1 millisecond

def test_evaluation_limit():
    config = StrictConfig()
    config.MAX_TIME_SEC = 100 # Give it plenty of time
    problem = MockAssemblyProblem(50)
    
    result = optimize(problem, config, random.Random(1))
    assert result["evaluations"] == 100, f"Failed: Did {result['evaluations']} evals instead of 100"
    print("Evaluation budget strictly respected")

def test_time_limit():
    config = StrictConfig()
    config.MAX_EVALUATIONS = 1_000_000 # Give it infinite evaluations
    config.MAX_TIME_SEC = 0.1 # 100 ms limit
    problem = MockAssemblyProblem(50)
    
    # Make the mock evaluation artificially slow to test the timeout
    original_eval = problem.evaluate
    problem.evaluate = lambda sol: time.sleep(0.001) or original_eval(sol)
    
    result = optimize(problem, config, random.Random(2))
    assert result["runtime_sec"] <= 0.15, f"Failed: Took too long ({result['runtime_sec']}s)"
    assert result["evaluations"] < config.MAX_EVALUATIONS, "Failed: Hit eval limit, not time limit"
    print(f"Time budget strictly respected. Stopped after {result['evaluations']} evals.")

if __name__ == "__main__":
    test_evaluation_limit()
    test_time_limit()