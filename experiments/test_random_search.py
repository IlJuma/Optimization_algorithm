import sys
import os
sys.path.append(os.path.abspath("../algorithms"))
import random
import random_search

# ---------------------------------------------------------
# Mock Config
# ---------------------------------------------------------
class MockConfig:
    MAX_EVALUATIONS = 5000
    MAX_TIME_SEC = 5

config = MockConfig()

# ---------------------------------------------------------
# Mock Problem Class
# ---------------------------------------------------------
class MockAssemblyProblem:
    def __init__(self, num_reads):
        self.num_reads = num_reads
        # We will pretend the "optimal" answer is just [0, 1, 2, ..., n-1]
        self.optimal_solution = list(range(num_reads))
        
    def random_solution(self, rng):
        """Returns a shuffled permutation of read indices."""
        sol = list(range(self.num_reads))
        rng.shuffle(sol)
        return sol
        
    def evaluate(self, solution):
        """
        Mock evaluation: cost is based on how far away each read 
        is from its 'optimal' position. (Lower is better).
        """
        cost = 0
        for i, val in enumerate(solution):
            cost += abs(i - val)
        return float(cost)

# ---------------------------------------------------------
# Run Algorithm
# ---------------------------------------------------------
if __name__ == "__main__":
    rng = random.Random(42)
    # Give it a tiny problem (e.g., 10 reads) so it has a chance to find a good score
    problem = MockAssemblyProblem(num_reads=10)
    
    print("Starting Random Search...")
    result = random_search.optimize(problem, config, rng)
    
    print("\n--- Results ---")
    print(f"Method:        {result['method']}")
    print(f"Evaluations:   {result['evaluations']}")
    print(f"Runtime:       {result['runtime_sec']:.4f} seconds")
    print(f"Best Score:    {result['best_score']}")
    print(f"Best Solution: {result['best_solution']}")