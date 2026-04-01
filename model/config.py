#
# GLOBAL CONFIGURATION FILE
#

# --- Data / Model Parameters ---
MIN_OVERLAP = 20
ALPHA_MISMATCH = 3.0
MIN_PAIR_GAP = -50
MAX_PAIR_GAP = 150
USE_MISMATCH_PENALTY = False
LAMBDA_MISMATCH = 1.0
BIG_M = 1e9

# P = 2 * L penalty for breaks
PARTIAL_OVERLAP_FACTOR = 1.0 

# --- Instance / Simulation Size ---
# Using Stage 3 scaling for final experiments
GENOME_LENGTH = 50_000 
N_GENOME_COPIES = 50
RANDOM_SEED = 123

# --- Shared Evaluation Budget ---
# Every algorithm MUST stop when it hits one of these limits
MAX_EVALUATIONS = 10_000
MAX_TIME_SEC = 60.0

# --- Random Search Specific ---
RS_NUM_SAMPLES = 10_000

# --- Simulated Annealing Specific ---
SA_INITIAL_TEMP = 10.0
SA_COOLING_RATE = 0.995
SA_MIN_TEMP = 1e-3
SA_NEIGHBORHOOD = "swap"

# --- Genetic Algorithm Specific ---
GA_POP_SIZE = 100
GA_NUM_GENERATIONS = 200
GA_MUTATION_RATE = 0.1
GA_CROSSOVER_RATE = 0.8
GA_ELITISM = 2