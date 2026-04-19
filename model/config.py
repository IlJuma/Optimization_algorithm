#
# GLOBAL CONFIGURATION FILE
#

# Global random seed (we currently use 2 different names across the project)
RANDOM_SEED = SEED = 123

# --- generate_chromosome ---
LENGTH = 5000                       # Total sequence length
    # LENGTH = 50_000 # Stage 3 scaling for final experiments, needs finetuning
MEAN_GC_CONTENT = 0.50              # Mean GC across the whole sequence
GC_WINDOW_SIZE = int(LENGTH / 20)   # GC is assigned per window
GC_STDDEV = 0.08                    # Variation in GC between windows
MIN_GC = 0.20                       # Clamp local GC to this minimum
MAX_GC = 0.80                       # Clamp local GC to this maximum
MAX_HOMOPOLYMER = 10                # Max allowed homopolymer length (None to disable)

# --- fragment_chromosome ---
N_GENOME_COPIES = 80                # Number of independent genome molecules to fragment
    # N_GENOME_COPIES = 50 # Stage 3 scaling for final experiments, needs finetuning
MEAN_CUT_SPACING = 350              # Random cut-site model
MIN_INSERT_SIZE = 250               # Size selection after fragmentation
MAX_INSERT_SIZE = 450               # Size selection after fragmentation
RECOVERY_FRACTION = 0.25            # Recovery / sampling after size selection
MAX_FRAGMENTS_TO_KEEP = None        # Optional additional sampling cap
FASTA_LINE_WIDTH = 80
RANDOMIZE_FRAGMENT_ORIENTATION = True   # Randomly flip about 50% of recovered fragments into reverse complements

# --- simulate_illumina ---
PAIRED_END = True
FRAGMENT_SAMPLING_FRACTION = 1.0
MAX_FRAGMENTS_TO_SEQUENCE = None
START_QUALITY = 36
END_QUALITY = 28
QUALITY_JITTER_STDDEV = 2.0
MIN_QUALITY = 2
MAX_QUALITY = 40
READ_NAME_PREFIX = "simread"

# --- Data / Model Parameters ---
# Edge scoring defaults
MIN_OVERLAP = 20  # Neighboring fragments must overlap by MIN_OVERLAP or more to be considered in the same contig
BREAK_FACTOR = 2  # Scales the break penalty: P = BREAK_FACTOR * L
PARTIAL_OVERLAP_FACTOR = 1  # Reduces break penalty when 0 < overlap < MIN_OVERLAP. Use 1 or 2

# Dense mode (n <= DENSE_THRESHOLD):
#   - All pairwise fragment edges (i, j) are precomputed and stored in a full matrix
#   - Fast O(1) lookup during optimization
#   - Higher upfront cost: O(n^2) time and memory

DENSE_THRESHOLD = 200  # Threshold for switching between dense and sparse edge evaluation.

# Sparse mode (n > DENSE_THRESHOLD):
#   - Edge information is computed lazily (on demand) and cached
#   - Only edges actually queried by the algorithm are evaluated
#   - Lower memory usage and avoids unnecessary computations

MAX_NEIGHBORS_PER_FRAGMENT = 50  # Maximum number of best (lowest-cost) feasible neighbors retained per fragment in sparse mode.
# Possible parameter names/values in case of paired reads alignment
    # ALPHA_MISMATCH = 3.0
    # MIN_PAIR_GAP = -50
    # MAX_PAIR_GAP = 150
    # USE_MISMATCH_PENALTY = False
    # LAMBDA_MISMATCH = 1.0
    # BIG_M = 1e9

# --- Shared Evaluation Budget ---
# Every algorithm MUST stop when it hits one of these limits
MAX_EVALUATIONS = 15_000
MAX_TIME_SEC = 60.0

# --- Random Search Specific ---
RS_NUM_SAMPLES = 15_000

# --- Simulated Annealing Specific ---
SA_INITIAL_TEMP = 1.0
SA_COOLING_RATE = 0.999
SA_MIN_TEMP = 1e-3
SA_NEIGHBORHOOD = "swap"

# Optimization Simulated Annealing Specific
SWEEP_OPT = False # this was for the comparison, it will "break" the comparison graph
OPT_SA_INITIAL_TEMP =[1.0, 10.0, 20.0]
OPT_SA_COOLING_RATE = [0.95, 0.99, 0.999]
OPT_SA_MIN_TEMP = 1e-12
OPT_SA_NEIGHBORHOOD = "swap"


# --- Genetic Algorithm Specific ---
GA_POP_SIZE = 30
GA_NUM_GENERATIONS = 500
GA_MUTATION_RATE = 0.5
GA_CROSSOVER_RATE = 0.8
GA_ELITISM = 2


# GA OPTIMIZATION (HYPERPARAMETER SEARCH)
GA_POP_SIZES = [100, 150]
GA_MUTATION_RATES = [0.05, 0.1]
GA_CROSSOVER_RATES = [0.7, 0.9]
GA_ELITISM_VALUES = [1, 3]
