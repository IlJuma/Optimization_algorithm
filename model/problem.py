from __future__ import annotations

from dataclasses import dataclass

from .data_loader_frag import (
    FragmentRecord,
    split_oriented_fragment_id,
    oriented_fragment_sequence,
)


# =========================
# CONFIG / DEFAULTS
# =========================

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


# =========================
# DATA CLASSES
# =========================

@dataclass
class EdgeInfo:
    overlap: int
    feasible: bool
    cost: float


# =========================
# HELPER FUNCTIONS
# =========================

def best_suffix_prefix_overlap(left, right):
    """
    Find the longest exact suffix-prefix overlap between two error-free sequences.
    """
    max_len = min(len(left), len(right))

    for ov in range(max_len, 0, -1):
        if left[-ov:] == right[:ov]:
            return ov

    return 0


def overlap_edge_info(left, right, min_overlap, break_factor, partial_overlap_factor):
    """
    Compute exact-overlap edge information for two error-free fragment sequences.

    Let:
        L = min(len(left), len(right))
        P = break_factor * L

    Cost:
        if overlap >= min_overlap:
            cost = L - overlap
        elif overlap > 0:
            cost = P - partial_overlap_factor * overlap
        else:
            cost = P

    Interpretation:
        - feasible=True means the two fragments can belong to the same contig
        - feasible=False means this adjacency is treated as a break
        - partial sub-threshold overlap is still a break, but better than zero overlap
    """
    L = min(len(left), len(right))
    P = break_factor * L

    overlap = best_suffix_prefix_overlap(left, right)

    if overlap >= min_overlap:
        return EdgeInfo(
            overlap=overlap,
            feasible=True,
            cost=L - overlap,
        )

    if overlap > 0:
        return EdgeInfo(
            overlap=overlap,
            feasible=False,
            cost=P - partial_overlap_factor * overlap,
        )

    return EdgeInfo(
        overlap=0,
        feasible=False,
        cost=P,
    )


# =========================
# MAIN PROBLEM CLASS
# =========================

class AssemblyProblem:
    """
    AssemblyProblem defines the scoring model for fragment ordering.

    Orientation convention:
    - A solution is a list of oriented fragment ids:
        <fragment_id>_F  -> stored fragment sequence as loaded from fragments.fasta
        <fragment_id>_R  -> reverse complement of the stored fragment sequence

    Key ideas:
    - Each base fragment must appear exactly once in a solution, with either _F or _R
    - Adjacent oriented fragments are scored using suffix-prefix overlap
    - Valid overlaps (>= MIN_OVERLAP) reduce cost
    - Invalid or weak overlaps are treated as contig breaks and penalized
      with P = BREAK_FACTOR * min(fragment lengths)
    - Partial overlaps below threshold receive a reduced penalty
    """

    def __init__(
        self,
        fragments,
        dense_threshold=DENSE_THRESHOLD,
        min_overlap=MIN_OVERLAP,
        break_factor=BREAK_FACTOR,
        partial_overlap_factor=PARTIAL_OVERLAP_FACTOR,
        max_neighbors_per_fragment=MAX_NEIGHBORS_PER_FRAGMENT,
    ):
        self.fragments = fragments

        self.base_fragment_ids = [f.fragment_id for f in fragments]
        self.fragment_ids = self.base_fragment_ids.copy()  # compatibility with older code
        self.oriented_fragment_ids = [
            f"{fragment_id}_{orientation}"
            for fragment_id in self.base_fragment_ids
            for orientation in ["F", "R"]
        ]

        self.n_fragments = len(fragments)
        self.n_oriented_fragments = len(self.oriented_fragment_ids)
        self.n = self.n_fragments  # compatibility with older code

        self.dense_threshold = dense_threshold
        self.min_overlap = min_overlap
        self.break_factor = break_factor
        self.partial_overlap_factor = partial_overlap_factor
        self.max_neighbors_per_fragment = max_neighbors_per_fragment

        self.use_dense = self.n_fragments <= self.dense_threshold

        self.fragment_by_id = {
            f.fragment_id: f for f in fragments
        }
        self.id_to_index = {
            oriented_id: i for i, oriented_id in enumerate(self.oriented_fragment_ids)
        }
        self.index_to_id = self.oriented_fragment_ids.copy()

        # Dense storage
        self.edge_matrix = None

        # Sparse / lazy storage
        self.edge_cache = {}
        self.neighbor_cache = {}

        if self.use_dense:
            self._precompute_dense()

    # =========================
    # ORIENTATION HELPERS
    # =========================

    def base_fragment_id(self, oriented_fragment_id):
        fragment_id, _ = split_oriented_fragment_id(oriented_fragment_id)
        return fragment_id

    def oriented_sequence(self, oriented_fragment_id):
        fragment_id, orientation = split_oriented_fragment_id(oriented_fragment_id)
        fragment = self.fragment_by_id[fragment_id]
        return oriented_fragment_sequence(fragment, orientation)

    def validate_solution(self, solution):
        """
        Check that a solution uses each base fragment exactly once.
        """
        if len(solution) != self.n_fragments:
            return False

        seen = set()
        for oriented_id in solution:
            base_id = self.base_fragment_id(oriented_id)
            if base_id in seen:
                return False
            seen.add(base_id)

        return len(seen) == self.n_fragments

    # =========================
    # EDGE COMPUTATION
    # =========================

    def _break_penalty(self, left_id, right_id):
        left_base = self.base_fragment_id(left_id)
        right_base = self.base_fragment_id(right_id)

        left = self.fragment_by_id[left_base]
        right = self.fragment_by_id[right_base]
        L = min(left.frag_len, right.frag_len)

        return self.break_factor * L

    def _compute_edge(self, left_id, right_id):
        left_base = self.base_fragment_id(left_id)
        right_base = self.base_fragment_id(right_id)

        if left_base == right_base:
            return EdgeInfo(
                overlap=0,
                feasible=False,
                cost=self._break_penalty(left_id, right_id),
            )

        left_seq = self.oriented_sequence(left_id)
        right_seq = self.oriented_sequence(right_id)

        return overlap_edge_info(
            left=left_seq,
            right=right_seq,
            min_overlap=self.min_overlap,
            break_factor=self.break_factor,
            partial_overlap_factor=self.partial_overlap_factor,
        )

    def _precompute_dense(self):
        self.edge_matrix = [
            [None for _ in range(self.n_oriented_fragments)]
            for _ in range(self.n_oriented_fragments)
        ]

        for i, left_id in enumerate(self.oriented_fragment_ids):
            for j, right_id in enumerate(self.oriented_fragment_ids):
                self.edge_matrix[i][j] = self._compute_edge(left_id, right_id)

    # =========================
    # PUBLIC EDGE INTERFACE
    # =========================

    def pair_info(self, left_id, right_id):
        if self.use_dense:
            i = self.id_to_index[left_id]
            j = self.id_to_index[right_id]
            return self.edge_matrix[i][j]

        key = (left_id, right_id)
        if key not in self.edge_cache:
            self.edge_cache[key] = self._compute_edge(left_id, right_id)

        return self.edge_cache[key]

    def pair_cost(self, left_id, right_id):
        return self.pair_info(left_id, right_id).cost

    def pair_overlap(self, left_id, right_id):
        return self.pair_info(left_id, right_id).overlap

    def is_edge_feasible(self, left_id, right_id):
        return self.pair_info(left_id, right_id).feasible

    def get_neighbors(self, oriented_fragment_id):
        """
        Return oriented fragments that can follow this oriented fragment in the
        same contig (i.e. overlap >= MIN_OVERLAP), excluding both orientations
        of the same base fragment.
        """
        if self.use_dense:
            neighbors = []
            left_base = self.base_fragment_id(oriented_fragment_id)

            for other_id in self.oriented_fragment_ids:
                if self.base_fragment_id(other_id) == left_base:
                    continue
                if self.is_edge_feasible(oriented_fragment_id, other_id):
                    neighbors.append(other_id)
            return neighbors

        if oriented_fragment_id in self.neighbor_cache:
            return self.neighbor_cache[oriented_fragment_id]

        candidates = []
        left_base = self.base_fragment_id(oriented_fragment_id)

        for other_id in self.oriented_fragment_ids:
            if self.base_fragment_id(other_id) == left_base:
                continue

            edge = self.pair_info(oriented_fragment_id, other_id)
            if edge.feasible:
                candidates.append((other_id, edge.cost))

        candidates.sort(key=lambda x: x[1])
        neighbors = [fragment_id for fragment_id, _ in candidates[:self.max_neighbors_per_fragment]]
        self.neighbor_cache[oriented_fragment_id] = neighbors
        return neighbors

    # =========================
    # SOLUTION EVALUATION
    # =========================

    def evaluate(self, solution):
        """
        Evaluate a full permutation of oriented fragment ids.

        Lower is better.
        """
        if len(solution) <= 1:
            return 0

        total = 0
        for i in range(len(solution) - 1):
            total += self.pair_cost(solution[i], solution[i + 1])
        return total

    def count_breaks(self, solution):
        breaks = 0
        for i in range(len(solution) - 1):
            if not self.is_edge_feasible(solution[i], solution[i + 1]):
                breaks += 1
        return breaks

    def count_contigs(self, solution):
        if not solution:
            return 0
        return 1 + self.count_breaks(solution)

    def total_overlap(self, solution):
        total = 0
        for i in range(len(solution) - 1):
            total += self.pair_overlap(solution[i], solution[i + 1])
        return total