from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Set

from .data_loader_frag import FragmentRecord


# =========================
# CONFIG / DEFAULTS
# =========================

DENSE_THRESHOLD = 200

# Edge scoring defaults
MIN_OVERLAP = 20
PARTIAL_OVERLAP_FACTOR = 1.0  # use 1.0 or 2.0
USE_CANDIDATE_FILTER = True
KMER_SIZE = 8

# Neighbor cap for sparse mode
MAX_NEIGHBORS_PER_FRAGMENT = 50


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

def best_suffix_prefix_overlap(left: str, right: str) -> int:
    """
    Find the longest exact suffix-prefix overlap between two error-free sequences.
    """
    max_len = min(len(left), len(right))

    for ov in range(max_len, 0, -1):
        if left[-ov:] == right[:ov]:
            return ov

    return 0


def kmer_set(seq: str, k: int) -> Set[str]:
    if k <= 0 or len(seq) < k:
        return set()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def quick_candidate_test(
    left_seq: str,
    right_seq: str,
    left_kmers: Set[str],
    right_kmers: Set[str],
) -> bool:
    """
    Cheap plausibility test before computing a full overlap.
    """
    if not left_seq or not right_seq:
        return False

    if left_kmers and right_kmers:
        if left_kmers.isdisjoint(right_kmers):
            return False

    return True


def overlap_edge_info(
    left: str,
    right: str,
    min_overlap: int,
    partial_overlap_factor: float,
) -> EdgeInfo:
    """
    Compute exact-overlap edge information for two error-free fragment sequences.

    Let:
        L = min(len(left), len(right))
        P = 2 * L

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
    P = 2 * L

    overlap = best_suffix_prefix_overlap(left, right)

    if overlap >= min_overlap:
        return EdgeInfo(
            overlap=overlap,
            feasible=True,
            cost=float(L - overlap),
        )

    if overlap > 0:
        return EdgeInfo(
            overlap=overlap,
            feasible=False,
            cost=float(P - partial_overlap_factor * overlap),
        )

    return EdgeInfo(
        overlap=0,
        feasible=False,
        cost=float(P),
    )


# =========================
# MAIN PROBLEM CLASS
# =========================

class AssemblyProblem:
    """
    Fragment ordering problem for overlap-based assembly experiments.

    A solution is represented as:
        List[str]  # ordered fragment_id values

    Design:
        - dense edge matrix for small n
        - lazy sparse caching for larger n
        - exact overlap scoring for fragment-only phase
        - breaks penalized using P = 2 * L
        - partial overlaps below MIN_OVERLAP are still breaks, but receive a
          slightly reduced penalty
    """

    def __init__(
        self,
        fragments: List[FragmentRecord],
        dense_threshold: int = DENSE_THRESHOLD,
        min_overlap: int = MIN_OVERLAP,
        partial_overlap_factor: float = PARTIAL_OVERLAP_FACTOR,
        use_candidate_filter: bool = USE_CANDIDATE_FILTER,
        kmer_size: int = KMER_SIZE,
        max_neighbors_per_fragment: int = MAX_NEIGHBORS_PER_FRAGMENT,
    ):
        self.fragments = fragments
        self.fragment_ids = [f.fragment_id for f in fragments]
        self.n = len(fragments)

        self.dense_threshold = dense_threshold
        self.min_overlap = min_overlap
        self.partial_overlap_factor = partial_overlap_factor
        self.use_candidate_filter = use_candidate_filter
        self.kmer_size = kmer_size
        self.max_neighbors_per_fragment = max_neighbors_per_fragment

        self.use_dense = self.n <= self.dense_threshold

        self.fragment_by_id: Dict[str, FragmentRecord] = {
            f.fragment_id: f for f in fragments
        }
        self.id_to_index: Dict[str, int] = {
            f.fragment_id: i for i, f in enumerate(fragments)
        }
        self.index_to_id: List[str] = self.fragment_ids.copy()

        self.kmer_index: Dict[str, Set[str]] = {
            f.fragment_id: kmer_set(f.sequence, self.kmer_size)
            for f in fragments
        }

        # Dense storage
        self.edge_matrix: Optional[List[List[Optional[EdgeInfo]]]] = None

        # Sparse / lazy storage
        self.edge_cache: Dict[Tuple[str, str], EdgeInfo] = {}
        self.neighbor_cache: Dict[str, List[str]] = {}

        if self.use_dense:
            self._precompute_dense()

    # =========================
    # SOLUTION HELPERS
    # =========================

    def initial_solution(self) -> List[str]:
        return self.fragment_ids.copy()

    def random_solution(self, rng) -> List[str]:
        sol = self.fragment_ids.copy()
        rng.shuffle(sol)
        return sol

    # =========================
    # EDGE COMPUTATION
    # =========================

    def _break_penalty(self, left_id: str, right_id: str) -> float:
        left = self.fragment_by_id[left_id]
        right = self.fragment_by_id[right_id]
        L = min(left.frag_len, right.frag_len)
        return float(2 * L)

    def _self_edge_info(self, fragment_id: str) -> EdgeInfo:
        frag = self.fragment_by_id[fragment_id]
        P = float(2 * frag.frag_len)
        return EdgeInfo(
            overlap=0,
            feasible=False,
            cost=P,
        )

    def _compute_edge(self, left_id: str, right_id: str) -> EdgeInfo:
        if left_id == right_id:
            return self._self_edge_info(left_id)

        left = self.fragment_by_id[left_id]
        right = self.fragment_by_id[right_id]

        if self.use_candidate_filter:
            if not quick_candidate_test(
                left.sequence,
                right.sequence,
                self.kmer_index[left_id],
                self.kmer_index[right_id],
            ):
                return EdgeInfo(
                    overlap=0,
                    feasible=False,
                    cost=self._break_penalty(left_id, right_id),
                )

        return overlap_edge_info(
            left=left.sequence,
            right=right.sequence,
            min_overlap=self.min_overlap,
            partial_overlap_factor=self.partial_overlap_factor,
        )

    def _precompute_dense(self) -> None:
        self.edge_matrix = [
            [None for _ in range(self.n)]
            for _ in range(self.n)
        ]

        for i, left_id in enumerate(self.fragment_ids):
            for j, right_id in enumerate(self.fragment_ids):
                self.edge_matrix[i][j] = self._compute_edge(left_id, right_id)

    # =========================
    # PUBLIC EDGE INTERFACE
    # =========================

    def pair_info(self, left_id: str, right_id: str) -> EdgeInfo:
        if self.use_dense:
            i = self.id_to_index[left_id]
            j = self.id_to_index[right_id]
            return self.edge_matrix[i][j]  # type: ignore[index]

        key = (left_id, right_id)
        if key not in self.edge_cache:
            self.edge_cache[key] = self._compute_edge(left_id, right_id)

        return self.edge_cache[key]

    def pair_cost(self, left_id: str, right_id: str) -> float:
        return self.pair_info(left_id, right_id).cost

    def pair_overlap(self, left_id: str, right_id: str) -> int:
        return self.pair_info(left_id, right_id).overlap

    def is_edge_feasible(self, left_id: str, right_id: str) -> bool:
        return self.pair_info(left_id, right_id).feasible

    def get_neighbors(self, fragment_id: str) -> List[str]:
        """
        Return plausible right-neighbor candidates for one fragment.

        Dense mode:
            returns all feasible neighbors

        Sparse mode:
            lazily computes candidates and caches top neighbors
        """
        if self.use_dense:
            neighbors = []
            for other_id in self.fragment_ids:
                if other_id == fragment_id:
                    continue
                if self.is_edge_feasible(fragment_id, other_id):
                    neighbors.append(other_id)
            return neighbors

        if fragment_id in self.neighbor_cache:
            return self.neighbor_cache[fragment_id]

        candidates = []
        for other_id in self.fragment_ids:
            if other_id == fragment_id:
                continue

            edge = self.pair_info(fragment_id, other_id)
            if edge.feasible:
                candidates.append((other_id, edge.cost))

        candidates.sort(key=lambda x: x[1])
        neighbors = [frag_id for frag_id, _ in candidates[:self.max_neighbors_per_fragment]]
        self.neighbor_cache[fragment_id] = neighbors
        return neighbors

    # =========================
    # SOLUTION EVALUATION
    # =========================

    def evaluate(self, solution: List[str]) -> float:
        """
        Evaluate a full permutation.

        Lower is better.
        """
        if len(solution) <= 1:
            return 0.0

        total = 0.0
        for i in range(len(solution) - 1):
            total += self.pair_cost(solution[i], solution[i + 1])
        return total

    def count_breaks(self, solution: List[str]) -> int:
        breaks = 0
        for i in range(len(solution) - 1):
            if not self.is_edge_feasible(solution[i], solution[i + 1]):
                breaks += 1
        return breaks

    def count_contigs(self, solution: List[str]) -> int:
        if not solution:
            return 0
        return 1 + self.count_breaks(solution)

    def total_overlap(self, solution: List[str]) -> int:
        total = 0
        for i in range(len(solution) - 1):
            total += self.pair_overlap(solution[i], solution[i + 1])
        return total

    def affected_indices_for_swap(self, n: int, a: int, b: int) -> List[int]:
        """
        Return adjacency indices affected by swapping positions a and b.

        Edge i means (solution[i], solution[i+1]).
        """
        affected = set()

        for idx in [a - 1, a, b - 1, b]:
            if 0 <= idx < n - 1:
                affected.add(idx)

        return sorted(affected)

    def delta_swap(self, solution: List[str], a: int, b: int) -> float:
        """
        Compute score change for swapping two positions.

        Returns:
            new_score - old_score
        """
        if a == b:
            return 0.0

        n = len(solution)
        old_indices = self.affected_indices_for_swap(n, a, b)

        old_cost = 0.0
        for idx in old_indices:
            old_cost += self.pair_cost(solution[idx], solution[idx + 1])

        new_solution = solution.copy()
        new_solution[a], new_solution[b] = new_solution[b], new_solution[a]

        new_cost = 0.0
        for idx in old_indices:
            new_cost += self.pair_cost(new_solution[idx], new_solution[idx + 1])

        return new_cost - old_cost

    # =========================
    # ALGORITHM-FACING VIEWS
    # =========================

    def fragments_algorithm_view(self) -> List[Dict[str, object]]:
        """
        Reduced fragment view without truth coordinates.
        """
        return [
            {
                "fragment_id": f.fragment_id,
                "sequence": f.sequence,
                "frag_len": f.frag_len,
            }
            for f in self.fragments
        ]

    def fragments_truth_view(self) -> List[Dict[str, object]]:
        """
        Full fragment view including truth metadata.
        """
        return [
            {
                "fragment_id": f.fragment_id,
                "source": f.source,
                "copy": f.copy,
                "frag_start": f.frag_start,
                "frag_end": f.frag_end,
                "frag_len": f.frag_len,
                "sequence": f.sequence,
            }
            for f in self.fragments
        ]

    # =========================
    # DEBUG / SUMMARY
    # =========================

    def summary(self) -> Dict[str, object]:
        return {
            "n_fragments": self.n,
            "use_dense": self.use_dense,
            "dense_threshold": self.dense_threshold,
            "min_overlap": self.min_overlap,
            "partial_overlap_factor": self.partial_overlap_factor,
            "kmer_size": self.kmer_size,
            "max_neighbors_per_fragment": self.max_neighbors_per_fragment,
            "cached_edges": len(self.edge_cache),
        }


# =========================
# EXAMPLE USAGE
# =========================

if __name__ == "__main__":
    from .data_loader_frag import load_fragments

    fragments = load_fragments()
    problem = AssemblyProblem(fragments)

    print(problem.summary())

    solution = problem.initial_solution()
    print("Initial score:", problem.evaluate(solution))
    print("Initial contigs:", problem.count_contigs(solution))
    print("Initial total overlap:", problem.total_overlap(solution))

    if len(solution) >= 2:
        print("Delta swap(0,1):", problem.delta_swap(solution, 0, 1))