"""
Graph Structure Benchmark Suite for CHILmesh Modernization

Compares four candidate data structure implementations:
1. NetworkX (reference)
2. Custom Compact Graph (recommended)
3. CSR Sparse Matrix + separate element storage
4. Half-Edge (half-mesh) representation

Measures: construction time, adjacency lookup, traversal speed, memory usage
"""

import numpy as np
import time
import sys
from typing import Tuple, List, Dict, Any
from pathlib import Path
from collections import defaultdict

# Try to import dependencies
try:
    import networkx as nx
except ImportError:
    nx = None

try:
    from scipy.sparse import csr_matrix
except ImportError:
    csr_matrix = None


class CompactGraph:
    """
    Recommended: Minimal, NumPy-friendly graph representation.

    Structure:
    - vertices: N×3 coordinate array
    - edge_list: List of (u, v) tuples (unique, undirected)
    - edge2elem: E×2 array [elem_left, elem_right] (-1 for boundary)
    - elem2vert: M×4 array (padded to 4 for mixed triangles/quads)
    - elem2edge: List[List[int]] element→edge indices
    - vert2edge: List[List[int]] vertex→edge indices
    - vert2elem: List[List[int]] vertex→element indices
    """

    def __init__(self, points: np.ndarray, connectivity: np.ndarray):
        self.vertices = points.copy()
        self.n_verts = points.shape[0]
        self.n_elems = connectivity.shape[0]
        self.elem2vert = connectivity.copy()

        # Initialize edge structure
        self._build_edges()
        self._build_adjacencies()

    def _build_edges(self):
        """O(n log n) edge discovery via sorting."""
        edges = set()
        for elem_id, elem in enumerate(self.elem2vert):
            n_verts = 3 if elem[2] == elem[3] else 4
            for i in range(n_verts):
                v1 = elem[i]
                v2 = elem[(i + 1) % n_verts]
                edge = tuple(sorted([v1, v2]))
                edges.add(edge)

        self.edge_list = sorted(list(edges))
        self.edge_id = {edge: idx for idx, edge in enumerate(self.edge_list)}
        self.n_edges = len(self.edge_list)

    def _build_adjacencies(self):
        """Build elem2edge, vert2edge, vert2elem maps."""
        self.elem2edge = [[] for _ in range(self.n_elems)]
        self.vert2edge = [set() for _ in range(self.n_verts)]
        self.vert2elem = [set() for _ in range(self.n_verts)]
        self.edge2elem = np.full((self.n_edges, 2), -1, dtype=int)

        for elem_id, elem in enumerate(self.elem2vert):
            n_verts = 3 if elem[2] == elem[3] else 4
            for i in range(n_verts):
                v1 = elem[i]
                v2 = elem[(i + 1) % n_verts]
                edge = tuple(sorted([v1, v2]))
                edge_id = self.edge_id[edge]

                self.elem2edge[elem_id].append(edge_id)
                self.vert2edge[v1].add(edge_id)
                self.vert2elem[v1].add(elem_id)

                # Track both elements for this edge
                if self.edge2elem[edge_id, 0] == -1:
                    self.edge2elem[edge_id, 0] = elem_id
                else:
                    self.edge2elem[edge_id, 1] = elem_id

        # Convert sets to sorted lists for deterministic iteration
        self.vert2edge = [sorted(list(s)) for s in self.vert2edge]
        self.vert2elem = [sorted(list(s)) for s in self.vert2elem]

    def neighbors(self, vert_id: int) -> List[int]:
        """O(degree) neighbor lookup."""
        neighbors = set()
        for edge_id in self.vert2edge[vert_id]:
            u, v = self.edge_list[edge_id]
            neighbors.add(v if u == vert_id else u)
        return sorted(list(neighbors))


class NetworkXGraph:
    """Reference implementation using NetworkX."""

    def __init__(self, points: np.ndarray, connectivity: np.ndarray):
        if nx is None:
            raise ImportError("NetworkX not installed")

        self.vertices = points.copy()
        self.n_verts = points.shape[0]
        self.n_elems = connectivity.shape[0]
        self.elem2vert = connectivity.copy()

        self.G = nx.Graph()
        self.G.add_nodes_from(range(self.n_verts))

        # Add edges
        edges = set()
        for elem in connectivity:
            n_verts = 3 if elem[2] == elem[3] else 4
            for i in range(n_verts):
                v1, v2 = elem[i], elem[(i + 1) % n_verts]
                edge = tuple(sorted([v1, v2]))
                edges.add(edge)

        self.G.add_edges_from(edges)
        self.n_edges = len(edges)

    def neighbors(self, vert_id: int) -> List[int]:
        """O(degree) neighbor lookup."""
        return list(self.G.neighbors(vert_id))


class CSRGraph:
    """Sparse matrix adjacency representation."""

    def __init__(self, points: np.ndarray, connectivity: np.ndarray):
        if csr_matrix is None:
            raise ImportError("scipy not installed")

        self.vertices = points.copy()
        self.n_verts = points.shape[0]
        self.n_elems = connectivity.shape[0]
        self.elem2vert = connectivity.copy()

        # Build adjacency matrix
        edges = set()
        for elem in connectivity:
            n_verts = 3 if elem[2] == elem[3] else 4
            for i in range(n_verts):
                v1, v2 = elem[i], elem[(i + 1) % n_verts]
                edges.add(tuple(sorted([v1, v2])))

        row, col = zip(*edges)
        data = np.ones(len(edges))
        self.adj_matrix = csr_matrix(
            (data, (row + col, col + row)),
            shape=(self.n_verts, self.n_verts)
        )
        self.n_edges = len(edges)

    def neighbors(self, vert_id: int) -> List[int]:
        """O(k) neighbor lookup where k = degree."""
        return sorted(self.adj_matrix[vert_id].nonzero()[1].tolist())


def create_test_meshes() -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Create test meshes of varying sizes."""

    # Try to load from CHILmesh examples
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
        import chilmesh

        return {
            "annulus": (chilmesh.examples.annulus().points,
                       chilmesh.examples.annulus().connectivity_list),
            "structured": (chilmesh.examples.structured().points,
                          chilmesh.examples.structured().connectivity_list),
            "donut": (chilmesh.examples.donut().points,
                     chilmesh.examples.donut().connectivity_list),
        }
    except Exception as e:
        print(f"Warning: Could not load CHILmesh examples: {e}")
        return {}


def benchmark_construction(name: str, GraphClass, points: np.ndarray,
                          connectivity: np.ndarray) -> float:
    """Measure construction time."""
    start = time.perf_counter()
    graph = GraphClass(points, connectivity)
    elapsed = time.perf_counter() - start
    return elapsed


def benchmark_neighbor_lookup(graph, n_samples: int = 100) -> float:
    """Measure average neighbor lookup time."""
    sample_verts = np.random.choice(graph.n_verts, size=min(n_samples, graph.n_verts),
                                    replace=False)
    start = time.perf_counter()
    for v in sample_verts:
        _ = graph.neighbors(v)
    elapsed = time.perf_counter() - start
    return elapsed / len(sample_verts)


def benchmark_bfs(graph, start_vert: int = 0) -> float:
    """Measure BFS traversal time."""
    visited = set()
    queue = [start_vert]

    start = time.perf_counter()
    while queue:
        v = queue.pop(0)
        if v in visited:
            continue
        visited.add(v)
        for u in graph.neighbors(v):
            if u not in visited:
                queue.append(u)
    elapsed = time.perf_counter() - start
    return elapsed


def run_benchmarks():
    """Run all benchmarks."""

    print("=" * 80)
    print("CHILmesh Graph Structure Benchmark Suite")
    print("=" * 80)

    meshes = create_test_meshes()
    if not meshes:
        print("Error: Could not load test meshes. Exiting.")
        return

    implementations = [
        ("CompactGraph", CompactGraph),
        ("NetworkX", NetworkXGraph if nx else None),
        ("CSR", CSRGraph if csr_matrix else None),
    ]

    results = defaultdict(dict)

    for mesh_name, (points, connectivity) in meshes.items():
        print(f"\n{'Mesh: ' + mesh_name:-^80}")
        print(f"  Vertices: {points.shape[0]}, Elements: {connectivity.shape[0]}")

        for impl_name, GraphClass in implementations:
            if GraphClass is None:
                print(f"  {impl_name:20s}: SKIPPED (dependency not installed)")
                continue

            try:
                # Construction
                construct_time = benchmark_construction(
                    impl_name, GraphClass, points, connectivity
                )

                # Create graph instance for further benchmarks
                graph = GraphClass(points, connectivity)

                # Neighbor lookup
                lookup_time = benchmark_neighbor_lookup(graph)

                # BFS
                bfs_time = benchmark_bfs(graph)

                print(f"  {impl_name:20s}:")
                print(f"    Construction:  {construct_time*1000:10.3f} ms")
                print(f"    Neighbor (avg): {lookup_time*1e6:10.3f} μs")
                print(f"    BFS:            {bfs_time*1000:10.3f} ms")

                results[mesh_name][impl_name] = {
                    "construct_ms": construct_time * 1000,
                    "lookup_us": lookup_time * 1e6,
                    "bfs_ms": bfs_time * 1000,
                    "edges": graph.n_edges,
                }

            except Exception as e:
                print(f"  {impl_name:20s}: ERROR - {e}")

    # Summary table
    print("\n" + "=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)

    for mesh_name, impl_results in results.items():
        print(f"\n{mesh_name}:")
        print(f"  {'Implementation':<20} {'Construct (ms)':<15} {'Lookup (μs)':<15} {'BFS (ms)':<15}")
        for impl, metrics in impl_results.items():
            print(f"  {impl:<20} {metrics['construct_ms']:<15.3f} "
                  f"{metrics['lookup_us']:<15.3f} {metrics['bfs_ms']:<15.3f}")


if __name__ == "__main__":
    run_benchmarks()
