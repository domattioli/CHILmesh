#include "halfedge.hpp"
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cassert>

namespace chilmesh {

// Pack two 32-bit ints into one 64-bit key (canonical: lo < hi)
static inline uint64_t edge_key(int32_t a, int32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

void HalfEdgeMesh::build(const double* pts, int n_pts,
                          const int32_t* conn, int n_el, int vpe) {
    n_verts = n_pts;
    n_elems = n_el;
    max_vpe = vpe;

    px.assign(pts, pts + n_pts);
    py.assign(pts + n_pts, pts + 2 * n_pts);  // caller packs x0..xN y0..yN
    // Actually caller passes interleaved [n,2] — adjust in bindings

    connectivity.assign(conn, conn + n_el * vpe);

    // Allocate half-edges: each elem edge = 1 half-edge
    // Each element with k real vertices has k half-edges
    // triangles in quad mesh have 3 real edges, quads have 4

    // Count real verts per elem
    std::vector<int32_t> vpe_arr(n_el);
    for (int e = 0; e < n_el; ++e) {
        const int32_t* row = conn + e * vpe;
        int k = vpe;
        // -1 padding means triangle
        while (k > 3 && row[k-1] == -1) --k;
        vpe_arr[e] = k;
    }

    // Total half-edges
    int n_he = 0;
    std::vector<int32_t> face_he_start(n_el);
    for (int e = 0; e < n_el; ++e) {
        face_he_start[e] = n_he;
        n_he += vpe_arr[e];
    }
    half_edges.resize(n_he);
    face_he = face_he_start;

    // Fill in vert, face, next, prev for each face loop
    for (int e = 0; e < n_el; ++e) {
        const int32_t* row = conn + e * vpe;
        int k = vpe_arr[e];
        int base = face_he_start[e];
        for (int i = 0; i < k; ++i) {
            HalfEdge& he = half_edges[base + i];
            he.vert = row[i];
            he.face = e;
            he.next = base + (i + 1) % k;
            he.prev = base + (i + k - 1) % k;
            he.twin = -1;  // to be filled
        }
    }

    // Build twin map: directed edge (a→b) matches twin (b→a)
    // Key: edge_key(a,b) with canonical order; store both directed HE ids
    // Use map from directed edge to HE id
    std::unordered_map<uint64_t, int32_t> dir_map;
    dir_map.reserve(n_he * 2);
    // directed key: a*N + b  (but N can be large, use 64-bit pack differently)
    auto dir_key = [&](int32_t a, int32_t b) -> uint64_t {
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };

    for (int h = 0; h < n_he; ++h) {
        int32_t a = half_edges[h].vert;
        int32_t b = half_edges[half_edges[h].next].vert;
        dir_map[dir_key(a, b)] = h;
    }

    // Match twins
    for (int h = 0; h < n_he; ++h) {
        if (half_edges[h].twin != -1) continue;
        int32_t a = half_edges[h].vert;
        int32_t b = half_edges[half_edges[h].next].vert;
        auto it = dir_map.find(dir_key(b, a));
        if (it != dir_map.end()) {
            int32_t t = it->second;
            half_edges[h].twin = t;
            half_edges[t].twin = h;
        }
    }

    // Build vert_he: first HE with vert[h] == v
    vert_he.assign(n_verts, -1);
    for (int h = 0; h < n_he; ++h) {
        int32_t v = half_edges[h].vert;
        if (vert_he[v] == -1) vert_he[v] = h;
    }
}

void HalfEdgeMesh::build_adjacency() {
    int n_he = (int)half_edges.size();

    // Assign edge IDs: canonical undirected edges
    // edge = min(h, twin(h)) for interior; boundary has twin=-1
    std::vector<int32_t> he_to_edge(n_he, -1);
    int edge_count = 0;
    for (int h = 0; h < n_he; ++h) {
        int32_t t = half_edges[h].twin;
        if (t == -1 || h < t) {
            he_to_edge[h] = edge_count++;
            if (t != -1) he_to_edge[t] = he_to_edge[h];
        }
    }
    n_edges = edge_count;

    // edge2vert: [n_edges, 2]
    edge2vert.assign(n_edges * 2, -1);
    for (int h = 0; h < n_he; ++h) {
        int32_t t = half_edges[h].twin;
        if (t == -1 || h < t) {
            int32_t eid = he_to_edge[h];
            int32_t a = half_edges[h].vert;
            int32_t b = half_edges[half_edges[h].next].vert;
            edge2vert[eid * 2] = std::min(a, b);
            edge2vert[eid * 2 + 1] = std::max(a, b);
        }
    }

    // elem2edge: [n_elems, max_vpe]
    elem2edge.assign(static_cast<size_t>(n_elems) * max_vpe, -1);
    for (int e = 0; e < n_elems; ++e) {
        int base_he = face_he[e];
        int k = max_vpe;
        while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
        for (int i = 0; i < k; ++i) {
            int h = base_he + i;
            elem2edge[e * max_vpe + i] = he_to_edge[h];
        }
    }

    // vert2edge: jagged
    vert2edge.assign(n_verts, {});
    for (int eid = 0; eid < n_edges; ++eid) {
        int32_t v0 = edge2vert[eid * 2];
        int32_t v1 = edge2vert[eid * 2 + 1];
        vert2edge[v0].push_back(eid);
        vert2edge[v1].push_back(eid);
    }

    adjacency_built = true;
}

void HalfEdgeMesh::skeletonize() {
    if (!adjacency_built) build_adjacency();

    // elem2edge is [n_elems * max_vpe]
    // edge2elem: build from half-edges
    std::vector<int32_t> edge2elem(n_edges * 2, -1);
    for (int h = 0; h < (int)half_edges.size(); ++h) {
        int32_t t = half_edges[h].twin;
        int32_t eid = -1;
        if (t == -1 || h < t) {
            eid = -1; // canonical
        }
        // simpler: directly from elem2edge
    }
    // Build edge2elem from elem2edge
    for (int e = 0; e < n_elems; ++e) {
        int k = max_vpe;
        while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
        for (int i = 0; i < k; ++i) {
            int32_t eid = elem2edge[e * max_vpe + i];
            if (eid < 0) continue;
            if (edge2elem[eid * 2] == -1) edge2elem[eid * 2] = e;
            else edge2elem[eid * 2 + 1] = e;
        }
    }

    // Identify boundary edges (twin == -1 → only one elem)
    // active set of elements
    std::vector<bool> active_elem(n_elems, true);
    std::vector<int32_t> elem_layer(n_elems, -1);

    // active edge boundary count per element
    std::vector<int32_t> boundary_edge_count(n_elems, 0);
    for (int e = 0; e < n_elems; ++e) {
        int k = max_vpe;
        while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
        for (int i = 0; i < k; ++i) {
            int32_t eid = elem2edge[e * max_vpe + i];
            if (eid < 0) continue;
            bool is_boundary = (edge2elem[eid * 2] == -1 || edge2elem[eid * 2 + 1] == -1);
            if (is_boundary) boundary_edge_count[e]++;
        }
    }

    layers.clear();
    int layer_idx = 0;
    int remaining = n_elems;

    // For vert active count
    std::vector<int32_t> vert_active_elem_count(n_verts, 0);
    for (int e = 0; e < n_elems; ++e) {
        int k = max_vpe;
        while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
        for (int i = 0; i < k; ++i) {
            int32_t v = connectivity[e * max_vpe + i];
            if (v >= 0) vert_active_elem_count[v]++;
        }
    }

    while (remaining > 0) {
        // Find outer elements: boundary_edge_count > 0
        std::vector<int32_t> to_peel;
        for (int e = 0; e < n_elems; ++e) {
            if (active_elem[e] && boundary_edge_count[e] > 0) {
                to_peel.push_back(e);
            }
        }
        if (to_peel.empty()) {
            // no boundary elements — remaining are fully interior, force peel all
            for (int e = 0; e < n_elems; ++e) {
                if (active_elem[e]) to_peel.push_back(e);
            }
        }

        LayerData ld;
        // OE = outer elements (those being peeled)
        ld.OE = to_peel;

        // Identify OV: vertices on the outer boundary of this layer
        // = vertices that belong to peeled elems AND have an active boundary edge
        std::vector<bool> is_peeled_vert(n_verts, false);
        for (int e : to_peel) {
            int k = max_vpe;
            while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
            for (int i = 0; i < k; ++i) {
                int32_t v = connectivity[e * max_vpe + i];
                if (v >= 0) is_peeled_vert[v] = true;
            }
        }

        // OV: vertices on the boundary frontier of this layer (boundary edges).
        // bEdgeIDs: the boundary edges themselves (Python parity).
        std::unordered_map<int32_t, bool> ov_set;
        std::unordered_map<int32_t, bool> bedge_set;
        for (int e : to_peel) {
            int k = max_vpe;
            while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
            for (int i = 0; i < k; ++i) {
                int32_t eid = elem2edge[e * max_vpe + i];
                if (eid < 0) continue;
                bool is_boundary = (edge2elem[eid * 2] == -1 || edge2elem[eid * 2 + 1] == -1);
                if (is_boundary) {
                    int32_t v0 = edge2vert[eid * 2];
                    int32_t v1 = edge2vert[eid * 2 + 1];
                    ov_set[v0] = true;
                    ov_set[v1] = true;
                    bedge_set[eid] = true;
                }
            }
        }
        for (auto& [v, _] : ov_set) ld.OV.push_back(v);
        for (auto& [eid, _] : bedge_set) ld.bEdgeIDs.push_back(eid);

        // IE: inner elements = active, not-peeled elements that share an OV vertex
        // (Python parity: IE is adjacent to OV only, NOT to every OE vertex)
        std::vector<bool> is_peel_elem(n_elems, false);
        for (int e : to_peel) is_peel_elem[e] = true;
        std::vector<bool> is_ov_vert(n_verts, false);
        for (int32_t v : ld.OV) is_ov_vert[v] = true;
        for (int e = 0; e < n_elems; ++e) {
            if (!active_elem[e] || is_peel_elem[e]) continue;
            int k = max_vpe;
            while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
            bool adjacent = false;
            for (int i = 0; i < k && !adjacent; ++i) {
                int32_t v = connectivity[e * max_vpe + i];
                if (v >= 0 && is_ov_vert[v]) adjacent = true;
            }
            if (adjacent) ld.IE.push_back(e);
        }

        // IV: vertices of (OE ∪ IE) not in OV (Python parity — IE vertices count).
        // Must run AFTER IE so vertices from inner elements are included.
        for (int e : ld.IE) {
            int k = max_vpe;
            while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
            for (int i = 0; i < k; ++i) {
                int32_t v = connectivity[e * max_vpe + i];
                if (v >= 0) is_peeled_vert[v] = true;
            }
        }
        for (int v = 0; v < n_verts; ++v) {
            if (is_peeled_vert[v] && ov_set.find(v) == ov_set.end()) {
                ld.IV.push_back(v);
            }
        }

        // Python parity: each skeletonization iteration consumes BOTH the outer
        // ring (OE) and the inner vertex-adjacent ring (IE). Build a single
        // peel list combining both, then update inactive state once.
        std::vector<int32_t> to_consume;
        to_consume.reserve(ld.OE.size() + ld.IE.size());
        for (int32_t e : ld.OE) to_consume.push_back(e);
        for (int32_t e : ld.IE) to_consume.push_back(e);

        layers.push_back(std::move(ld));

        // Mark all consumed elements inactive, update boundary counts
        for (int e : to_consume) {
            active_elem[e] = false;
            elem_layer[e] = layer_idx;
            remaining--;

            // Update vertex active counts
            int k = max_vpe;
            while (k > 3 && connectivity[e * max_vpe + k - 1] == -1) --k;
            for (int i = 0; i < k; ++i) {
                int32_t v = connectivity[e * max_vpe + i];
                if (v >= 0) vert_active_elem_count[v]--;
            }

            // For each edge of e, find the neighbor and update its boundary count
            for (int i = 0; i < k; ++i) {
                int32_t eid = elem2edge[e * max_vpe + i];
                if (eid < 0) continue;
                int32_t nb = -1;
                if (edge2elem[eid * 2] == e) nb = edge2elem[eid * 2 + 1];
                else if (edge2elem[eid * 2 + 1] == e) nb = edge2elem[eid * 2];
                if (nb >= 0 && active_elem[nb]) {
                    boundary_edge_count[nb]++;
                    // Mark this edge as now boundary
                    edge2elem[eid * 2] = -1;
                    edge2elem[eid * 2 + 1] = -1;
                }
            }
        }

        layer_idx++;
    }
}

std::vector<double> HalfEdgeMesh::signed_areas() const {
    std::vector<double> areas(n_elems, 0.0);
    for (int e = 0; e < n_elems; ++e) {
        int k = max_vpe;
        const int32_t* row = connectivity.data() + e * max_vpe;
        while (k > 3 && row[k-1] == -1) --k;
        // Shoelace formula for polygon
        double area = 0.0;
        for (int i = 0; i < k; ++i) {
            int32_t v0 = row[i];
            int32_t v1 = row[(i + 1) % k];
            area += px[v0] * py[v1] - px[v1] * py[v0];
        }
        areas[e] = area * 0.5;
    }
    return areas;
}

} // namespace chilmesh
