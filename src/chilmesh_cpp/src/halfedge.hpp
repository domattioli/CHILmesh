#pragma once
#include <vector>
#include <cstdint>
#include <string>
#include <stdexcept>

namespace chilmesh {

// Half-edge structure — flat array storage for cache efficiency
struct HalfEdge {
    int32_t twin;  // opposite half-edge (-1 if boundary)
    int32_t next;  // next half-edge in face loop
    int32_t prev;  // prev half-edge in face loop
    int32_t vert;  // origin vertex
    int32_t face;  // owning face (-1 if boundary)
};

struct LayerData {
    std::vector<int32_t> OE;         // outer elements
    std::vector<int32_t> IE;         // inner elements
    std::vector<int32_t> OV;         // outer vertices
    std::vector<int32_t> IV;         // inner vertices
    std::vector<int32_t> bEdgeIDs;   // boundary edges defining this layer's outer frontier
};

class HalfEdgeMesh {
public:
    // Geometry
    std::vector<double> px, py;  // vertex coords (flat, cache-friendly)
    int32_t n_verts = 0;
    int32_t n_elems = 0;
    int32_t n_edges = 0;

    // Half-edge topology
    std::vector<HalfEdge> half_edges;

    // CSR-style: first half-edge out of each vertex
    std::vector<int32_t> vert_he;   // size n_verts, -1 if none

    // first half-edge of each face
    std::vector<int32_t> face_he;   // size n_elems

    // elem2vert (connectivity): n_elems x max_vpe  (-1 for triangles in 4-col)
    std::vector<int32_t> connectivity;  // row-major [n_elems * max_vpe]
    int32_t max_vpe = 3;  // vertices per element (3 or 4)

    // Adjacency (populated in full_init)
    std::vector<int32_t> edge2vert;   // [n_edges * 2]
    std::vector<int32_t> elem2edge;   // [n_elems * max_vpe]
    std::vector<std::vector<int32_t>> vert2edge;  // jagged

    // Skeletonization layers
    std::vector<LayerData> layers;

    bool adjacency_built = false;

    // Build from points+connectivity arrays
    void build(const double* pts, int n_pts,
               const int32_t* conn, int n_el, int vpe);

    // Build adjacency structures
    void build_adjacency();

    // Skeletonize (layer peeling)
    void skeletonize();

    // Quality analysis: returns signed area for each element
    std::vector<double> signed_areas() const;
};

} // namespace chilmesh
