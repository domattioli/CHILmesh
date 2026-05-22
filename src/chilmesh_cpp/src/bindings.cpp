#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "halfedge.hpp"

namespace py = pybind11;
using namespace chilmesh;

// Helper: build mesh from numpy arrays
static HalfEdgeMesh* _build(
        py::array_t<double, py::array::c_style | py::array::forcecast> points,
        py::array_t<int32_t, py::array::c_style | py::array::forcecast> connectivity) {

    auto pts = points.unchecked<2>();
    auto conn = connectivity.unchecked<2>();

    int n_pts = (int)pts.shape(0);
    int n_el  = (int)conn.shape(0);
    int vpe   = (int)conn.shape(1);

    // Pack points as x0..xN y0..yN for build()
    std::vector<double> flat_pts(n_pts * 2);
    for (int i = 0; i < n_pts; ++i) {
        flat_pts[i] = pts(i, 0);
        flat_pts[n_pts + i] = pts(i, 1);
    }

    // Flatten connectivity
    std::vector<int32_t> flat_conn(n_el * vpe);
    for (int e = 0; e < n_el; ++e)
        for (int j = 0; j < vpe; ++j)
            flat_conn[e * vpe + j] = conn(e, j);

    auto* mesh = new HalfEdgeMesh();
    mesh->build(flat_pts.data(), n_pts, flat_conn.data(), n_el, vpe);
    return mesh;
}

PYBIND11_MODULE(chilmesh_cpp, m) {
    m.doc() = "CHILmesh C++ half-edge backend";

    py::class_<HalfEdgeMesh>(m, "CppMesh")
        .def_readonly("n_verts", &HalfEdgeMesh::n_verts)
        .def_readonly("n_elems", &HalfEdgeMesh::n_elems)
        .def_readonly("n_edges", &HalfEdgeMesh::n_edges)
        .def_readonly("adjacency_built", &HalfEdgeMesh::adjacency_built)
        .def_property_readonly("n_layers", [](const HalfEdgeMesh& m) {
            return (int)m.layers.size();
        })
        .def_property_readonly("edge2vert", [](const HalfEdgeMesh& m) {
            if (m.edge2vert.empty()) {
                std::vector<py::ssize_t> sh = {0, 2};
                return py::array_t<int32_t>(sh);
            }
            return py::array_t<int32_t>(
                {(py::ssize_t)m.n_edges, (py::ssize_t)2},
                m.edge2vert.data());
        })
        .def_property_readonly("elem2edge", [](const HalfEdgeMesh& m) {
            if (m.elem2edge.empty()) {
                std::vector<py::ssize_t> sh = {0, (py::ssize_t)m.max_vpe};
                return py::array_t<int32_t>(sh);
            }
            return py::array_t<int32_t>(
                {(py::ssize_t)m.n_elems, (py::ssize_t)m.max_vpe},
                m.elem2edge.data());
        })
        .def_property_readonly("vert2edge", [](const HalfEdgeMesh& m) {
            py::list result;
            for (const auto& ve : m.vert2edge) {
                result.append(py::array_t<int32_t>(ve.size(), ve.data()));
            }
            return result;
        })
        .def_property_readonly("layers", [](const HalfEdgeMesh& m) {
            py::list result;
            for (const auto& ld : m.layers) {
                py::dict layer;
                layer["OE"] = py::array_t<int32_t>(ld.OE.size(), ld.OE.data());
                layer["IE"] = py::array_t<int32_t>(ld.IE.size(), ld.IE.data());
                layer["OV"] = py::array_t<int32_t>(ld.OV.size(), ld.OV.data());
                layer["IV"] = py::array_t<int32_t>(ld.IV.size(), ld.IV.data());
                layer["bEdgeIDs"] = py::array_t<int32_t>(ld.bEdgeIDs.size(), ld.bEdgeIDs.data());
                result.append(layer);
            }
            return result;
        });

    // fast_init: build half-edge + CCW check, no adjacency
    m.def("fast_init", [](
            py::array_t<double, py::array::c_style | py::array::forcecast> points,
            py::array_t<int32_t, py::array::c_style | py::array::forcecast> connectivity) {
        return _build(points, connectivity);
    }, py::return_value_policy::take_ownership,
       "Build half-edge mesh (no adjacency). Returns CppMesh.");

    // full_init: build + adjacency + skeletonize
    m.def("full_init", [](
            py::array_t<double, py::array::c_style | py::array::forcecast> points,
            py::array_t<int32_t, py::array::c_style | py::array::forcecast> connectivity) {
        auto* mesh = _build(points, connectivity);
        mesh->build_adjacency();
        mesh->skeletonize();
        return mesh;
    }, py::return_value_policy::take_ownership,
       "Build half-edge mesh with adjacency and skeletonization. Returns CppMesh.");

    // quality_analysis: signed areas
    m.def("quality_analysis", [](const HalfEdgeMesh& mesh) {
        auto areas = mesh.signed_areas();
        return py::array_t<double>(areas.size(), areas.data());
    }, "Compute signed area for each element. Returns numpy float64 array.");

    // get_vertex_edges
    m.def("get_vertex_edges", [](const HalfEdgeMesh& mesh, int vertex_id) {
        if (!mesh.adjacency_built)
            throw std::runtime_error("Adjacency not built. Use full_init.");
        if (vertex_id < 0 || vertex_id >= mesh.n_verts)
            throw std::out_of_range("vertex_id out of range");
        const auto& ve = mesh.vert2edge[vertex_id];
        return py::array_t<int32_t>(ve.size(), ve.data());
    }, "Return edge IDs incident to vertex.");
}
