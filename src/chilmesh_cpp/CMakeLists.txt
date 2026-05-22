cmake_minimum_required(VERSION 3.15)
project(chilmesh_cpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 REQUIRED)

pybind11_add_module(chilmesh_cpp
    src/halfedge.cpp
    src/bindings.cpp
)

target_compile_options(chilmesh_cpp PRIVATE -O3 -march=native)
target_include_directories(chilmesh_cpp PRIVATE src)
