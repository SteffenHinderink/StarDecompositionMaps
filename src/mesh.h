#pragma once

#include <Eigen/Dense>

// >--- Miscellaneous utility functions
const char SEP =
#ifdef _WIN32
    '\\'
#else
    '/'
#endif
    ;

Eigen::Vector3d hsv_to_rgb(const Eigen::Vector3d& hsv);

Eigen::Vector3d normalize_color(const Eigen::Vector3d& color);
// ---<

#if !defined(TRI) && !defined(TET)
#define TRI
#define TET
#endif

#ifdef TRI

#define _USE_MATH_DEFINES
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
// TriMesh_ArrayKernelT must be included before EigenVectorT
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#ifdef TET
using SVertex = OpenMesh::VertexHandle;
using SHalfedge = OpenMesh::HalfedgeHandle;
using SEdge = OpenMesh::EdgeHandle;
using SFace = OpenMesh::FaceHandle;
template <typename E, typename T>
using SProperty = OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type>;
#else
using Vertex = OpenMesh::VertexHandle;
using Halfedge = OpenMesh::HalfedgeHandle;
using Edge = OpenMesh::EdgeHandle;
using Face = OpenMesh::FaceHandle;
template <typename E, typename T>
using Property = OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type>;
#endif

struct Traits : public OpenMesh::DefaultTraits {
    using Point = Eigen::Vector3d;
};

#ifdef TET
class TriMesh
#else
class Mesh
#endif
    : public OpenMesh::TriMesh_ArrayKernelT<Traits> {
public:
    Eigen::Vector3d position(OpenMesh::VertexHandle v) { return point(v); }

    void set_position(OpenMesh::VertexHandle v, const Eigen::Vector3d& p) { set_point(v, p); }

    OpenMesh::PolyConnectivity::ConstFaceVertexCCWRange face_vertices(OpenMesh::FaceHandle f) { return fv_ccw_range(f); }

    using OpenMesh::TriMesh_ArrayKernelT<Traits>::property;
    template <typename E, typename T>
    OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type> property(std::string key = "", const T& init = T()) {
        if (key == "") {
            OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type> p = OpenMesh::makeTemporaryProperty<E, T>(*this);
            if (init != T()) {
                fill_property<E, T>(p, init);
            }
            return p;
        }
        if (OpenMesh::hasProperty<E, T>(*this, key.c_str())) {
            return OpenMesh::getProperty<E, T>(*this, key.c_str());
        }
        OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type> p = OpenMesh::getOrMakeProperty<E, T>(*this, key.c_str());
        if (init != T()) {
            fill_property<E, T>(p, init);
        }
        return p;
    }

    std::array<OpenMesh::VertexHandle, 2> edge_vertices(OpenMesh::EdgeHandle e);

private:
    template <typename E, typename T>
    void fill_property(OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type>& p, const T& value) {
        if constexpr (std::is_same<E, OpenMesh::VertexHandle>::value) {
            for (auto v : vertices()) {
                p[v] = value;
            }
        }
        if constexpr (std::is_same<E, OpenMesh::HalfedgeHandle>::value) {
            for (auto h : halfedges()) {
                p[h] = value;
            }
        }
        if constexpr (std::is_same<E, OpenMesh::EdgeHandle>::value) {
            for (auto e : edges()) {
                p[e] = value;
            }
        }
        if constexpr (std::is_same<E, OpenMesh::FaceHandle>::value) {
            for (auto f : faces()) {
                p[f] = value;
            }
        }
    }
};

class PolyMesh : public OpenMesh::PolyMesh_ArrayKernelT<Traits> {
public:
    Eigen::Vector3d position(OpenMesh::VertexHandle v) { return point(v); }

    void set_position(OpenMesh::VertexHandle v, const Eigen::Vector3d& p) { set_point(v, p); }

    OpenMesh::PolyConnectivity::ConstFaceVertexCCWRange face_vertices(OpenMesh::FaceHandle f) { return fv_ccw_range(f); }

    using OpenMesh::PolyMesh_ArrayKernelT<Traits>::property;
    template <typename E, typename T>
    OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type> property(std::string key = "", const T& init = T()) {
        if (key == "") {
            OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type> p = OpenMesh::makeTemporaryProperty<E, T>(*this);
            if (init != T()) {
                fill_property<E, T>(p, init);
            }
            return p;
        }
        if (OpenMesh::hasProperty<E, T>(*this, key.c_str())) {
            return OpenMesh::getProperty<E, T>(*this, key.c_str());
        }
        OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type> p = OpenMesh::getOrMakeProperty<E, T>(*this, key.c_str());
        if (init != T()) {
            fill_property<E, T>(p, init);
        }
        return p;
    }

    std::array<OpenMesh::VertexHandle, 2> edge_vertices(OpenMesh::EdgeHandle e);

private:
    template <typename E, typename T>
    void fill_property(OpenMesh::PropertyManager<typename OpenMesh::HandleToPropHandle<E, T>::type>& p, const T& value) {
        if constexpr (std::is_same<E, OpenMesh::VertexHandle>::value) {
            for (auto v : vertices()) {
                p[v] = value;
            }
        }
        if constexpr (std::is_same<E, OpenMesh::HalfedgeHandle>::value) {
            for (auto h : halfedges()) {
                p[h] = value;
            }
        }
        if constexpr (std::is_same<E, OpenMesh::EdgeHandle>::value) {
            for (auto e : edges()) {
                p[e] = value;
            }
        }
        if constexpr (std::is_same<E, OpenMesh::FaceHandle>::value) {
            for (auto f : faces()) {
                p[f] = value;
            }
        }
    }
};

#endif

#ifdef TET

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#ifdef tetgenH
#include <tetgen.h>
#else
#include <predicates.h>
#endif

double orient3d(const std::vector<Eigen::Vector3d>& tetrahedron);

using Vertex = OpenVolumeMesh::VertexHandle;
using Halfedge = OpenVolumeMesh::HalfEdgeHandle;
using Edge = OpenVolumeMesh::EdgeHandle;
using Halfface = OpenVolumeMesh::HalfFaceHandle;
using Face = OpenVolumeMesh::FaceHandle;
using Cell = OpenVolumeMesh::CellHandle;
template <typename E, typename T>
using Property = OpenVolumeMesh::PropertyPtr<T, typename E::EntityTag>;

class Vector3dWrapper : public Eigen::Vector3d {
    using Eigen::Vector3d::Vector3d;

public:
    Vector3dWrapper(int) : Eigen::Vector3d() {}
};

class Mesh : public OpenVolumeMesh::TetrahedralGeometryKernel<Vector3dWrapper> {
public:
    Eigen::Vector3d position(Vertex v) { return vertex(v); }

    void set_position(Vertex v, const Eigen::Vector3d& p) { set_vertex(v, p); }

    template <typename E, typename T>
    Property<E, T> property(std::string key = "", const T& init = T()) {
        Property<E, T> p;
        if (key == "") {
            p = request_property<T, typename E::EntityTag>(std::to_string(n_tmp_propertys++), init);
        } else {
            if (property_exists<T, typename E::EntityTag>(key)) {
                return request_property<T, typename E::EntityTag>(key);
            }
            p = request_property<T, typename E::EntityTag>(key, init);
            set_persistent(p);
        }
        return p;
    }

    Halfedge halfedge_opposite_halfedge(Halfedge e, ::Cell c);

    Halfface vertex_opposite_halfface(Vertex v, ::Cell c);

    void remove_cell(::Cell c);

    void remove_cells(std::vector<::Cell> cells);

    Vertex split_tet(::Cell c, const Eigen::Vector3d& p);

    std::vector<Vertex> boundary_vertices();

    std::vector<::Edge> boundary_edges();

    std::vector<Halfface> boundary_halffaces();

    Halfface boundary_halfface(Halfedge e);

    Halfedge boundary_next(Halfedge e);

    Halfedge boundary_out(Vertex v);

    bool degenerate_or_inverted(::Cell c);

private:
    int n_tmp_propertys;
};

Mesh read(std::string filename, int& orientation);

Mesh read(std::string filename);

#endif
