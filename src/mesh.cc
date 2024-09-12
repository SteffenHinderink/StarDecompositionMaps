#include "mesh.h"

Eigen::Vector3d hsv_to_rgb(const Eigen::Vector3d& hsv) {
    int h = (int) (hsv(0) / (M_PI / 3)) % 6;
    double f = hsv(0) / (M_PI / 3) - h;
    double p = hsv(2) * (1 - hsv(1));
    double q = hsv(2) * (1 - hsv(1) * f);
    double t = hsv(2) * (1 - hsv(1) * (1 - f));
    switch (h) {
    case 0:
        return {hsv(2), t, p};
    case 1:
        return {q, hsv(2), p};
    case 2:
        return {p, hsv(2), t};
    case 3:
        return {p, q, hsv(2)};
    case 4:
        return {t, p, hsv(2)};
    case 5:
        return {hsv(2), p, q};
    }
    return {0, 0, 0};
}

Eigen::Vector3d normalize_color(const Eigen::Vector3d& color) {
    return ((color * 2 - Eigen::Vector3d::Ones()).normalized() + Eigen::Vector3d::Ones()) / 2;
}

#ifdef TRI

std::array<OpenMesh::VertexHandle, 2>
#ifdef TET
TriMesh
#else
Mesh
#endif
    ::edge_vertices(OpenMesh::EdgeHandle e) {
    OpenMesh::HalfedgeHandle h = halfedge_handle(e, 0);
    return {from_vertex_handle(h), to_vertex_handle(h)};
}

std::array<OpenMesh::VertexHandle, 2> PolyMesh::edge_vertices(OpenMesh::EdgeHandle e) {
    OpenMesh::HalfedgeHandle h = halfedge_handle(e, 0);
    return {from_vertex_handle(h), to_vertex_handle(h)};
}

#endif

#ifdef TET

bool exactinitialized = false;

double orient3d(const std::vector<Eigen::Vector3d>& _tetrahedron) {
    std::vector<Eigen::Vector3d> tetrahedron = _tetrahedron;
    if (!exactinitialized) {
        exactinit();
        exactinitialized = true;
    }
    return orient3d(tetrahedron[0].data(), tetrahedron[1].data(), tetrahedron[2].data(), tetrahedron[3].data());
}

Halfedge Mesh::halfedge_opposite_halfedge(Halfedge e, ::Cell c) {
    Vertex a = from_vertex_handle(e);
    Vertex b = to_vertex_handle(e);
    std::vector<Vertex> opposite_vertices;
    for (auto cv : tet_vertices(c)) {
        if (cv != a && cv != b) {
            opposite_vertices.push_back(cv);
        }
    }
    if (halfface_opposite_vertex(find_halfface({a, b, opposite_vertices[0]})) != opposite_vertices[1]) {
        return find_halfedge(opposite_vertices[1], opposite_vertices[0]);
    }
    return find_halfedge(opposite_vertices[0], opposite_vertices[1]);
}

Halfface Mesh::vertex_opposite_halfface(Vertex v, ::Cell c) {
    for (auto ch : cell_halffaces(c)) {
        bool has_v = false;
        for (auto chv : halfface_vertices(ch)) {
            if (chv == v) {
                has_v = true;
                break;
            }
        }
        if (!has_v) {
            return ch;
        }
    }
    return Halfface(-1);
}

void Mesh::remove_cell(::Cell c) {
    std::vector<Vertex> vertices;
    for (auto cv : tet_vertices(c)) {
        vertices.push_back(cv);
    }
    std::vector<::Edge> edges;
    for (auto ce : cell_edges(c)) {
        edges.push_back(ce);
    }
    std::vector<::Face> faces;
    for (auto cf : cell_faces(c)) {
        faces.push_back(cf);
    }
    delete_cell(c);
    for (auto f : faces) {
        if (is_boundary(halfface_handle(f, 0)) && is_boundary(halfface_handle(f, 1))) {
            delete_face(f);
        }
    }
    for (auto e : edges) {
        if (!ef_iter(e).is_valid()) {
            delete_edge(e);
        }
    }
    for (auto v : vertices) {
        if (!ve_iter(v).is_valid()) {
            delete_vertex(v);
        }
    }
}

void Mesh::remove_cells(std::vector<::Cell> cells) {
    std::set<::Face> faces;
    for (auto c : cells) {
        for (auto f : cell_faces(c)) {
            faces.insert(f);
        }
        delete_cell(c);
    }
    std::set<::Edge> edges;
    for (auto f : faces) {
        if (is_boundary(halfface_handle(f, 0)) && is_boundary(halfface_handle(f, 1))) {
            for (auto e : face_edges(f)) {
                edges.insert(e);
            }
            delete_face(f);
        }
    }
    std::set<Vertex> vertices;
    for (auto e : edges) {
        if (!ef_iter(e).is_valid()) {
            for (auto v : edge_vertices(e)) {
                vertices.insert(v);
            }
            delete_edge(e);
        }
    }
    for (auto v : vertices) {
        if (!ve_iter(v).is_valid()) {
            delete_vertex(v);
        }
    }
}

Vertex Mesh::split_tet(::Cell c, const Eigen::Vector3d& p) {
    Vertex v = add_vertex(p);
    std::vector<std::vector<Vertex>> tets;
    for (auto ch : cell_halffaces(c)) {
        std::vector<Vertex> tet;
        for (auto chv : halfface_vertices(ch)) {
            tet.push_back(chv);
        }
        tet.push_back(v);
        tets.push_back(tet);
    }
    delete_cell(c);
    for (auto tet : tets) {
        add_cell(tet, true);
    }
    return v;
}

std::vector<Vertex> Mesh::boundary_vertices() {
    std::vector<Vertex> boundary_vertices;
    for (auto v : vertices()) {
        if (is_boundary(v)) {
            boundary_vertices.push_back(v);
        }
    }
    return boundary_vertices;
}

std::vector<::Edge> Mesh::boundary_edges() {
    std::vector<::Edge> boundary_edges;
    for (auto e : edges()) {
        if (is_boundary(e)) {
            boundary_edges.push_back(e);
        }
    }
    return boundary_edges;
}

std::vector<Halfface> Mesh::boundary_halffaces() {
    std::vector<Halfface> boundary_halffaces;
    for (auto h : halffaces()) {
        if (is_boundary(h)) {
            boundary_halffaces.push_back(h);
        }
    }
    return boundary_halffaces;
}

Halfface Mesh::boundary_halfface(Halfedge e) {
    for (auto eh : halfedge_halffaces(e)) {
        if (is_boundary(eh)) {
            return eh;
        }
    }
    return Halfface(-1);
}

Halfedge Mesh::boundary_next(Halfedge e) {
    return next_halfedge_in_halfface(e, boundary_halfface(e));
}

Halfedge Mesh::boundary_out(Vertex v) {
    for (auto ve : vertex_edges(v)) {
        if (is_boundary(ve)) {
            Halfedge e = halfedge_handle(ve, 0);
            if (from_vertex_handle(e) == v) {
                return e;
            } else {
                return opposite_halfedge_handle(e);
            }
        }
    }
    return Halfedge(-1);
}

bool Mesh::degenerate_or_inverted(::Cell c) {
    std::vector<Eigen::Vector3d> tetrahedron;
    for (auto cv : tet_vertices(c)) {
        tetrahedron.push_back(position(cv));
    }
    return orient3d(tetrahedron) >= 0;
}

Mesh read(std::string filename, int& orientation) {
    Mesh mesh;
    std::string extension = filename.substr(filename.rfind('.') + 1);
    if (extension == "ovm") {
        OpenVolumeMesh::IO::FileManager fm;
        fm.setVerbosityLevel(1);
        fm.readFile(filename, mesh);
        return mesh;
    }
    std::vector<std::vector<std::string>> input;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> split_line;
        std::istringstream ss(line);
        std::string s;
        while (std::getline(ss, s, ' ')) {
            if (s != "") {
                split_line.push_back(s);
            }
        }
        input.push_back(split_line);
    }
    std::vector<Eigen::Vector3d> positions;
    std::vector<std::vector<int>> indices;
    if (extension == "vtk") {
        int vertices_idx, n_vertices, cells_idx, n_cells;
        for (int i = 0; i < input.size(); i++) {
            if (input[i].size() > 1 && input[i][0] == "POINTS") {
                vertices_idx = i + 1;
                n_vertices = std::stoi(input[i][1]);
            }
            if (input[i].size() > 1 && input[i][0] == "CELLS") {
                cells_idx = i + 1;
                n_cells = std::stoi(input[i][1]);
            }
        }
        if (input[vertices_idx].size() == 3) {
            for (int i = vertices_idx; i < vertices_idx + n_vertices; i++) {
                positions.push_back(
                    {std::stod(input[i][0]),
                     std::stod(input[i][1]),
                     std::stod(input[i][2])});
            }
        } else {
            for (int i = 0; i < n_vertices; i++) {
                positions.push_back(
                    {std::stod(input[vertices_idx][3 * i]),
                     std::stod(input[vertices_idx][3 * i + 1]),
                     std::stod(input[vertices_idx][3 * i + 2])});
            }
        }
        if (input[cells_idx][0] == "OFFSETS") {
            n_cells--;
            int offset = cells_idx + n_cells + 3;
            for (int i = 0; i < n_cells; i++) {
                indices.push_back(
                    {std::stoi(input[offset + 4 * i][0]),
                     std::stoi(input[offset + 4 * i + 1][0]),
                     std::stoi(input[offset + 4 * i + 2][0]),
                     std::stoi(input[offset + 4 * i + 3][0])});
            }
        } else {
            if (input[cells_idx].size() == 5) {
                for (int i = cells_idx; i < cells_idx + n_cells; i++) {
                    indices.push_back(
                        {std::stoi(input[i][1]),
                         std::stoi(input[i][2]),
                         std::stoi(input[i][3]),
                         std::stoi(input[i][4])});
                }
            } else {
                for (int i = 0; i < n_cells; i++) {
                    indices.push_back(
                        {std::stoi(input[cells_idx + 5 * i + 1][0]),
                         std::stoi(input[cells_idx + 5 * i + 2][0]),
                         std::stoi(input[cells_idx + 5 * i + 3][0]),
                         std::stoi(input[cells_idx + 5 * i + 4][0])});
                }
            }
        }
    } else if (extension == "mesh") {
        int vertices_idx, n_vertices, cells_idx, n_cells;
        for (int i = 0; i < input.size(); i++) {
            if (input[i].size() > 0 && input[i][0] == "Vertices") {
                vertices_idx = i + 2;
                n_vertices = std::stoi(input[i + 1][0]);
            }
            if (input[i].size() > 0 && input[i][0] == "Tetrahedra") {
                cells_idx = i + 2;
                n_cells = std::stoi(input[i + 1][0]);
            }
        }
        for (int i = vertices_idx; i < vertices_idx + n_vertices; i++) {
            positions.push_back(
                {std::stod(input[i][0]),
                 std::stod(input[i][1]),
                 std::stod(input[i][2])});
        }
        for (int i = cells_idx; i < cells_idx + n_cells; i++) {
            indices.push_back(
                {std::stoi(input[i][0]),
                 std::stoi(input[i][1]),
                 std::stoi(input[i][2]),
                 std::stoi(input[i][3])});
        }
    } else if (extension == "msh") {
        int vertices_idx, n_vertices, cells_idx, n_cells;
        for (int i = 0; i < input.size(); i++) {
            if (input[i].size() > 0 && input[i][0] == "$Nodes") {
                vertices_idx = i + 2;
                n_vertices = std::stoi(input[i + 1][0]);
            }
            if (input[i].size() > 0 && input[i][0] == "$Elements") {
                cells_idx = i + 2;
                n_cells = std::stoi(input[i + 1][0]);
            }
        }
        for (int i = vertices_idx; i < vertices_idx + n_vertices; i++) {
            positions.push_back(
                {std::stod(input[i][1]),
                 std::stod(input[i][2]),
                 std::stod(input[i][3])});
        }
        for (int i = cells_idx; i < cells_idx + n_cells; i++) {
            indices.push_back(
                {std::stoi(input[i][3]),
                 std::stoi(input[i][4]),
                 std::stoi(input[i][5]),
                 std::stoi(input[i][6])});
        }
    } else if (extension == "tet") {
        int n_vertices, n_cells;
        for (int i = 0; i < input.size(); i++) {
            if (input[i].size() > 1 && input[i][1] == "vertices") {
                n_vertices = std::stoi(input[i][0]);
            }
            if (input[i].size() > 1 && input[i][1] == "tets") {
                n_cells = std::stoi(input[i][0]);
            }
        }
        for (int i = 2; i < 2 + n_vertices; i++) {
            positions.push_back(
                {std::stod(input[i][0]),
                 std::stod(input[i][1]),
                 std::stod(input[i][2])});
        }
        for (int i = 2 + n_vertices; i < 2 + n_vertices + n_cells; i++) {
            indices.push_back(
                {std::stoi(input[i][1]),
                 std::stoi(input[i][2]),
                 std::stoi(input[i][3]),
                 std::stoi(input[i][4])});
        }
    }
    int min_idx = 1;
    for (int i = 0; i < indices.size() && min_idx > 0; i++) {
        for (int j = 0; j < 4 && min_idx > 0; j++) {
            if (indices[i][j] == 0) {
                min_idx = 0;
            }
        }
    }
    if (orientation == 0) {
        std::vector<Eigen::Vector3d> tetrahedron(4);
        for (int i = 0; i < 4; i++) {
            tetrahedron[i] = positions[indices[0][i] - min_idx];
        }
        orientation = orient3d(tetrahedron) < 0 ? 1 : -1;
    }
    std::vector<Vertex> vertices;
    for (int i = 0; i < positions.size(); i++) {
        vertices.push_back(mesh.add_vertex(positions[i]));
    }
    for (int i = 0; i < indices.size(); i++) {
        mesh.add_cell(
            vertices[indices[i][0] - min_idx],
            vertices[indices[i][1] - min_idx],
            vertices[indices[i][orientation > 0 ? 2 : 3] - min_idx],
            vertices[indices[i][orientation > 0 ? 3 : 2] - min_idx], true);
    }
    return mesh;
}

Mesh read(std::string filename) {
    int orientation = 0;
    return read(filename, orientation);
}

#endif
