#include "sdm.h"

void forbid(Mesh& mesh, std::vector<Vertex>& loop) {
    // Forbid boundary
    auto f_forbidden = mesh.property<Face, bool>("f_forbidden");
    auto e_forbidden = mesh.property<Edge, bool>("e_forbidden");
    auto v_forbidden = mesh.property<Vertex, bool>("v_forbidden");
    auto bdr = mesh.property<Halfface, bool>("bdr");
    for (auto f : mesh.faces()) {
        f_forbidden[f] = false;
        if (mesh.is_boundary(f) || bdr[mesh.halfface_handle(f, 0)] || bdr[mesh.halfface_handle(f, 1)]) {
            f_forbidden[f] = true;
        }
    }
    for (auto e : mesh.edges()) {
        e_forbidden[e] = false;
        for (auto ef : mesh.edge_faces(e)) {
            if (f_forbidden[ef]) {
                e_forbidden[e] = true;
                break;
            }
        }
    }
    for (auto v : mesh.vertices()) {
        v_forbidden[v] = false;
        for (auto vf : mesh.vertex_faces(v)) {
            if (f_forbidden[vf]) {
                v_forbidden[v] = true;
                break;
            }
        }
    }
    // Allow loop
    for (int i = 0; i < loop.size(); i++) {
        e_forbidden[mesh.edge_handle(mesh.find_halfedge(loop[i], loop[(i + 1) % loop.size()]))] = false;
        v_forbidden[loop[i]] = false;
    }
}

void fill_sides(Mesh& mesh, std::vector<Vertex>& loop, std::vector<std::set<Halfface>>& sides) {
    auto f_forbidden = mesh.property<Face, bool>("f_forbidden");
    auto e_forbidden = mesh.property<Edge, bool>("e_forbidden");
    auto v_forbidden = mesh.property<Vertex, bool>("v_forbidden");
    auto border = mesh.property<Edge, bool>("", false);
    for (int i = 0; i < 2; i++) {
        sides.push_back(std::set<Halfface>());
    }
    for (int i = 0; i < loop.size(); i++) {
        Halfedge e = mesh.find_halfedge(loop[i], loop[(i + 1) % loop.size()]);
        border[mesh.edge_handle(e)] = true;
        for (auto eh : mesh.halfedge_halffaces(e)) {
            if (mesh.is_boundary(eh)) {
                sides[0].insert(eh);
                break;
            }
        }
        e = mesh.opposite_halfedge_handle(e);
        for (auto eh : mesh.halfedge_halffaces(e)) {
            if (mesh.is_boundary(eh)) {
                sides[1].insert(eh);
                break;
            }
        }
    }
    auto fill = [&](std::set<Halfface>& side) {
        std::queue<Halfface> queue;
        auto enqueued = mesh.property<Halfface, bool>("", false);
        for (auto h : side) {
            queue.push(h);
            enqueued[h] = true;
        }
        while (!queue.empty()) {
            Halfface h = queue.front();
            queue.pop();
            side.insert(h);
            for (auto he : mesh.halfface_halfedges(h)) {
                if (!border[mesh.edge_handle(he)]) {
                    // Prioritize inner
                    Halfface hh(-1);
                    for (auto heh : mesh.halfedge_halffaces(he)) {
                        Face f = mesh.face_handle(heh);
                        if (f != mesh.face_handle(h) && !mesh.is_boundary(f) && f_forbidden[f]) {
                            hh = mesh.opposite_halfface_handle(heh);
                            break;
                        }
                    }
                    if (!mesh.is_valid(hh) && mesh.is_boundary(h)) {
                        hh = mesh.boundary_halfface(mesh.opposite_halfedge_handle(he));
                    }
                    if (mesh.is_valid(hh) && !enqueued[hh]) {
                        queue.push(hh);
                        enqueued[hh] = true;
                    }
                }
            }
        }
    };
    for (int i = 0; i < 2; i++) {
        fill(sides[i]);
    }
}

void pre_split(Mesh& mesh, std::vector<std::set<Halfface>>& sides) {
    auto f_forbidden = mesh.property<Face, bool>("f_forbidden");
    auto e_forbidden = mesh.property<Edge, bool>("e_forbidden");
    auto v_forbidden = mesh.property<Vertex, bool>("v_forbidden");
    auto Q = mesh.property<Vertex, Vector3q>("Q");
    auto vside = mesh.property<Vertex, int>("", 0);
    auto eside = mesh.property<Edge, int>("", 0);
    auto fside = mesh.property<Face, int>("", 0);
    for (auto h : sides[0]) {
        for (auto hv : mesh.halfface_vertices(h)) {
            if (v_forbidden[hv]) {
                vside[hv] = -1;
            }
        }
        for (auto he : mesh.halfface_edges(h)) {
            if (e_forbidden[he]) {
                eside[he] = -1;
            }
        }
        Face f = mesh.face_handle(h);
        ASSERT(f_forbidden[f]);
        fside[f] = -1;
    }
    for (auto h : sides[1]) {
        for (auto hv : mesh.halfface_vertices(h)) {
            if (v_forbidden[hv]) {
                vside[hv] = 1;
            }
        }
        for (auto he : mesh.halfface_edges(h)) {
            if (e_forbidden[he]) {
                eside[he] = 1;
            }
        }
        Face f = mesh.face_handle(h);
        ASSERT(f_forbidden[f] && fside[f] == 0);
        fside[f] = 1;
    }
    std::vector<Edge> edges;
    for (auto e : mesh.edges()) {
        if (e_forbidden[e]) {
            continue;
        }
        auto [v0, v1] = mesh.edge_vertices(e);
        if (vside[v0] * vside[v1] == -1) {
            edges.push_back(e);
        }
    }
    for (auto e : edges) {
        Vector3q q = Vector3q::Zero();
        for (auto ev : mesh.edge_vertices(e)) {
            q += Q[ev];
        }
        q /= 2;
        Vertex v = mesh.split_edge(e);
        Q[v] = q;
        vside[v] = 0;
        for (auto ve : mesh.vertex_edges(v)) {
            eside[ve] = 0;
        }
        for (auto vf : mesh.vertex_faces(v)) {
            fside[vf] = 0;
        }
    }
    std::vector<Face> faces;
    for (auto f : mesh.faces()) {
        if (f_forbidden[f]) {
            continue;
        }
        bool incident0 = false;
        bool incident1 = false;
        for (auto fv : mesh.face_vertices(f)) {
            if (vside[fv] == -1) {
                incident0 = true;
            }
            if (vside[fv] == 1) {
                incident1 = true;
            }
        }
        for (auto fe : mesh.face_edges(f)) {
            if (eside[fe] == -1) {
                incident0 = true;
            }
            if (eside[fe] == 1) {
                incident1 = true;
            }
        }
        if (incident0 && incident1) {
            faces.push_back(f);
        }
    }
    for (auto f : faces) {
        Eigen::Vector3d p = Eigen::Vector3d::Zero();
        Vector3q q = Vector3q::Zero();
        for (auto fv : mesh.face_vertices(f)) {
            p += mesh.position(fv);
            q += Q[fv];
        }
        p /= 3;
        q /= 3;
        Vertex v = mesh.split_face(f, p);
        Q[v] = q;
        vside[v] = 0;
        for (auto ve : mesh.vertex_edges(v)) {
            eside[ve] = 0;
        }
        for (auto vf : mesh.vertex_faces(v)) {
            fside[vf] = 0;
        }
    }
    std::vector<Cell> cells;
    for (auto c : mesh.cells()) {
        bool incident0 = false;
        bool incident1 = false;
        for (auto cv : mesh.tet_vertices(c)) {
            if (vside[cv] == -1) {
                incident0 = true;
            }
            if (vside[cv] == 1) {
                incident1 = true;
            }
        }
        for (auto ce : mesh.cell_edges(c)) {
            if (eside[ce] == -1) {
                incident0 = true;
            }
            if (eside[ce] == 1) {
                incident1 = true;
            }
        }
        for (auto cf : mesh.cell_faces(c)) {
            if (fside[cf] == -1) {
                incident0 = true;
            }
            if (fside[cf] == 1) {
                incident1 = true;
            }
        }
        if (incident0 && incident1) {
            cells.push_back(c);
        }
    }
    for (auto c : cells) {
        Eigen::Vector3d p = Eigen::Vector3d::Zero();
        Vector3q q = Vector3q::Zero();
        for (auto cv : mesh.tet_vertices(c)) {
            p += mesh.position(cv);
            q += Q[cv];
        }
        p /= 4;
        q /= 4;
        Vertex v = mesh.split_tet(c, p);
        Q[v] = q;
    }
}

int euler_characteristic(Mesh& mesh, std::set<Halfface>& surface) {
    int n_vertices = 0;
    auto vertex_counted = mesh.property<Vertex, bool>("", false);
    int n_edges = 0;
    auto edge_counted = mesh.property<Edge, bool>("", false);
    for (auto h : surface) {
        for (auto hv : mesh.halfface_vertices(h)) {
            if (!vertex_counted[hv]) {
                n_vertices++;
                vertex_counted[hv] = true;
            }
        }
        for (auto he : mesh.halfface_edges(h)) {
            if (!edge_counted[he]) {
                n_edges++;
                edge_counted[he] = true;
            }
        }
    }
    return n_vertices - n_edges + (int) surface.size();
}

void laplace_surface(Mesh& mesh, std::vector<std::set<Halfface>>& sides, std::set<Halfface>& surface) {
    auto f_forbidden = mesh.property<Face, bool>("f_forbidden");
    auto e_forbidden = mesh.property<Edge, bool>("e_forbidden");
    auto v_forbidden = mesh.property<Vertex, bool>("v_forbidden");
    double invalid = 2;
    double l = -1;
    double u = 1;
    double mean = (l + u) / 2;
    auto value = mesh.property<Cell, double>("", invalid);
    std::queue<Cell> queue;
    auto enqueued = mesh.property<Cell, bool>("", false);
    for (auto h : sides[0]) {
        Cell hc = mesh.incident_cell(mesh.opposite_halfface_handle(h));
        queue.push(hc);
        enqueued[hc] = true;
    }
    for (auto h : sides[1]) {
        Cell hc = mesh.incident_cell(mesh.opposite_halfface_handle(h));
        queue.push(hc);
        enqueued[hc] = true;
    }
    while (!queue.empty()) {
        Cell c = queue.front();
        queue.pop();
        value[c] = mean;
        for (auto cc : mesh.cell_cells(c)) {
            if (!enqueued[cc]) {
                queue.push(cc);
                enqueued[cc] = true;
            }
        }
    }
    for (auto h : sides[0]) {
        Cell hc = mesh.incident_cell(mesh.opposite_halfface_handle(h));
        ASSERT(value[hc] != invalid);
        value[hc] = l;
        for (auto he : mesh.halfface_edges(h)) {
            if (e_forbidden[he]) {
                for (auto hec : mesh.edge_cells(he)) {
                    if (value[hec] != invalid) {
                        value[hec] = l;
                    }
                }
            }
        }
        for (auto hv : mesh.halfface_vertices(h)) {
            if (v_forbidden[hv]) {
                for (auto hvc : mesh.vertex_cells(hv)) {
                    if (value[hvc] != invalid) {
                        value[hvc] = l;
                    }
                }
            }
        }
    }
    for (auto h : sides[1]) {
        Cell hc = mesh.incident_cell(mesh.opposite_halfface_handle(h));
        ASSERT(value[hc] != invalid && value[hc] != l);
        value[hc] = u;
        for (auto he : mesh.halfface_edges(h)) {
            if (e_forbidden[he]) {
                for (auto hec : mesh.edge_cells(he)) {
                    if (value[hec] != invalid) {
                        ASSERT(value[hec] != l);
                        value[hec] = u;
                    }
                }
            }
        }
        for (auto hv : mesh.halfface_vertices(h)) {
            if (v_forbidden[hv]) {
                for (auto hvc : mesh.vertex_cells(hv)) {
                    if (value[hvc] != invalid) {
                        ASSERT(value[hvc] != l);
                        value[hvc] = u;
                    }
                }
            }
        }
    }
    auto idx = mesh.property<Cell, int>();
    int n_unknown = 0;
    for (auto c : mesh.cells()) {
        if (value[c] == mean) {
            idx[c] = n_unknown;
            n_unknown++;
        }
    }
    if (n_unknown > 0) {
        std::vector<Eigen::Triplet<double>> triplets;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(n_unknown);
        for (auto c : mesh.cells()) {
            if (value[c] == mean) {
                int n_neighbors = 0;
                for (auto cc : mesh.cell_cells(c)) {
                    if (value[cc] == mean) {
                        triplets.push_back(Eigen::Triplet<double>(idx[c], idx[cc], -1));
                    } else {
                        b(idx[c]) += value[cc];
                    }
                    n_neighbors++;
                }
                triplets.push_back(Eigen::Triplet<double>(idx[c], idx[c], n_neighbors));
            }
        }
        Eigen::SparseMatrix<double> a(n_unknown, n_unknown);
        a.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
        solver.compute(a);
        ASSERT(solver.info() == Eigen::Success);
        Eigen::VectorXd x = solver.solve(b);
        for (auto c : mesh.cells()) {
            if (value[c] == mean) {
                value[c] = x(idx[c]);
            }
        }
    }
    for (auto h : mesh.halffaces()) {
        if (!mesh.is_boundary(mesh.face_handle(h)) && value[mesh.incident_cell(h)] >= mean && value[mesh.incident_cell(mesh.opposite_halfface_handle(h))] < mean) {
            surface.insert(h);
        }
    }
    { // 3 -> 1 shifts
        std::queue<Cell> queue;
        auto enqueued = mesh.property<Cell, bool>("", false);
        for (auto h : surface) {
            Cell c = mesh.incident_cell(h);
            if (!enqueued[c]) {
                queue.push(c);
                enqueued[c] = true;
            }
            c = mesh.incident_cell(mesh.opposite_halfface_handle(h));
            if (!enqueued[c]) {
                queue.push(c);
                enqueued[c] = true;
            }
        }
        while (!queue.empty()) {
            Cell c = queue.front();
            queue.pop();
            Halfface h;
            int cnt = 0;
            for (auto ch : mesh.cell_halffaces(c)) {
                if (surface.count(ch) == 1) {
                    cnt++;
                } else if (surface.count(mesh.opposite_halfface_handle(ch)) == 1) {
                    cnt--;
                } else {
                    h = ch;
                }
            }
            if ((cnt == -3 || cnt == 3) && !mesh.is_boundary(mesh.opposite_halfface_handle(h))) {
                for (auto cf : mesh.cell_faces(c)) {
                    surface.erase(mesh.halfface_handle(cf, 0));
                    surface.erase(mesh.halfface_handle(cf, 1));
                }
                if (cnt == -3) {
                    surface.insert(h);
                } else {
                    surface.insert(mesh.opposite_halfface_handle(h));
                }
                Cell hc = mesh.incident_cell(mesh.opposite_halfface_handle(h));
                if (!enqueued[hc]) {
                    queue.push(hc);
                    enqueued[hc] = true;
                }
            }
        }
    }
    { // Naive non-manifold shifts
        std::queue<Edge> queue;
        auto enqueued = mesh.property<Edge, bool>("", false);
        for (auto h : surface) {
            for (auto he : mesh.halfface_edges(h)) {
                if (!enqueued[he]) {
                    queue.push(he);
                    enqueued[he] = true;
                }
            }
        }
        while (!queue.empty()) {
            Edge e = queue.front();
            queue.pop();
            int cnt = 0;
            for (auto eh : mesh.edge_halffaces(e)) {
                if (surface.count(eh) == 1) {
                    cnt++;
                }
            }
            if (cnt > 2) {
                double avg_value = 0;
                for (auto ec : mesh.edge_cells(e)) {
                    avg_value += value[ec];
                }
                avg_value /= mesh.valence(e);
                for (auto ec : mesh.edge_cells(e)) {
                    value[ec] = avg_value;
                    for (auto ecf : mesh.cell_faces(ec)) {
                        surface.erase(mesh.halfface_handle(ecf, 0));
                        surface.erase(mesh.halfface_handle(ecf, 1));
                    }
                }
                for (auto ec : mesh.edge_cells(e)) {
                    for (auto ecf : mesh.cell_faces(ec)) {
                        for (int i = 0; i < 2; i++) {
                            Halfface h = mesh.halfface_handle(ecf, i);
                            if (!mesh.is_boundary(mesh.face_handle(h)) && value[mesh.incident_cell(h)] >= mean && value[mesh.incident_cell(mesh.opposite_halfface_handle(h))] < mean) {
                                surface.insert(h);
                            }
                        }
                    }
                }
            }
        }
    }
}

void shift(Mesh& mesh, std::set<Halfface>& surface) {
    auto f_forbidden = mesh.property<Face, bool>("f_forbidden");
    auto e_forbidden = mesh.property<Edge, bool>("e_forbidden");
    auto v_forbidden = mesh.property<Vertex, bool>("v_forbidden");
    auto f_surface = [&](Face f) {
        for (int i = 0; i < 2; i++) {
            if (surface.count(mesh.halfface_handle(f, i)) == 1) {
                return true;
            }
        }
        return false;
    };
    auto e_surface = [&](Edge e) {
        for (auto eh : mesh.edge_halffaces(e)) {
            if (surface.count(eh) == 1) {
                return true;
            }
        }
        return false;
    };
    auto v_surface = [&](Vertex v) {
        for (auto vh : mesh.vertex_halffaces(v)) {
            if (surface.count(vh) == 1) {
                return true;
            }
        }
        return false;
    };
    auto Q = mesh.property<Vertex, Vector3q>("Q");
    auto split_edge = [&](Edge e) {
        ASSERT(!e_forbidden[e] && !e_surface(e));
        Vector3q q = Vector3q::Zero();
        for (auto ev : mesh.edge_vertices(e)) {
            q += Q[ev];
        }
        q /= 2;
        Vertex v = mesh.split_edge(e);
        Q[v] = q;
    };
    auto split_face = [&](Face f) {
        ASSERT(!f_forbidden[f] && !f_surface(f));
        Eigen::Vector3d p = Eigen::Vector3d::Zero();
        Vector3q q = Vector3q::Zero();
        for (auto fv : mesh.face_vertices(f)) {
            p += mesh.position(fv);
            q += Q[fv];
        }
        p /= 3;
        q /= 3;
        Vertex v = mesh.split_face(f, p);
        Q[v] = q;
    };
    auto split_cell = [&](Cell c) {
        Eigen::Vector3d p = Eigen::Vector3d::Zero();
        Vector3q q = Vector3q::Zero();
        for (auto cv : mesh.tet_vertices(c)) {
            p += mesh.position(cv);
            q += Q[cv];
        }
        p /= 4;
        q /= 4;
        Vertex v = mesh.split_tet(c, p);
        Q[v] = q;
    };
    auto mark_to = [&](std::set<Cell>& block, std::set<Halfface>& seeds, Property<Face, bool>& f_to, Property<Edge, bool>& e_to, Property<Vertex, bool>& v_to) {
        for (auto c : block) {
            for (auto ch : mesh.cell_halffaces(c)) {
                Halfface h = mesh.opposite_halfface_handle(ch);
                if (mesh.is_boundary(h) || block.count(mesh.incident_cell(h)) == 0) {
                    Face f = mesh.face_handle(h);
                    f_to[f] = true;
                    for (auto fe : mesh.face_edges(f)) {
                        e_to[fe] = true;
                    }
                    for (auto fv : mesh.face_vertices(f)) {
                        v_to[fv] = true;
                    }
                }
            }
        }
        std::queue<Halfface> queue;
        auto enqueued = mesh.property<Halfface, bool>("", false);
        for (auto h : seeds) {
            queue.push(h);
            enqueued[h] = true;
        }
        while (!queue.empty()) {
            Halfface h = queue.front();
            queue.pop();
            Face f = mesh.face_handle(h);
            ASSERT(surface.count(h) == 1 && f_to[f]);
            int n_filled_edges = 0;
            for (auto fe : mesh.face_edges(f)) {
                if (!e_to[fe]) {
                    n_filled_edges++;
                }
            }
            int n_filled_vertices = 0;
            for (auto fv : mesh.face_vertices(f)) {
                if (!v_to[fv]) {
                    n_filled_vertices++;
                }
            }
            ASSERT((n_filled_edges >= 1 && n_filled_vertices >= 2) || seeds.count(h) == 1);
            if (n_filled_edges == 1 && n_filled_vertices == 3) {
                continue; // Keep floor manifold and therefore of disk topology
            }
            f_to[f] = false;
            for (auto fe : mesh.face_edges(f)) {
                e_to[fe] = false;
            }
            for (auto fv : mesh.face_vertices(f)) {
                v_to[fv] = false;
            }
            for (auto he : mesh.halfface_halfedges(h)) {
                for (auto heh : mesh.halfedge_halffaces(he)) {
                    Halfface hh = mesh.opposite_halfface_handle(heh);
                    if (!enqueued[hh] && f_to[mesh.face_handle(hh)] && surface.count(hh) == 1) {
                        queue.push(hh);
                        enqueued[hh] = true;
                    }
                }
            }
        }
    };
    auto get_dome = [&](Vertex v) {
        std::set<Cell> dome;
        std::queue<Cell> queue;
        auto enqueued = mesh.property<Cell, bool>("", false);
        for (auto vh : mesh.vertex_halffaces(v)) {
            if (surface.count(vh) == 1) {
                Cell c = mesh.incident_cell(mesh.opposite_halfface_handle(vh));
                queue.push(c);
                enqueued[c] = true;
            }
        }
        while (!queue.empty()) {
            Cell c = queue.front();
            queue.pop();
            dome.insert(c);
            for (auto ch : mesh.cell_halffaces(c)) {
                bool incident_to_v = false;
                for (auto chv : mesh.halfface_vertices(ch)) {
                    if (chv == v) {
                        incident_to_v = true;
                        break;
                    }
                }
                if (!incident_to_v) {
                    continue;
                }
                Halfface h = mesh.opposite_halfface_handle(ch);
                if (surface.count(h) == 1) {
                    continue;
                }
                ASSERT(!mesh.is_boundary(h));
                Cell cc = mesh.incident_cell(h);
                if (!enqueued[cc]) {
                    queue.push(cc);
                    enqueued[cc] = true;
                }
            }
        }
        return dome;
    };
    auto get_fan = [&](Edge e) {
        std::set<Cell> fan;
        Halfface h;
        for (auto eh : mesh.edge_halffaces(e)) {
            if (surface.count(eh) == 1) {
                h = mesh.opposite_halfface_handle(eh);
                break;
            }
        }
        while (surface.count(h) == 0) {
            ASSERT(!mesh.is_boundary(h));
            Cell c = mesh.incident_cell(h);
            fan.insert(c);
            h = mesh.opposite_halfface_handle(mesh.adjacent_halfface_in_cell(h, mesh.halfedge_handle(e, 0))); // Any halfedge can be used
        }
        return fan;
    };
    struct Elements {
        std::set<Edge> edges;
        std::set<Face> faces;
        std::set<Cell> cells;
    };
    auto dome_split_elements = [&](Vertex v) {
        std::set<Cell> dome = get_dome(v);
        std::set<Halfface> seeds;
        for (auto vh : mesh.vertex_halffaces(v)) {
            if (surface.count(vh) == 1) {
                seeds.insert(vh);
            }
        }
        auto f_to = mesh.property<Face, bool>("", false);
        auto e_to = mesh.property<Edge, bool>("", false);
        auto v_to = mesh.property<Vertex, bool>("", false);
        mark_to(dome, seeds, f_to, e_to, v_to);
        Elements split_elements;
        std::set<Face> handled_faces;
        std::set<Cell> handled_cells;
        // VV edges
        for (auto ve : mesh.vertex_edges(v)) {
            auto [v0, v1] = mesh.edge_vertices(ve);
            Vertex opp = v0 == v ? v1 : v0;
            if (!v_to[opp]) {
                continue;
            }
            ASSERT(!e_forbidden[ve] && !e_surface(ve));
            if (v_forbidden[opp] || v_surface(opp)) {
                split_elements.edges.insert(ve);
                for (auto vef : mesh.edge_faces(ve)) {
                    handled_faces.insert(vef);
                }
                for (auto vec : mesh.edge_cells(ve)) {
                    handled_cells.insert(vec);
                }
            }
        }
        // VE triangles
        for (auto vf : mesh.vertex_faces(v)) {
            if (handled_faces.count(vf) == 1) {
                continue;
            }
            Edge opp;
            for (auto vfe : mesh.face_edges(vf)) {
                auto [v0, v1] = mesh.edge_vertices(vfe);
                if (v0 != v && v1 != v) {
                    opp = vfe;
                    break;
                }
            }
            if (!e_to[opp]) {
                continue;
            }
            ASSERT(!f_forbidden[vf] && !f_surface(vf));
            if (e_forbidden[opp] || e_surface(opp)) {
                split_elements.faces.insert(vf);
                for (auto vfc : mesh.face_cells(vf)) {
                    if (mesh.is_valid(vfc)) {
                        handled_cells.insert(vfc);
                    }
                }
            }
        }
        // VF tetrahedra
        for (auto vc : mesh.vertex_cells(v)) {
            if (handled_cells.count(vc) == 1) {
                continue;
            }
            Face opp = mesh.face_handle(mesh.vertex_opposite_halfface(v, vc));
            if (!f_to[opp]) {
                continue;
            }
            ASSERT(dome.count(vc) == 1);
            if (f_forbidden[opp] || f_surface(opp)) {
                split_elements.cells.insert(vc);
            }
        }
        return split_elements;
    };
    auto shift = [&](std::set<Cell>& block) {
        for (auto c : block) {
            for (auto ch : mesh.cell_halffaces(c)) {
                Face f = mesh.face_handle(ch);
                Halfface h = mesh.opposite_halfface_handle(ch);
                if (mesh.is_boundary(h) || block.count(mesh.incident_cell(h)) == 0) {
                    if (surface.count(h) == 1) {
                        surface.erase(h);
                    } else {
                        ASSERT(!f_forbidden[f]);
                        for (auto fe : mesh.face_edges(f)) {
                            if (!e_surface(fe)) {
                                ASSERT(!e_forbidden[fe]);
                            }
                        }
                        for (auto fv : mesh.face_vertices(f)) {
                            if (!v_surface(fv)) {
                                ASSERT(!v_forbidden[fv]);
                            }
                        }
                        surface.insert(ch);
                    }
                }
            }
        }
    };

    // Free forbidden vertices
    std::set<Vertex> forbidden_surface_vertices;
    for (auto h : surface) {
        for (auto hv : mesh.halfface_vertices(h)) {
            if (v_forbidden[hv]) {
                forbidden_surface_vertices.insert(hv);
            }
        }
    }
    int factor = mesh.n_vertices();
    auto cmp = [factor](const std::pair<Vertex, Elements>& a, const std::pair<Vertex, Elements>& b) {
        return factor * (a.second.edges.size() + a.second.faces.size() + a.second.cells.size()) + a.first.idx() >
               factor * (b.second.edges.size() + b.second.faces.size() + b.second.cells.size()) + b.first.idx();
    };
    PriorityQueue<std::pair<Vertex, Elements>, decltype(cmp)> priority_queue(cmp);
    for (auto v : forbidden_surface_vertices) {
        priority_queue.push_or_update({v, dome_split_elements(v)}, v.idx());
    }
    while (!priority_queue.empty()) {
        auto [v, split_elements] = priority_queue.pop();
        if (!v_surface(v)) {
            continue;
        }
        // Split
        for (auto e : split_elements.edges) {
            split_edge(e);
        }
        for (auto f : split_elements.faces) {
            split_face(f);
        }
        for (auto c : split_elements.cells) {
            split_cell(c);
        }
        // Shift dome
        std::set<Cell> dome = get_dome(v);
        shift(dome);
        // Update
        std::set<Vertex> dome_vertices;
        for (auto c : dome) {
            for (auto cv : mesh.tet_vertices(c)) {
                dome_vertices.insert(cv);
            }
        }
        std::set<Vertex> dome_vertex_vertices;
        for (auto v : dome_vertices) {
            for (auto vv : mesh.vertex_vertices(v)) {
                dome_vertex_vertices.insert(vv);
            }
        }
        for (auto v : dome_vertex_vertices) {
            if (v_forbidden[v] && v_surface(v)) {
                priority_queue.push_or_update({v, dome_split_elements(v)}, v.idx());
            }
        }
    }

    // Free forbidden edges
    std::set<Edge> forbidden_surface_edges;
    for (auto h : surface) {
        for (auto he : mesh.halfface_edges(h)) {
            if (e_forbidden[he]) {
                forbidden_surface_edges.insert(he);
            }
        }
    }
    for (auto e : forbidden_surface_edges) {
        if (!e_surface(e)) {
            continue;
        }
        std::set<Cell> fan = get_fan(e);
        std::set<Halfface> seeds;
        for (auto eh : mesh.edge_halffaces(e)) {
            if (surface.count(eh) == 1) {
                seeds.insert(eh);
            }
        }
        auto f_to = mesh.property<Face, bool>("", false);
        auto e_to = mesh.property<Edge, bool>("", false);
        auto v_to = mesh.property<Vertex, bool>("", false);
        mark_to(fan, seeds, f_to, e_to, v_to);
        auto [v0, v1] = mesh.edge_vertices(e);
        // EV / EE triangles
        std::vector<Face> faces;
        for (auto ef : mesh.edge_faces(e)) {
            Vertex opp_v;
            for (auto efv : mesh.face_vertices(ef)) {
                if (efv != v0 && efv != v1) {
                    opp_v = efv;
                    break;
                }
            }
            std::vector<Edge> opp_e;
            for (auto efe : mesh.face_edges(ef)) {
                if (efe != e) {
                    opp_e.push_back(efe);
                }
            }
            ASSERT(opp_e.size() == 2);
            if (!v_to[opp_v] && !e_to[opp_e[0]] && !e_to[opp_e[1]]) {
                continue;
            }
            ASSERT(!f_forbidden[ef] && !f_surface(ef));
            if (v_forbidden[opp_v] || v_surface(opp_v) ||
                e_forbidden[opp_e[0]] || e_surface(opp_e[0]) ||
                e_forbidden[opp_e[1]] || e_surface(opp_e[1])) {
                faces.push_back(ef);
            }
        }
        for (auto f : faces) {
            split_face(f);
        }
        // EE / EF tetrahedra
        std::vector<Cell> cells;
        for (auto ec : mesh.edge_cells(e)) {
            Edge opp_e = mesh.edge_handle(mesh.halfedge_opposite_halfedge(mesh.halfedge_handle(e, 0), ec));
            std::vector<Face> opp_f = {
                mesh.face_handle(mesh.vertex_opposite_halfface(v0, ec)),
                mesh.face_handle(mesh.vertex_opposite_halfface(v1, ec))};
            if (!e_to[opp_e] && !f_to[opp_f[0]] && !f_to[opp_f[1]]) {
                continue;
            }
            ASSERT(fan.count(ec) == 1);
            if (e_forbidden[opp_e] || e_surface(opp_e) ||
                f_forbidden[opp_f[0]] || f_surface(opp_f[0]) ||
                f_forbidden[opp_f[1]] || f_surface(opp_f[1])) {
                cells.push_back(ec);
            }
        }
        for (auto c : cells) {
            split_cell(c);
        }
        // Shift fan
        fan = get_fan(e);
        shift(fan);
    }

    // Free forbidden faces
    std::set<Face> forbidden_surface_faces;
    for (auto h : surface) {
        Face f = mesh.face_handle(h);
        if (f_forbidden[f]) {
            forbidden_surface_faces.insert(f);
        }
    }
    for (auto f : forbidden_surface_faces) {
        if (!f_surface(f)) {
            continue;
        }
        Halfface h;
        for (auto fh : mesh.face_halffaces(f)) {
            if (surface.count(fh) == 1) {
                h = fh;
                break;
            }
        }
        Cell c = mesh.incident_cell(mesh.opposite_halfface_handle(h));
        std::set<Cell> tet;
        tet.insert(c);
        std::set<Halfface> seed;
        seed.insert(h);
        auto f_to = mesh.property<Face, bool>("", false);
        auto e_to = mesh.property<Edge, bool>("", false);
        auto v_to = mesh.property<Vertex, bool>("", false);
        mark_to(tet, seed, f_to, e_to, v_to);
        // FV / FE / FF tetrahedron
        bool split = false;
        for (auto cv : mesh.tet_vertices(c)) {
            if (v_to[cv] && (v_forbidden[cv] || v_surface(cv))) {
                split = true;
                break;
            }
        }
        if (!split) {
            for (auto ce : mesh.cell_edges(c)) {
                if (e_to[ce] && (e_forbidden[ce] || e_surface(ce))) {
                    split = true;
                    break;
                }
            }
        }
        if (!split) {
            for (auto cf : mesh.cell_faces(c)) {
                if (f_to[cf] && f_forbidden[cf]) {
                    split = true;
                    break;
                }
            }
        }
        if (split) {
            split_cell(c);
            std::cout << "Split FV / FE / FF tetrahedron" << std::endl;
        }
        // Shift tet
        tet.clear();
        c = mesh.incident_cell(mesh.opposite_halfface_handle(h));
        tet.insert(c);
        shift(tet);
    }

    // 3 -> 1 shifts
    std::queue<Cell> queue;
    auto enqueued = mesh.property<Cell, bool>("", false);
    for (auto h : surface) {
        Cell c = mesh.incident_cell(mesh.opposite_halfface_handle(h));
        if (!enqueued[c]) {
            queue.push(c);
            enqueued[c] = true;
        }
    }
    while (!queue.empty()) {
        Cell c = queue.front();
        queue.pop();
        enqueued[c] = false;
        int n_surface_faces = 0;
        Halfface h;
        for (auto ch : mesh.cell_halffaces(c)) {
            if (surface.count(mesh.opposite_halfface_handle(ch)) == 1) {
                n_surface_faces++;
            } else {
                h = ch;
            }
        }
        if (n_surface_faces < 3 || f_forbidden[mesh.face_handle(h)]) {
            continue;
        }
        std::set<Cell> tet;
        tet.insert(c);
        shift(tet);
        Cell cc = mesh.incident_cell(mesh.opposite_halfface_handle(h));
        if (!enqueued[cc]) {
            queue.push(cc);
            enqueued[cc] = true;
        }
    }
}

bool to_mesh(Mesh& mesh, std::set<Halfface>& surface, TriMesh& result) {
    auto f_forbidden = mesh.property<Face, bool>("f_forbidden");
    auto e_forbidden = mesh.property<Edge, bool>("e_forbidden");
    auto v_forbidden = mesh.property<Vertex, bool>("v_forbidden");
    std::stringstream tmp_buffer;
    std::streambuf* cerr_buffer = std::cerr.rdbuf(tmp_buffer.rdbuf());
    auto M2m = mesh.property<Vertex, SVertex>("M2m");
    for (auto v : mesh.vertices()) {
        M2m[v] = SVertex(-1);
    }
    auto m2M = result.property<SVertex, Vertex>("m2M", Vertex(-1));
    for (auto h : surface) {
        std::vector<SVertex> vertices;
        for (auto hv : mesh.halfface_vertices(h)) {
            if (!result.is_valid_handle(M2m[hv])) {
                SVertex v = result.add_vertex(mesh.position(hv));
                M2m[hv] = v;
                m2M[v] = hv;
            }
            vertices.push_back(M2m[hv]);
        }
        if (!result.is_valid_handle(result.add_face(vertices))) {
            for (auto v : vertices) {
                std::cerr.rdbuf(cerr_buffer);
                return false;
            }
        }
    }
    std::cerr.rdbuf(cerr_buffer);
    int euler_characteristic = (int) result.n_vertices() - (int) result.n_edges() + (int) result.n_faces();
    if (euler_characteristic != 1) {
        return false;
    }
    for (auto f : result.faces()) {
        std::vector<Vertex> vertices;
        for (auto fv : result.face_vertices(f)) {
            vertices.push_back(m2M[fv]);
        }
        Face mf = mesh.face_handle(mesh.find_halfface(vertices));
        if (f_forbidden[mf]) {
            return false;
        }
    }
    for (auto e : result.edges()) {
        SHalfedge h = result.halfedge_handle(e, 0);
        Edge me = mesh.edge_handle(mesh.find_halfedge(m2M[result.from_vertex_handle(h)], m2M[result.to_vertex_handle(h)]));
        if (e_forbidden[me]) {
            return false;
        }
    }
    for (auto v : result.vertices()) {
        if (v_forbidden[m2M[v]]) {
            return false;
        }
    }
    return true;
}

bool cut(Mesh& mesh, std::vector<Vertex> loop, TriMesh& result) {
    forbid(mesh, loop);
    std::vector<std::set<Halfface>> sides;
    fill_sides(mesh, loop, sides);
    pre_split(mesh, sides);
    int chi0 = euler_characteristic(mesh, sides[0]);
    int chi1 = euler_characteristic(mesh, sides[1]);
    if (chi0 == 1 && chi1 == 1) {
        std::set<Halfface> surface;
        laplace_surface(mesh, sides, surface);
        if (to_mesh(mesh, surface, result)) {
            return true;
        }
    }
    result = TriMesh();
    std::set<Halfface>& surface = sides[0];
    int chi = chi0;
    bool flip = true;
    if (chi != 1 || (sides[1].size() < surface.size() && chi1 == 1)) {
        surface = sides[1];
        chi = chi1;
        flip = false;
    }
    if (chi != 1) {
        return false;
    }
    shift(mesh, surface);
    if (flip) {
        std::set<Halfface> flipped_surface;
        for (auto h : surface) {
            flipped_surface.insert(mesh.opposite_halfface_handle(h));
        }
        surface = flipped_surface;
    }
    ASSERT(to_mesh(mesh, surface, result));
    return true;
}
