#include "sdm.h"

std::vector<Vertex> link(Mesh& mesh, TriMesh& mcut, SProperty<SVertex, Vertex>& mc2m, TriMesh& cut, SProperty<SVertex, BarycentricVector>& c2m, std::function<Vertex(Face)> split_face) {
    auto Q = mesh.property<Vertex, Vector3q>("Q");
    auto cmp = mesh.property<Cell, int>("cmp");
    auto bdr = mesh.property<Halfface, bool>("bdr", false);

    auto is_cut_halfface = mesh.property<Halfface, bool>("", false);
    for (auto mcf : mcut.faces()) {
        std::vector<Vertex> vertices;
        for (auto mcfv : mcut.fv_ccw_range(mcf)) {
            vertices.push_back(mc2m[mcfv]);
        }
        is_cut_halfface[mesh.find_halfface({vertices})] = true;
    }
    auto is_cut_face = [&](Face f) {
        Halfface h = mesh.halfface_handle(f, 0);
        return is_cut_halfface[h] || is_cut_halfface[mesh.opposite_halfface_handle(h)];
    };

    // Recursively splits a tet until there are no two cutout edges remaining and saves the children
    auto children = mesh.property<Cell, std::vector<Cell>>();
    std::function<void(Cell)> split = [&](Cell c) {
        if (children[c].size() > 0) {
            for (auto child : children[c]) {
                split(child);
            }
            return;
        }
        int n_1_cutout_edge_faces = 0;
        for (auto ch : mesh.cell_halffaces(c)) {
            if (!is_cut_face(mesh.face_handle(ch))) {
                int n_cutout_edges = 0;
                for (auto che : mesh.halfface_edges(ch)) {
                    bool cutout_edge = false;
                    for (auto cheh : mesh.edge_halffaces(che)) {
                        if (is_cut_face(mesh.face_handle(cheh))) {
                            cutout_edge = true;
                            break;
                        }
                    }
                    if (cutout_edge) {
                        n_cutout_edges++;
                    }
                }
                if (n_cutout_edges > 1) {
                    // Split non-cutout face with > 1 cutout edges
                    Cell c2 = mesh.incident_cell(mesh.opposite_halfface_handle(ch));
                    Vertex o = mesh.halfface_opposite_vertex(ch);
                    Vertex v = split_face(mesh.face_handle(ch)); // Cell properties are copied
                    for (auto vc : mesh.vertex_cells(v)) {
                        bool has_o = false;
                        for (auto vcv : mesh.tet_vertices(vc)) {
                            if (vcv == o) {
                                has_o = true;
                                break;
                            }
                        }
                        if (has_o) {
                            children[c].push_back(vc);
                        } else if (mesh.is_valid(c2)) {
                            children[c2].push_back(vc);
                        }
                    }
                    // Split children recursively
                    for (auto child : children[c]) {
                        split(child);
                    }
                    return;
                }
                if (n_cutout_edges == 1) {
                    n_1_cutout_edge_faces++;
                }
            }
        }
        if (n_1_cutout_edge_faces == 4) {
            // Split tet with two opposite cutout edges (Case 19)
            Eigen::Vector3d p = Eigen::Vector3d::Zero();
            Vector3q q = Vector3q::Zero();
            for (auto cv : mesh.tet_vertices(c)) {
                p += mesh.position(cv);
                q += Q[cv];
            }
            p /= 4;
            q /= 4;
            Vertex v = mesh.split_tet(c, p); // Cell properties are not copied
            Q[v] = q;
            for (auto vc : mesh.vertex_cells(v)) {
                cmp[vc] = cmp[c];
                children[c].push_back(vc);
            }
        }
    };
    std::set<Cell> link_cells;
    for (auto mcv : mcut.vertices()) {
        Vertex v = mc2m[mcv];
        for (auto vc : mesh.vertex_cells(v)) {
            link_cells.insert(vc);
        }
    }
    for (auto c : link_cells) {
        split(c);
    }

    // Add new vertices
    std::vector<Vertex> new_vertices;
    std::map<Vertex, Vertex> m2m_new;
    auto N2N = mesh.property<Vertex, Vertex>("N2N", Vertex(-1));
    auto c2m_new = cut.property<SVertex, Vertex>("", Vertex(-1));
    for (auto cv : cut.vertices()) {
        BarycentricVector b = c2m[cv];
        Eigen::Vector3d p = Eigen::Vector3d::Zero();
        Vector3q q = Vector3q::Zero();
        for (int i = 0; i < b.vertices.size(); i++) {
            p += b.coordinates[i].get_d() * mesh.position(b.vertices[i]);
            q += b.coordinates[i] * Q[b.vertices[i]];
        }
        Vertex v = mesh.add_vertex(p);
        Q[v] = q;
        new_vertices.push_back(v);
        if (b.vertices.size() == 1) {
            m2m_new[b.vertices[0]] = v;
            N2N[b.vertices[0]] = v;
        }
        c2m_new[cv] = v;
    }

    // Bouquet
    for (auto cf : cut.faces()) {
        std::vector<Vertex> vertices;
        std::vector<Vertex> super_vertices;
        for (auto cfv : cut.fv_ccw_range(cf)) {
            vertices.push_back(c2m_new[cfv]);
            for (auto v : c2m[cfv].vertices) {
                if (std::find(super_vertices.begin(), super_vertices.end(), v) == super_vertices.end()) {
                    super_vertices.push_back(v);
                }
            }
        }
        Halfface h = mesh.find_halfface(super_vertices);
        ASSERT(mesh.is_valid(h));
        if (!is_cut_halfface[h]) {
            h = mesh.opposite_halfface_handle(h);
        }
        Vertex v = mesh.halfface_opposite_vertex(h);
        if (m2m_new.count(v) == 1) {
            v = m2m_new[v];
        }
        vertices.push_back(v);
        Cell c = mesh.add_cell(vertices, true);
        cmp[c] = cmp[mesh.incident_cell(h)];
        h = mesh.opposite_halfface_handle(h);
        v = mesh.halfface_opposite_vertex(h);
        if (m2m_new.count(v) == 1) {
            v = m2m_new[v];
        }
        vertices[3] = vertices[2];
        vertices[2] = v;
        c = mesh.add_cell(vertices, true);
        cmp[c] = cmp[mesh.incident_cell(h)];
        bdr[mesh.vertex_opposite_halfface(v, c)] = true;
    }

    // Fan
    for (auto ce : cut.edges()) {
        std::vector<Vertex> vertices;
        std::vector<bool> on_super;
        std::vector<Vertex> super_vertices;
        std::vector<mpq_class> alphas;
        for (auto cev : cut.edge_vertices(ce)) {
            vertices.push_back(c2m_new[cev]);
            BarycentricVector b = c2m[cev];
            on_super.push_back(b.vertices.size() == 1);
            mpq_class alpha = 0;
            for (int i = 0; i < b.vertices.size(); i++) {
                Vertex v = b.vertices[i];
                if (std::find(super_vertices.begin(), super_vertices.end(), v) == super_vertices.end()) {
                    super_vertices.push_back(v);
                }
                if (v == super_vertices[0]) {
                    alpha = b.coordinates[i];
                }
            }
            alphas.push_back(alpha);
        }
        if (super_vertices.size() != 2) {
            continue;
        }
        if (alphas[0] < alphas[1]) {
            Vertex tmp = super_vertices[0];
            super_vertices[0] = super_vertices[1];
            super_vertices[1] = tmp;
        }
        Halfedge e = mesh.find_halfedge(super_vertices[0], super_vertices[1]);
        ASSERT(mesh.is_valid(e));
        for (int i = 0; i < 2; i++) {
            vertices.push_back(Vertex(-1));
            super_vertices.push_back(Vertex(-1));
        }
        for (auto ec : mesh.halfedge_cells(e)) {
            bool face_connected = false;
            for (auto ecf : mesh.cell_faces(ec)) {
                if (is_cut_face(ecf)) {
                    face_connected = true;
                    break;
                }
            }
            if (face_connected) {
                continue;
            }
            Halfedge e2 = mesh.halfedge_opposite_halfedge(e, ec);
            Vertex v0 = mesh.from_vertex_handle(e2);
            super_vertices[2] = v0;
            if (m2m_new.count(v0) == 1) {
                v0 = m2m_new[v0];
            }
            vertices[2] = v0;
            on_super.push_back(true);
            Vertex v1 = mesh.to_vertex_handle(e2);
            super_vertices[3] = v1;
            if (m2m_new.count(v1) == 1) {
                v1 = m2m_new[v1];
            }
            vertices[3] = v1;
            on_super.push_back(true);
            Cell c = mesh.add_cell(vertices, true);
            cmp[c] = cmp[ec];
            for (int i0 = 0; i0 < 4; i0++) {
                for (int i1 = i0 + 1; i1 < 4; i1++) {
                    for (int i2 = i1 + 1; i2 < 4; i2++) {
                        for (int i3 = 0; i3 < 2; i3++) {
                            int a = i0;
                            int b = i3 == 0 ? i1 : i2;
                            int c = i3 == 0 ? i2 : i1;
                            if (on_super[a] && on_super[b] && on_super[c]) {
                                bdr[mesh.find_halfface({vertices[a], vertices[b], vertices[c]})] = bdr[mesh.find_halfface({super_vertices[a], super_vertices[b], super_vertices[c]})];
                            }
                        }
                    }
                }
            }
        }
    }

    // Vertex tetrahedra
    std::vector<std::vector<Vertex>> tets;
    std::vector<int> cmps;
    std::vector<std::vector<bool>> bdrs;
    for (auto cv : cut.vertices()) {
        std::vector<Vertex> super_vertices = c2m[cv].vertices;
        if (super_vertices.size() != 1) {
            continue;
        }
        Vertex v = super_vertices[0];
        for (auto vc : mesh.vertex_cells(v)) {
            bool edge_connected = false;
            for (auto vce : mesh.cell_edges(vc)) {
                for (auto vcef : mesh.edge_faces(vce)) {
                    if (is_cut_face(vcef)) {
                        edge_connected = true;
                        break;
                    }
                }
                if (edge_connected) {
                    break;
                }
            }
            if (edge_connected) {
                continue;
            }
            Halfface h = mesh.vertex_opposite_halfface(v, vc);
            std::vector<Vertex> vertices;
            super_vertices.clear();
            for (auto hv : mesh.halfface_vertices(h)) {
                super_vertices.push_back(hv);
                if (m2m_new.count(hv) == 1) {
                    hv = m2m_new[hv];
                }
                vertices.push_back(hv);
            }
            vertices.push_back(c2m_new[cv]);
            super_vertices.push_back(v);
            tets.push_back(vertices);
            cmps.push_back(cmp[vc]);
            std::vector<bool> tet_bdrs;
            for (int i0 = 0; i0 < 4; i0++) {
                for (int i1 = i0 + 1; i1 < 4; i1++) {
                    for (int i2 = i1 + 1; i2 < 4; i2++) {
                        for (int i3 = 0; i3 < 2; i3++) {
                            int a = i0;
                            int b = i3 == 0 ? i1 : i2;
                            int c = i3 == 0 ? i2 : i1;
                            tet_bdrs.push_back(bdr[mesh.find_halfface({super_vertices[a], super_vertices[b], super_vertices[c]})]);
                        }
                    }
                }
            }
            bdrs.push_back(tet_bdrs);
        }
    }
    for (auto mcv : mcut.vertices()) {
        mesh.delete_vertex(mc2m[mcv]);
    }
    for (int i = 0; i < tets.size(); i++) {
        Cell c = mesh.add_cell(tets[i], true);
        if (mesh.is_valid(c)) {
            cmp[c] = cmps[i];
            int j = 0;
            for (int i0 = 0; i0 < 4; i0++) {
                for (int i1 = i0 + 1; i1 < 4; i1++) {
                    for (int i2 = i1 + 1; i2 < 4; i2++) {
                        for (int i3 = 0; i3 < 2; i3++) {
                            int a = i0;
                            int b = i3 == 0 ? i1 : i2;
                            int c = i3 == 0 ? i2 : i1;
                            bdr[mesh.find_halfface({tets[i][a], tets[i][b], tets[i][c]})] = bdrs[i][j];
                            j++;
                        }
                    }
                }
            }
        }
    }

    return new_vertices;
}

std::vector<std::pair<Mesh, Mesh>> sdm(Mesh& M, Mesh& N) {
    for (auto c : M.cells()) {
        if (M.degenerate_or_inverted(c)) {
            std::cout << "Degenerate or inverted tetrahedron in M" << std::endl;
            return {};
        }
    }
    for (auto c : N.cells()) {
        if (N.degenerate_or_inverted(c)) {
            std::cout << "Degenerate or inverted tetrahedron in N" << std::endl;
            return {};
        }
    }
    auto M2N = M.property<Vertex, Vertex>("M2N", Vertex(-1));
    auto N2M = N.property<Vertex, Vertex>("N2M", Vertex(-1));
    auto hit = N.property<Vertex, bool>("", false);
    for (auto Mv : M.vertices()) {
        if (!M.is_boundary(Mv)) {
            continue;
        }
        Vertex Nv(Mv.idx());
        M2N[Mv] = Nv;
        N2M[Nv] = Mv;
        if (!N.is_valid(Nv) || !N.is_boundary(Nv) || hit[Nv]) {
            std::cout << "M and N not compatible" << std::endl;
            return {};
        }
        hit[Nv] = true;
    }
    if (M.boundary_vertices().size() != N.boundary_vertices().size()) {
        std::cout << "M and N not compatible" << std::endl;
        return {};
    }
    for (auto Mh : M.halffaces()) {
        if (!M.is_boundary(Mh)) {
            continue;
        }
        std::vector<Vertex> pvertices;
        for (auto Mhv : M.halfface_vertices(Mh)) {
            pvertices.push_back(M2N[Mhv]);
        }
        Halfface Nh = N.find_halfface(pvertices);
        if (!N.is_valid(Nh) || !N.is_boundary(Nh)) {
            std::cout << "M and N not compatible" << std::endl;
            return {};
        }
    }
    auto MQ = M.property<Vertex, Vector3q>("Q");
    for (auto v : M.vertices()) {
        MQ[v] = M.position(v).cast<mpq_class>();
    }
    auto NQ = N.property<Vertex, Vector3q>("Q");
    for (auto v : N.vertices()) {
        NQ[v] = N.position(v).cast<mpq_class>();
    }

    std::cout << "Decompose" << std::endl;
    std::vector<Vector3q> centers = decompose(N);
    auto cmp = N.property<Cell, int>("cmp");
    int n_cuts = 0;
    for (auto c : N.cells()) {
        int index = cmp[c];
        ASSERT(index >= 0);
        if (n_cuts < index) {
            n_cuts = index;
        }
    }
    std::cout << n_cuts + 1 << " components" << std::endl;
    std::vector<std::pair<Mesh, Mesh>> components(n_cuts + 1);

    std::queue<int> indices;
    for (int i = 0; i <= n_cuts; i++) {
        indices.push(i);
    }
    int n_skips = 0;
    int cnt = 0;
    int i = 0;
    while (!indices.empty()) {
        int index = indices.front();
        indices.pop();
        auto bdr = M.property<Halfface, bool>("bdr");
        for (auto h : M.halffaces()) {
            bdr[h] = false;
        }
        TriMesh n_union;
        auto N2u = N.property<Vertex, SVertex>("", SVertex(-1));
        auto u2N = n_union.property<SVertex, Vertex>("", Vertex(-1));
        for (auto c : N.cells()) {
            if (cmp[c] != index) {
                continue;
            }
            for (auto ch : N.cell_halffaces(c)) {
                Halfface h = N.opposite_halfface_handle(ch);
                if (N.is_boundary(h)) {
                    std::vector<Vertex> vertices;
                    for (auto chv : N.halfface_vertices(ch)) {
                        vertices.push_back(N2M[chv]);
                    }
                    bdr[M.find_halfface(vertices)] = true;
                } else if (cmp[N.incident_cell(h)] != index) {
                    std::vector<SVertex> vertices;
                    for (auto hv : N.halfface_vertices(h)) {
                        if (!n_union.is_valid_handle(N2u[hv])) {
                            SVertex v = n_union.add_vertex(N.position(hv));
                            N2u[hv] = v;
                            u2N[v] = hv;
                        }
                        vertices.push_back(N2u[hv]);
                    }
                    SFace f = n_union.add_face(vertices);
                    ASSERT(n_union.is_valid_handle(f));
                }
            }
        }

        bool disk_topology = true;
        std::vector<TriMesh> ns;
        auto handled = n_union.property<SFace, bool>("", false);
        int j = 0;
        for (auto face : n_union.faces()) {
            if (handled[face]) {
                continue;
            }
            TriMesh n;
            std::stringstream ss;
            ss << "N2n" << j;
            auto N2n = N.property<Vertex, SVertex>(ss.str());
            for (auto v : N.vertices()) {
                N2n[v] = SVertex(-1);
            }
            auto n2N = n.property<SVertex, Vertex>("n2N", Vertex(-1));
            std::queue<SFace> queue;
            queue.push(face);
            handled[queue.front()] = true;
            while (!queue.empty()) {
                SFace f = queue.front();
                queue.pop();
                std::vector<SVertex> vertices;
                for (auto fv : n_union.fv_ccw_range(f)) {
                    if (!n.is_valid_handle(N2n[u2N[fv]])) {
                        SVertex v = n.add_vertex(n_union.position(fv));
                        N2n[u2N[fv]] = v;
                        n2N[v] = u2N[fv];
                    }
                    vertices.push_back(N2n[u2N[fv]]);
                }
                n.add_face(vertices);
                for (auto ff : n_union.ff_range(f)) {
                    if (!handled[ff]) {
                        queue.push(ff);
                        handled[ff] = true;
                    }
                }
            }
            if ((int) n.n_vertices() - (int) n.n_edges() + (int) n.n_faces() != 1) {
                disk_topology = false;
                break;
            }
            for (auto v : n.vertices()) {
                if (!n.is_boundary(v) && N.is_boundary(n2N[v])) {
                    disk_topology = false;
                    break;
                }
            }
            if (!disk_topology) {
                break;
            }
            for (auto e : n.edges()) {
                if (n.is_boundary(e)) {
                    continue;
                }
                SHalfedge h = n.halfedge_handle(e, 0);
                Vertex u = n2N[n.from_vertex_handle(h)];
                Vertex v = n2N[n.to_vertex_handle(h)];
                if (N.is_boundary(N.find_halfedge(u, v))) {
                    disk_topology = false;
                    break;
                }
            }
            if (!disk_topology) {
                break;
            }
            ns.push_back(n);
            j++;
        }

        if (!disk_topology) {
            indices.push(index);
            n_skips++;
            if (n_skips > n_cuts) {
                std::cout << "No disk topology order was found" << std::endl;
                return {};
            }
            continue;
        }

        for (auto n : ns) {
            auto n2N = n.property<SVertex, Vertex>("n2N");
            for (auto f : n.faces()) {
                std::vector<Vertex> vertices;
                for (auto fv : n.face_vertices(f)) {
                    vertices.push_back(n2N[fv]);
                }
            }
        }

        for (j = 0; j < ns.size(); j++) {
            cnt++;
            std::cout << "Cut " << cnt << " / " << n_cuts << std::endl;

            TriMesh n = ns[j];
            std::stringstream ss;
            ss << "N2n" << j;
            auto N2n = N.property<Vertex, SVertex>(ss.str());
            auto n2N = n.property<SVertex, Vertex>("n2N");
            for (auto v : n.vertices()) {
                auto N2N = N.property<Vertex, Vertex>("N2N", Vertex(-1));
                while (N.is_deleted(n2N[v])) {
                    n2N[v] = N2N[n2N[v]];
                    N2n[n2N[v]] = v;
                }
                ASSERT(n.is_boundary(v) == N.is_boundary(n2N[v]));
            }
            SHalfedge s;
            for (auto h : n.halfedges()) {
                if (h.is_boundary()) {
                    s = h;
                    break;
                }
            }
            std::vector<Vertex> loop;
            SHalfedge h = s;
            do {
                loop.push_back(N2M[n2N[n.to_vertex_handle(h)]]);
                h = n.next_halfedge_handle(h);
            } while (h != s);

            TriMesh m;
            if (!cut(M, loop, m)) {
                std::cout << "No ball topology" << std::endl;
                return {};
            }

            std::cout << "Align " << m.n_faces() << " x " << n.n_faces() << std::endl;
            auto m2M = m.property<SVertex, Vertex>("m2M", Vertex(-1));
            auto m2n = m.property<SVertex, SVertex>("m2n", SVertex(-1));
            auto n2m = n.property<SVertex, SVertex>("n2m", SVertex(-1));
            for (auto mv : m.vertices()) {
                if (!m.is_boundary(mv)) {
                    continue;
                }
                SVertex nv = N2n[M2N[m2M[mv]]];
                m2n[mv] = nv;
                n2m[nv] = mv;
            }
            TriMesh x;
            align(m, n, x);

            auto x2M = x.property<SVertex, BarycentricVector>("x2M");
            auto x2N = x.property<SVertex, BarycentricVector>("x2N");
            auto centroid = [&](Mesh& mesh, Face f) {
                auto Q = mesh.property<Vertex, Vector3q>("Q");
                Eigen::Vector3d p = Eigen::Vector3d::Zero();
                Vector3q q = Vector3q::Zero();
                for (auto fv : mesh.face_vertices(f)) {
                    p += mesh.position(fv);
                    q += Q[fv];
                }
                p /= 3;
                q /= 3;
                std::pair<Eigen::Vector3d, Vector3q> pq = {p, q};
                return pq;
            };
            std::vector<Vertex> Mvertices = link(M, m, m2M, x, x2M, [&](Face f) {
                Vertex Nv(-1);
                if (M.is_boundary(f)) {
                    std::vector<Vertex> vertices;
                    for (auto fv : M.face_vertices(f)) {
                        vertices.push_back(M2N[fv]);
                    }
                    Face Nf = N.face_handle(N.find_halfface(vertices));
                    auto [p, q] = centroid(N, Nf);
                    Nv = N.split_face(Nf, p);
                    NQ[Nv] = q;
                }
                Halfface h0 = M.halfface_handle(f, 0);
                bool bdr0 = bdr[h0];
                Vertex v0 = M.halfface_opposite_vertex(h0);
                Halfface h1 = M.halfface_handle(f, 1);
                bool bdr1 = bdr[h1];
                Vertex v1 = M.halfface_opposite_vertex(h1);
                auto [p, q] = centroid(M, f);
                Vertex Mv = M.split_face(f, p);
                MQ[Mv] = q;
                for (auto vh : M.vertex_halffaces(Mv)) {
                    if (M.halfface_opposite_vertex(vh) == v0) {
                        bdr[vh] = bdr0;
                    }
                    if (M.halfface_opposite_vertex(vh) == v1) {
                        bdr[vh] = bdr1;
                    }
                }
                if (N.is_valid(Nv)) {
                    M2N[Mv] = Nv;
                    N2M[Nv] = Mv;
                }
                return Mv;
            });
            std::vector<Vertex> Nvertices = link(N, n, n2N, x, x2N, [&](Face f) {
                Vertex Mv(-1);
                if (N.is_boundary(f)) {
                    std::vector<Vertex> vertices;
                    for (auto fv : N.face_vertices(f)) {
                        vertices.push_back(N2M[fv]);
                    }
                    Face Mf = M.face_handle(M.find_halfface(vertices));
                    Halfface h0 = M.halfface_handle(Mf, 0);
                    bool bdr0 = bdr[h0];
                    Vertex v0 = M.halfface_opposite_vertex(h0);
                    Halfface h1 = M.halfface_handle(Mf, 1);
                    bool bdr1 = bdr[h1];
                    Vertex v1 = M.halfface_opposite_vertex(h1);
                    auto [p, q] = centroid(M, Mf);
                    Mv = M.split_face(Mf, p);
                    MQ[Mv] = q;
                    for (auto vh : M.vertex_halffaces(Mv)) {
                        if (M.halfface_opposite_vertex(vh) == v0) {
                            bdr[vh] = bdr0;
                        }
                        if (M.halfface_opposite_vertex(vh) == v1) {
                            bdr[vh] = bdr1;
                        }
                    }
                }
                auto [p, q] = centroid(N, f);
                Vertex Nv = N.split_face(f, p);
                NQ[Nv] = q;
                if (M.is_valid(Mv)) {
                    N2M[Nv] = Mv;
                    M2N[Mv] = Nv;
                }
                return Nv;
            });
            ASSERT(Mvertices.size() == Nvertices.size());
            for (int k = 0; k < Mvertices.size(); k++) {
                Vertex Mv = Mvertices[k];
                Vertex Nv = Nvertices[k];
                M2N[Mv] = Nv;
                N2M[Nv] = Mv;
            }
        }

        Mesh& S = components[i].first;
        Mesh& S2 = components[i].second;
        auto M2S = M.property<Vertex, Vertex>("", Vertex(-1));
        auto S2M = S.property<Vertex, Vertex>("", Vertex(-1));
        auto SQ = S.property<Vertex, Vector3q>("Q");
        auto S2Q = S2.property<Vertex, Vector3q>("Q");
        std::queue<Cell> queue;
        auto enqueued = M.property<Cell, bool>("", false);
        for (auto h : M.halffaces()) {
            if (bdr[h]) {
                Cell c = M.incident_cell(h);
                queue.push(c);
                enqueued[c] = true;
                break;
            }
        }
        while (!queue.empty()) {
            Cell c = queue.front();
            queue.pop();
            std::vector<Vertex> vertices;
            for (auto cv : M.tet_vertices(c)) {
                if (!S.is_valid(M2S[cv])) {
                    Vertex v = S.add_vertex(M.position(cv));
                    M2S[cv] = v;
                    S2M[v] = cv;
                    SQ[v] = MQ[cv];
                    S2Q[S2.add_vertex(Eigen::Vector3d::Zero())] = Vector3q::Zero();
                }
                vertices.push_back(M2S[cv]);
            }
            S.add_cell(vertices, true);
            S2.add_cell(vertices, true);
            M.remove_cell(c);
            for (auto ch : M.cell_halffaces(c)) {
                Halfface h = M.opposite_halfface_handle(ch);
                if (!bdr[ch] && !M.is_boundary(h)) {
                    Cell cc = M.incident_cell(h);
                    if (!enqueued[cc]) {
                        queue.push(cc);
                        enqueued[cc] = true;
                    }
                }
            }
        }
        for (auto v : S2.vertices()) {
            if (S2.is_boundary(v)) {
                S2.set_position(v, N.position(M2N[S2M[v]]));
                S2Q[v] = NQ[M2N[S2M[v]]];
            }
        }

        std::vector<Cell> cells;
        for (auto c : N.cells()) {
            if (cmp[c] == index) {
                cells.push_back(c);
            }
        }
        N.remove_cells(cells);

        i++;
    }

    std::cout << "Finished" << std::endl;
    return components;
}
