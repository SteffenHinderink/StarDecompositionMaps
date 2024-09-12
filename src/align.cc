#include "sdm.h"

void align(TriMesh& m, TriMesh& n, TriMesh& x) {
    auto m2M = m.property<SVertex, Vertex>("m2M");
    auto n2N = n.property<SVertex, Vertex>("n2N");
    PolyMesh mp;
    PolyMesh np;
    auto mbary = mp.property<SVertex, BarycentricVector>();
    auto nbary = np.property<SVertex, BarycentricVector>();
    auto to_poly = [](TriMesh& trimesh, SProperty<SVertex, Vertex>& tri2tet,
                      PolyMesh& polymesh, SProperty<SVertex, BarycentricVector>& bary) {
        auto tri2poly = trimesh.property<SVertex, SVertex>();
        for (auto v : trimesh.vertices()) {
            SVertex pv = polymesh.add_vertex(trimesh.position(v));
            bary[pv] = {{tri2tet[v]}, {1}};
            tri2poly[v] = pv;
        }
        for (auto f : trimesh.faces()) {
            std::vector<SVertex> vertices;
            for (auto fv : trimesh.face_vertices(f)) {
                vertices.push_back(tri2poly[fv]);
            }
            polymesh.add_face(vertices);
        }
    };
    to_poly(m, m2M, mp, mbary);
    to_poly(n, n2N, np, nbary);
    std::set<SFace> mpart;
    for (auto f : mp.faces()) {
        mpart.insert(f);
    }
    std::set<SFace> npart;
    for (auto f : np.faces()) {
        npart.insert(f);
    }
    auto m2n = m.property<SVertex, SVertex>("m2n");
    auto mp2np = mp.property<SVertex, SVertex>();
    auto np2mp = np.property<SVertex, SVertex>();
    for (auto mv : mp.vertices()) {
        if (!mp.is_boundary(mv)) {
            continue;
        }
        SVertex nv = m2n[mv];
        ASSERT(np.is_valid_handle(nv));
        mp2np[mv] = nv;
        np2mp[nv] = mv;
    }

    auto linear_dependence = [](SHalfedge h0, SHalfedge h1, PolyMesh& mesh, SProperty<SVertex, BarycentricVector>& bary) {
        BarycentricVector v0 = bary[mesh.to_vertex_handle(h0)] - bary[mesh.from_vertex_handle(h0)];
        BarycentricVector v1 = bary[mesh.to_vertex_handle(h1)] - bary[mesh.from_vertex_handle(h1)];
        mpq_class s = 0;
        for (int i = 0; i < v0.vertices.size(); i++) {
            bool checked = false;
            for (int j = 0; j < v1.vertices.size(); j++) {
                if (v1.vertices[j] == v0.vertices[i]) {
                    ASSERT(v0.coordinates[i] != 0 && v1.coordinates[j] != 0);
                    mpq_class factor = v1.coordinates[j] / v0.coordinates[i];
                    if (s == 0) {
                        s = factor;
                    } else if (factor != s) {
                        return mpq_class(0);
                    }
                    checked = true;
                }
            }
            if (!checked) {
                return mpq_class(0);
            }
        }
        return s;
    };

    auto split_edge = [](SHalfedge h, PolyMesh& mesh, SProperty<SVertex, BarycentricVector>& bary, mpq_class alpha = 0.5) {
        SVertex v0 = mesh.from_vertex_handle(h);
        SVertex v1 = mesh.to_vertex_handle(h);
        Eigen::Vector3d p = alpha.get_d() * mesh.position(v0) + (1 - alpha.get_d()) * mesh.position(v1);
        BarycentricVector b = alpha * bary[v0] + (1 - alpha) * bary[v1];
        ASSERT(b.vertices.size() == b.coordinates.size() && b.vertices.size() <= 3);
        mpq_class sum = 0;
        for (int i = 0; i < b.vertices.size(); i++) {
            sum += b.coordinates[i];
        }
        ASSERT(sum == 1);
        SVertex v = mesh.add_vertex(p);
        bary[v] = b;
        mesh.split(mesh.edge_handle(h), v);
        return v;
    };

    auto split_face = [](SFace f, PolyMesh& mesh, SProperty<SVertex, BarycentricVector>& bary) {
        Eigen::Vector3d p = Eigen::Vector3d::Zero();
        BarycentricVector b;
        int k = 0;
        for (auto fv : mesh.face_vertices(f)) {
            p += mesh.position(fv);
            b += bary[fv];
            k++;
        }
        p /= k;
        b /= k;
        ASSERT(b.vertices.size() == b.coordinates.size() && b.vertices.size() == 3);
        mpq_class sum = 0;
        for (int i = 0; i < b.vertices.size(); i++) {
            sum += b.coordinates[i];
        }
        ASSERT(sum == 1);
        SVertex v = mesh.add_vertex(p);
        bary[v] = b;
        mesh.split(f, v);
        return v;
    };

    std::function<void(std::set<SFace>&, std::set<SFace>&)> align = [&align, &linear_dependence, &split_edge, &split_face, &m = mp, &mbary, &m2n = mp2np, &n = np, &nbary, &n2m = np2mp](std::set<SFace>& mpart, std::set<SFace>& npart) {
        // Initialize
        auto mf = m.property<SFace, bool>("", false);
        auto mh_boundary = m.property<SHalfedge, bool>("", false);
        auto mv_boundary = m.property<SVertex, bool>("", false);
        auto nf = n.property<SFace, bool>("", false);
        auto nh_boundary = n.property<SHalfedge, bool>("", false);
        auto nv_boundary = n.property<SVertex, bool>("", false);
        for (auto f : mpart) {
            mf[f] = true;
        }
        for (auto f : npart) {
            nf[f] = true;
        }
        int k = 0;
        for (auto h : m.halfedges()) {
            if (!m.is_boundary(m.opposite_halfedge_handle(h)) && mf[m.face_handle(m.opposite_halfedge_handle(h))] && (m.is_boundary(h) || !mf[m.face_handle(h)])) {
                mh_boundary[h] = true;
                SVertex v = m.to_vertex_handle(h);
                mv_boundary[v] = true;
                ASSERT(n2m[m2n[v]] == v);
                k++;
            }
        }
        for (auto h : n.halfedges()) {
            if (!n.is_boundary(n.opposite_halfedge_handle(h)) && nf[n.face_handle(n.opposite_halfedge_handle(h))] && (n.is_boundary(h) || !nf[n.face_handle(h)])) {
                nh_boundary[h] = true;
                SVertex v = n.to_vertex_handle(h);
                nv_boundary[v] = true;
                ASSERT(m2n[n2m[v]] == v);
                k--;
            }
        }
        ASSERT(k == 0);

        // Assign roles
        bool am = mpart.size() >= npart.size();
        PolyMesh& a = am ? m : n;
        SProperty<SFace, bool>& af = am ? mf : nf;
        SProperty<SHalfedge, bool>& ah_boundary = am ? mh_boundary : nh_boundary;
        SProperty<SVertex, bool>& av_boundary = am ? mv_boundary : nv_boundary;
        SProperty<SVertex, BarycentricVector>& abary = am ? mbary : nbary;
        SProperty<SVertex, SVertex>& a2b = am ? m2n : n2m;
        PolyMesh& b = am ? n : m;
        SProperty<SFace, bool>& bf = am ? nf : mf;
        SProperty<SHalfedge, bool>& bh_boundary = am ? nh_boundary : mh_boundary;
        SProperty<SVertex, bool>& bv_boundary = am ? nv_boundary : mv_boundary;
        SProperty<SVertex, BarycentricVector>& bbary = am ? nbary : mbary;
        SProperty<SVertex, SVertex>& b2a = am ? n2m : m2n;

        // Fiedler vector
        auto idx = a.property<SFace, int>();
        int dim = 0;
        for (auto f : a.faces()) {
            if (af[f]) {
                idx[f] = dim;
                dim++;
            }
        }
        if (dim == 1) {
            // Recursion anchor
            ASSERT(m.valence(*mpart.begin()) == n.valence(*npart.begin()));
            return;
        }
        std::vector<Eigen::Triplet<double>> triplets;
        for (auto f : a.faces()) {
            if (!af[f]) {
                continue;
            }
            int n_neighbors = 0;
            for (auto ff : a.ff_range(f)) {
                if (!af[ff]) {
                    continue;
                }
                triplets.push_back(Eigen::Triplet<double>(idx[f], idx[ff], -1));
                n_neighbors++;
            }
            triplets.push_back(Eigen::Triplet<double>(idx[f], idx[f], n_neighbors));
        }
        Eigen::SparseMatrix<double> dl(dim, dim);
        dl.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::VectorXd field(dim);
        {
            // Algorithm for Fiedler vector from Wu et al. 2014
            Eigen::VectorXd u(dim);
            u.fill(1);
            u(0) = 1 + std::sqrt(dim);
            double alpha = dim + std::sqrt(dim);
            Eigen::VectorXd h = dl * (u / alpha);
            double gamma = u.transpose() * (h / alpha);
            Eigen::VectorXd v = h - (gamma / 2.0) * u;
            Eigen::VectorXd r = u.bottomRows(dim - 1);
            Eigen::VectorXd s = v.bottomRows(dim - 1);
            Eigen::VectorXd t(dim - 1);
            t.fill(1);
            for (int i = 0; i < 20; i++) { // Inverse power method
                t = t / t.norm();
                Eigen::VectorXd z = t;
                // Conjugate gradient to solve (L' - rsT - srT) * z = t
                Eigen::VectorXd R = t - (dl.bottomRightCorner(dim - 1, dim - 1) * z - double(s.transpose() * z) * r - double(r.transpose() * z) * s); // t - L2 * z (z as initial solution guess)
                Eigen::VectorXd p = R;
                double rsold = R.transpose() * R;
                for (int j = 0; j < 50; j++) {
                    Eigen::VectorXd Ap = dl.bottomRightCorner(dim - 1, dim - 1) * p - double(s.transpose() * p) * r - double(r.transpose() * p) * s; // L2 * p
                    double alpha = rsold / (p.transpose() * Ap);
                    z = z + alpha * p;
                    R = R - alpha * Ap;
                    double rsnew = R.transpose() * R;
                    if (std::sqrt(rsnew) < 1e-8) {
                        break;
                    }
                    p = R + (rsnew / rsold) * p;
                    rsold = rsnew;
                }
                t = z;
            }
            field(0) = 0;
            field.bottomRows(dim - 1) = t;
            double beta = u.transpose() * (field / alpha);
            field -= beta * u;
        }

        // Path along median
        std::vector<double> field_entries;
        for (int i = 0; i < dim; i++) {
            field_entries.push_back(field(i));
        }
        std::sort(field_entries.begin(), field_entries.end());
        double epsilon = 0.001;
        double threshold = field_entries[0];
        for (int i = dim / 2; i < dim && threshold <= field_entries[0] + epsilon; i++) {
            threshold = field_entries[i];
        }
        std::vector<SVertex> apath;
        auto in_path = a.property<SVertex, bool>("", false);
        SHalfedge path_h;
        for (auto h : a.halfedges()) {
            if (ah_boundary[h]) {
                SHalfedge next;
                for (auto vh : a.voh_range(a.to_vertex_handle(h))) {
                    if (ah_boundary[vh]) {
                        next = vh;
                        break;
                    }
                }
                if (field(idx[a.face_handle(a.opposite_halfedge_handle(next))]) < threshold &&
                    field(idx[a.face_handle(a.opposite_halfedge_handle(h))]) >= threshold) {
                    apath.push_back(a.to_vertex_handle(h));
                    in_path[a.to_vertex_handle(h)] = true;
                    path_h = a.opposite_halfedge_handle(next);
                    break;
                }
            }
        }
        if (apath.size() == 1) {
            do {
                path_h = a.next_halfedge_handle(path_h);
                while (field(idx[a.face_handle(a.opposite_halfedge_handle(path_h))]) < threshold) {
                    path_h = a.next_halfedge_handle(a.opposite_halfedge_handle(path_h));
                }
                ASSERT(field(idx[a.face_handle(path_h)]) < threshold &&
                       field(idx[a.face_handle(a.opposite_halfedge_handle(path_h))]) >= threshold);
                if (in_path[a.to_vertex_handle(path_h)]) {
                    apath.clear();
                    break;
                }
                apath.push_back(a.to_vertex_handle(path_h));
                in_path[a.to_vertex_handle(path_h)] = true;
            } while (!av_boundary[apath.back()]);
        }
        if (apath.size() == 0) {
            // Fallback solution: Boundary-boundary Laplace
            std::vector<SFace> boundary_faces;
            for (auto f : a.faces()) {
                if (!af[f]) {
                    continue;
                }
                for (auto fh : a.fh_range(f)) {
                    if (ah_boundary[a.opposite_halfedge_handle(fh)]) {
                        for (auto bf : boundary_faces) {
                            ASSERT(bf != f);
                        }
                        boundary_faces.push_back(f);
                        break;
                    }
                }
            }
            ASSERT(boundary_faces.size() > 1);
            double min = std::numeric_limits<double>::max();
            SFace fmin;
            double max = std::numeric_limits<double>::lowest();
            SFace fmax;
            for (auto f : boundary_faces) {
                if (field(idx[f]) < min) {
                    min = field(idx[f]);
                    fmin = f;
                }
                if (field(idx[f]) > max) {
                    max = field(idx[f]);
                    fmax = f;
                }
            }
            if (fmin == fmax) {
                fmin = boundary_faces[0];
                fmax = boundary_faces[1];
            }
            triplets.clear();
            for (auto f : a.faces()) {
                if (!af[f]) {
                    continue;
                }
                if (f == fmin || f == fmax) {
                    triplets.push_back(Eigen::Triplet<double>(idx[f], idx[f], 1));
                    continue;
                }
                int n_neighbors = 0;
                for (auto ff : a.ff_range(f)) {
                    if (!af[ff]) {
                        continue;
                    }
                    triplets.push_back(Eigen::Triplet<double>(idx[f], idx[ff], -1));
                    n_neighbors++;
                }
                triplets.push_back(Eigen::Triplet<double>(idx[f], idx[f], n_neighbors));
            }
            dl.setFromTriplets(triplets.begin(), triplets.end());
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dim);
            rhs(idx[fmin]) = -1;
            rhs(idx[fmax]) = 1;
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(dl);
            ASSERT(solver.info() == Eigen::Success);
            field = solver.solve(rhs);
            field_entries.clear();
            for (int i = 0; i < dim; i++) {
                field_entries.push_back(field(i));
            }
            std::sort(field_entries.begin(), field_entries.end());
            threshold = field_entries[0];
            for (int i = dim / 2; i < dim && threshold <= field_entries[0] + epsilon; i++) {
                threshold = field_entries[i];
            }
            apath.clear();
            in_path = a.property<SVertex, bool>("", false);
            path_h = SHalfedge(-1);
            for (auto h : a.halfedges()) {
                if (ah_boundary[h]) {
                    SHalfedge next;
                    for (auto vh : a.voh_range(a.to_vertex_handle(h))) {
                        if (ah_boundary[vh]) {
                            next = vh;
                            break;
                        }
                    }
                    if (field(idx[a.face_handle(a.opposite_halfedge_handle(next))]) < threshold &&
                        field(idx[a.face_handle(a.opposite_halfedge_handle(h))]) >= threshold) {
                        apath.push_back(a.to_vertex_handle(h));
                        in_path[a.to_vertex_handle(h)] = true;
                        path_h = a.opposite_halfedge_handle(next);
                        break;
                    }
                }
            }
            ASSERT(apath.size() == 1);
            bool island = false;
            do {
                path_h = a.next_halfedge_handle(path_h);
                while (field(idx[a.face_handle(a.opposite_halfedge_handle(path_h))]) < threshold) {
                    path_h = a.next_halfedge_handle(a.opposite_halfedge_handle(path_h));
                }
                ASSERT(field(idx[a.face_handle(path_h)]) < threshold &&
                       field(idx[a.face_handle(a.opposite_halfedge_handle(path_h))]) >= threshold);
                if (in_path[a.to_vertex_handle(path_h)]) {
                    island = true;
                    break;
                }
                apath.push_back(a.to_vertex_handle(path_h));
                in_path[a.to_vertex_handle(path_h)] = true;
            } while (!av_boundary[apath.back()]);
            if (island) {
                // Fallback solution 2: Any path
                SVertex s = apath[0];
                auto prev = a.property<SVertex, SVertex>("", SVertex(-1));
                std::queue<SVertex> queue;
                prev[s] = s;
                queue.push(s);
                while (!queue.empty()) {
                    SVertex v = queue.front();
                    queue.pop();
                    for (auto vh : a.voh_range(v)) {
                        if (a.is_boundary(a.edge_handle(vh)) || !af[a.face_handle(vh)] || !af[a.face_handle(a.opposite_halfedge_handle(vh))]) {
                            continue;
                        }
                        SVertex vv = a.to_vertex_handle(vh);
                        if (!a.is_valid_handle(prev[vv])) {
                            prev[vv] = v;
                            queue.push(vv);
                        }
                    }
                }
                SVertex t(-1);
                SVertex vcw = s;
                SVertex vccw = s;
                while (true) {
                    for (auto vh : a.voh_range(vcw)) {
                        if (ah_boundary[vh]) {
                            vcw = a.to_vertex_handle(vh);
                            break;
                        }
                    }
                    if (vcw == vccw) {
                        break;
                    }
                    if (a.is_valid_handle(prev[vcw])) {
                        t = vcw;
                    }
                    for (auto vh : a.voh_range(vccw)) {
                        if (ah_boundary[a.opposite_halfedge_handle(vh)]) {
                            vccw = a.to_vertex_handle(vh);
                            break;
                        }
                    }
                    if (vccw == vcw) {
                        break;
                    }
                    if (a.is_valid_handle(prev[vccw])) {
                        t = vcw;
                    }
                }
                apath = {t};
                do {
                    ASSERT(a.is_valid_handle(apath.back()));
                    apath.push_back(prev[apath.back()]);
                } while (apath.back() != s);
                std::reverse(apath.begin(), apath.end());
            }
        }
        ASSERT(apath.size() > 1 && apath.front() != apath.back());

        SVertex s = a2b[apath.front()];
        SVertex t = a2b[apath.back()];

        // Pre split
        // Mark boundary
        //       t
        //   -1     1
        // ...       ...
        //   -1     1
        //       s
        auto side = b.property<SVertex, int>("", 0);
        SHalfedge h;
        for (auto sh : b.voh_range(s)) {
            if (bh_boundary[sh]) {
                h = sh;
                break;
            }
        }
        while (b.to_vertex_handle(h) != t) {
            side[b.to_vertex_handle(h)] = -1;
            for (auto vh : b.voh_range(b.to_vertex_handle(h))) {
                if (bh_boundary[vh]) {
                    h = vh;
                    break;
                }
            }
        }
        for (auto th : b.voh_range(t)) {
            if (bh_boundary[th]) {
                h = th;
                break;
            }
        }
        while (b.to_vertex_handle(h) != s) {
            side[b.to_vertex_handle(h)] = 1;
            for (auto vh : b.voh_range(b.to_vertex_handle(h))) {
                if (bh_boundary[vh]) {
                    h = vh;
                    break;
                }
            }
        }
        // Edge splits
        std::vector<SEdge> edges;
        for (auto e : b.edges()) {
            if (b.is_boundary(e) || !bf[b.face_handle(b.halfedge_handle(e, 0))] || !bf[b.face_handle(b.halfedge_handle(e, 1))]) {
                continue;
            }
            auto [v0, v1] = b.edge_vertices(e);
            if (side[v0] * side[v1] == -1) {
                edges.push_back(e);
            }
        }
        for (auto e : edges) {
            SVertex v = split_edge(b.halfedge_handle(e, 0), b, bbary);
            bv_boundary[v] = false;
            side[v] = 0;
        }
        // Face splits
        std::vector<std::pair<SHalfedge, SHalfedge>> halfedge_pairs;
        std::vector<SFace> faces;
        for (auto f : b.faces()) {
            if (!bf[f]) {
                continue;
            }
            SHalfedge hl1(-1);
            SHalfedge hr1(-1);
            for (auto fh : b.fh_range(f)) {
                if (side[b.to_vertex_handle(fh)] == -1) {
                    hl1 = fh;
                }
                if (side[b.to_vertex_handle(fh)] == 1) {
                    hr1 = fh;
                }
            }
            if (!b.is_valid_handle(hl1) || !b.is_valid_handle(hr1)) {
                continue;
            }
            while (side[b.to_vertex_handle(hl1)] == -1) {
                hl1 = b.next_halfedge_handle(hl1);
            }
            while (side[b.to_vertex_handle(hr1)] == 1) {
                hr1 = b.next_halfedge_handle(hr1);
            }
            while (side[b.to_vertex_handle(hl1)] != 1) {
                hl1 = b.next_halfedge_handle(hl1);
            }
            while (side[b.to_vertex_handle(hr1)] != -1) {
                hr1 = b.next_halfedge_handle(hr1);
            }
            SHalfedge hl0 = b.prev_halfedge_handle(hl1);
            SHalfedge hr0 = b.prev_halfedge_handle(hr1);
            if (linear_dependence(hl0, hr1, b, bbary) > 0 || linear_dependence(hr0, hl1, b, bbary) > 0) {
                faces.push_back(f);
            } else {
                ASSERT(side[b.to_vertex_handle(hl0)] == 0 && side[b.from_vertex_handle(hr1)] == 0);
                halfedge_pairs.push_back({hl0, hr1});
            }
        }
        SHalfedge st = b.find_halfedge(s, t);
        if (b.is_valid_handle(st)) {
            for (int i = 0; i < 2; i++) {
                if (bh_boundary[b.opposite_halfedge_handle(st)] && bf[b.face_handle(st)]) {
                    for (auto v : b.face_vertices(b.face_handle(st))) {
                        if (v != s && v != t && bv_boundary[v]) {
                            faces.push_back(b.face_handle(st));
                            break;
                        }
                    }
                }
                st = b.opposite_halfedge_handle(st);
            }
        }
        // Face splits with existing vertices
        for (auto hh : halfedge_pairs) {
            SHalfedge h = b.insert_edge(hh.first, hh.second);
            bh_boundary[h] = false;
            bh_boundary[b.opposite_halfedge_handle(h)] = false;
            bf[b.face_handle(h)] = true;
            bf[b.face_handle(b.opposite_halfedge_handle(h))] = true;
        }
        // Face splits inserting new vertex
        for (auto f : faces) {
            SVertex v = split_face(f, b, bbary);
            bv_boundary[v] = false;
            side[v] = 0;
            for (auto vh : b.voh_range(v)) {
                bh_boundary[vh] = false;
                bh_boundary[b.opposite_halfedge_handle(vh)] = false;
            }
            for (auto vf : b.vf_range(v)) {
                bf[vf] = true;
            }
        }

        auto prev = b.property<SVertex, SVertex>("", SVertex(-1));
        std::queue<SVertex> queue;
        prev[s] = s;
        queue.push(s);
        while (!queue.empty()) {
            SVertex v = queue.front();
            queue.pop();
            if (v == t) {
                break;
            }
            for (auto vh : b.voh_range(v)) {
                if (b.is_boundary(b.edge_handle(vh)) || !bf[b.face_handle(vh)] || !bf[b.face_handle(b.opposite_halfedge_handle(vh))]) {
                    continue;
                }
                SVertex vv = b.to_vertex_handle(vh);
                if (bv_boundary[vv] && vv != s && vv != t) {
                    continue;
                }
                if (!b.is_valid_handle(prev[vv])) {
                    prev[vv] = v;
                    queue.push(vv);
                }
            }
        }
        std::vector<SVertex> bpath = {t};
        do {
            ASSERT(b.is_valid_handle(bpath.back()));
            bpath.push_back(prev[bpath.back()]);
        } while (bpath.back() != s);
        std::reverse(bpath.begin(), bpath.end());
        ASSERT(apath.front() == b2a[bpath.front()] && apath.back() == b2a[bpath.back()] &&
               a2b[apath.front()] == bpath.front() && a2b[apath.back()] == bpath.back());

        // Subdivide shorter path
        auto subdivide = [&split_edge](std::vector<SVertex> path, int length, PolyMesh& mesh, SProperty<SVertex, BarycentricVector>& bary) {
            ASSERT(path.size() <= length);
            std::vector<int> n_insert(path.size() - 1, 0);
            for (int i = 0; i < length - path.size(); i++) {
                n_insert[i % n_insert.size()]++;
            }
            std::vector<SVertex> subdivided_path;
            for (int i = 0; i < n_insert.size(); i++) {
                subdivided_path.push_back(path[i]);
                for (int j = 0; j < n_insert[i]; j++) {
                    subdivided_path.push_back(split_edge(mesh.find_halfedge(path[i + 1], subdivided_path.back()), mesh, bary, mpq_class(1, n_insert[i] - j + 1)));
                }
            }
            subdivided_path.push_back(path.back());
            return subdivided_path;
        };
        if (apath.size() < bpath.size()) {
            apath = subdivide(apath, bpath.size(), a, abary);
        } else if (bpath.size() < apath.size()) {
            bpath = subdivide(bpath, apath.size(), b, bbary);
        }
        ASSERT(apath.size() == bpath.size());
        for (int i = 0; i < apath.size(); i++) {
            a2b[apath[i]] = bpath[i];
            b2a[bpath[i]] = apath[i];
        }

        // Divide into 2 halves
        auto fill = [](std::set<SFace>& part, SProperty<SEdge, bool>& border, PolyMesh& mesh) {
            std::queue<SFace> queue;
            auto enqueued = mesh.property<SFace, bool>("", false);
            for (auto f : part) {
                if (!enqueued[f]) {
                    queue.push(f);
                    enqueued[f] = true;
                }
            }
            while (!queue.empty()) {
                SFace f = queue.front();
                queue.pop();
                part.insert(f);
                for (auto fh : mesh.fh_range(f)) {
                    if (border[mesh.edge_handle(fh)]) {
                        continue;
                    }
                    SFace ff = mesh.face_handle(mesh.opposite_halfedge_handle(fh));
                    if (!enqueued[ff]) {
                        queue.push(ff);
                        enqueued[ff] = true;
                    }
                }
            }
        };
        auto mborder = m.property<SEdge, bool>("", false);
        for (auto h : m.halfedges()) {
            if (mh_boundary[h]) {
                mborder[m.edge_handle(h)] = true;
            }
        }
        auto nborder = n.property<SEdge, bool>("", false);
        for (auto h : n.halfedges()) {
            if (nh_boundary[h]) {
                nborder[n.edge_handle(h)] = true;
            }
        }
        std::set<SFace> mxpart;
        std::set<SFace> mypart;
        std::set<SFace> nxpart;
        std::set<SFace> nypart;
        for (int i = 0; i < apath.size() - 1; i++) {
            SHalfedge mh = m.find_halfedge(am ? apath[i] : bpath[i], am ? apath[i + 1] : bpath[i + 1]);
            mborder[m.edge_handle(mh)] = true;
            mxpart.insert(m.face_handle(mh));
            mypart.insert(m.face_handle(m.opposite_halfedge_handle(mh)));
            SHalfedge nh = n.find_halfedge(am ? bpath[i] : apath[i], am ? bpath[i + 1] : apath[i + 1]);
            nborder[n.edge_handle(nh)] = true;
            nxpart.insert(n.face_handle(nh));
            nypart.insert(n.face_handle(n.opposite_halfedge_handle(nh)));
        }
        fill(mxpart, mborder, m);
        fill(mypart, mborder, m);
        fill(nxpart, nborder, n);
        fill(nypart, nborder, n);

        align(mxpart, nxpart);
        align(mypart, nypart);
    };
    align(mpart, npart);

    // Triangulate
    std::queue<std::pair<SFace, SFace>> queue;
    for (auto mf : mp.faces()) {
        std::vector<SVertex> vertices;
        for (auto fv : mp.face_vertices(mf)) {
            vertices.push_back(mp2np[fv]);
        }
        if (vertices.size() > 3) {
            SFace nf = np.face_handle(np.find_halfedge(vertices[0], vertices[1]));
            ASSERT(np.is_valid_handle(nf) && np.valence(nf) == vertices.size());
            queue.push({mf, nf});
        }
    }
    while (!queue.empty()) {
        auto [mf, nf] = queue.front();
        queue.pop();
        int k = mp.valence(mf);
        ASSERT(k == np.valence(nf) && k > 3);
        std::vector<SHalfedge> mhalfedges;
        std::vector<SHalfedge> nhalfedges;
        for (auto h : mp.fh_ccw_range(mf)) {
            mhalfedges.push_back(h);
            nhalfedges.push_back(np.find_halfedge(mp2np[mp.from_vertex_handle(h)], mp2np[mp.to_vertex_handle(h)]));
        }
        ASSERT(mhalfedges.size() == k && nhalfedges.size() == k);
        bool split = false;
        for (int i = 0; i < k && !split; i++) {
            for (int j = i + 1; j < k && !split; j++) {
                if (i == (j + 1) % k || j == (i + 1) % k) {
                    continue;
                }
                if (linear_dependence(mhalfedges[i], mhalfedges[(j + 1) % k], mp, mbary) <= 0 &&
                    linear_dependence(mhalfedges[j], mhalfedges[(i + 1) % k], mp, mbary) <= 0 &&
                    linear_dependence(nhalfedges[i], nhalfedges[(j + 1) % k], np, nbary) <= 0 &&
                    linear_dependence(nhalfedges[j], nhalfedges[(i + 1) % k], np, nbary) <= 0) {
                    SHalfedge mh = mp.insert_edge(mhalfedges[i], mhalfedges[(j + 1) % k]);
                    SHalfedge nh = np.insert_edge(nhalfedges[i], nhalfedges[(j + 1) % k]);
                    split = true;
                    ASSERT(mp.valence(mp.face_handle(mh)) == np.valence(np.face_handle(nh)));
                    if (mp.valence(mp.face_handle(mh)) > 3) {
                        queue.push({mp.face_handle(mh), np.face_handle(nh)});
                    }
                    ASSERT(mp.valence(mp.face_handle(mp.opposite_halfedge_handle(mh))) == np.valence(np.face_handle(np.opposite_halfedge_handle(nh))));
                    if (mp.valence(mp.face_handle(mp.opposite_halfedge_handle(mh))) > 3) {
                        queue.push({mp.face_handle(mp.opposite_halfedge_handle(mh)), np.face_handle(np.opposite_halfedge_handle(nh))});
                    }
                }
            }
        }
        if (!split) {
            SVertex mv = split_face(mf, mp, mbary);
            SVertex nv = split_face(nf, np, nbary);
            mp2np[mv] = nv;
            np2mp[nv] = mv;
        }
    }
    for (auto f : mp.faces()) {
        ASSERT(mp.valence(f) == 3);
    }
    for (auto f : np.faces()) {
        ASSERT(np.valence(f) == 3);
    }

    ASSERT(mp.n_vertices() == np.n_vertices() && mp.n_faces() == np.n_faces());
    x = TriMesh();
    auto x2M = x.property<SVertex, BarycentricVector>("x2M");
    auto x2N = x.property<SVertex, BarycentricVector>("x2N");
    auto mp2x = mp.property<SVertex, SVertex>("", SVertex(-1));
    for (auto f : mp.faces()) {
        std::vector<SVertex> vertices;
        for (auto fv : mp.face_vertices(f)) {
            if (!x.is_valid_handle(mp2x[fv])) {
                ASSERT(np2mp[mp2np[fv]] == fv);
                SVertex v = x.add_vertex(mp.position(fv));
                mp2x[fv] = v;
                x2M[v] = mbary[fv];
                x2N[v] = nbary[mp2np[fv]];
            }
            ASSERT(x.is_valid_handle(mp2x[fv]));
            vertices.push_back(mp2x[fv]);
        }
        ASSERT(x.is_valid_handle(x.add_face(vertices)));
    }
}
