#include "sdm.h"

std::vector<Vector3q> decompose(Mesh& mesh) {
    auto Q = mesh.property<Vertex, Vector3q>("Q");
    auto normal = mesh.property<Halfface, Vector3q>();
    for (auto f : mesh.faces()) {
        Halfface h = mesh.halfface_handle(f, 0);
        std::vector<Vector3q> triangle;
        for (auto hv : mesh.halfface_vertices(h)) {
            triangle.push_back(Q[hv]);
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        normal[h] = n;
        normal[mesh.opposite_halfface_handle(h)] = -n;
    }
    std::vector<Cell> cells;
    for (auto c : mesh.cells()) {
        cells.push_back(c);
    }
    std::default_random_engine gen(0);
    std::shuffle(cells.begin(), cells.end(), gen);
    std::vector<Vector3q> centers;
    auto cmp = mesh.property<Cell, int>("cmp", -1);
    int i = 0;
    for (auto cell : cells) {
        if (cmp[cell] != -1) {
            continue;
        }
        Vector3q center = Vector3q::Zero();
        for (auto cv : mesh.tet_vertices(cell)) {
            center += Q[cv];
        }
        center /= 4;
        std::set<Halfface> surface;
        std::queue<Cell> queue;
        auto enqueued = mesh.property<Cell, bool>("", false);
        queue.push(cell);
        enqueued[cell] = true;
        while (!queue.empty()) {
            Cell c = queue.front();
            queue.pop();
            enqueued[c] = false;
            if (cmp[c] != -1) {
                continue;
            }
            std::set<Halfface> surface_extension;
            for (auto ch : mesh.cell_halffaces(c)) {
                Halfface h = mesh.opposite_halfface_handle(ch);
                if (mesh.is_boundary(h) || cmp[mesh.incident_cell(h)] != i) {
                    surface_extension.insert(h);
                }
            }
            bool visible = true;
            for (auto h : surface_extension) {
                if ((Q[*(mesh.halfface_vertices(h).first)] - center).dot(normal[h]) <= 0) { // (a - c)^T n
                    visible = false;
                    break;
                }
            }
            if (!visible) {
                continue;
            }
            cmp[c] = i;
            for (auto ch : mesh.cell_halffaces(c)) {
                surface.erase(ch);
            }
            for (auto h : surface_extension) {
                surface.insert(h);
            }
            std::vector<Eigen::Vector3d> positions;
            std::vector<Eigen::Vector3d> normals;
            for (auto h : surface) {
                positions.push_back(mesh.position(*(mesh.halfface_vertices(h).first)));
                normals.push_back(normal[h].unaryExpr([](mpq_class x) { return x.get_d(); }).normalized()); // (0, 0, 0) stays (0, 0, 0)
            }
            Vector3q new_center = kernel_chebyshev_center(positions, normals).cast<mpq_class>();
            bool new_center_valid = true;
            for (auto h : surface) {
                std::vector<Vector3q> triangle;
                for (auto hv : mesh.halfface_vertices(h)) {
                    triangle.push_back(Q[hv]);
                }
                Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
                if ((triangle[0] - new_center).dot(n) <= 0) {
                    new_center_valid = false;
                    break;
                }
            }
            if (new_center_valid) {
                center = new_center;
            }
            for (auto h : surface) {
                if (!mesh.is_boundary(h)) {
                    Cell hc = mesh.incident_cell(h);
                    if (!enqueued[hc]) {
                        queue.push(hc);
                        enqueued[hc] = true;
                    }
                }
            }
        }
        centers.push_back(center);
        i++;
    }
    return centers;
}
