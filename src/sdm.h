#pragma once

#include <random>

#include <Eigen/Sparse>

#include "assertion.h"
#include "bary.h"
#include "lp.h"
#include "mesh.h"
#include "priority_queue.h"
#include "vectorq.h"

std::vector<std::pair<Mesh, Mesh>> sdm(Mesh& M, Mesh& N);

std::vector<Vector3q> decompose(Mesh& mesh);

bool cut(Mesh& mesh, std::vector<Vertex> loop, TriMesh& result);

void align(TriMesh& m, TriMesh& n, TriMesh& x);
