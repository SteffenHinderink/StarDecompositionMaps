#pragma once

#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>
#include <Eigen/Dense>

Eigen::Vector3d kernel_chebyshev_center(const std::vector<Eigen::Vector3d>& positions, const std::vector<Eigen::Vector3d>& normals);
