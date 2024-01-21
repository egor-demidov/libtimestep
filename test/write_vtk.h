//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_WRITE_VTK_H
#define INTEGRATORS_WRITE_VTK_H

#include <vector>
#include <Eigen/Eigen>

void write_vtk(std::vector<Eigen::Vector3d> const & xs, double r_part, size_t suffix, std::string const & path);

#endif //INTEGRATORS_WRITE_VTK_H
