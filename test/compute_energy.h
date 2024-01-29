//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_COMPUTE_ENERGY_H
#define INTEGRATORS_COMPUTE_ENERGY_H

#include <vector>

#include <Eigen/Eigen>

double compute_kinetic_energy(std::vector<Eigen::Vector3d> const & vs, double m);

double compute_linear_momentum(std::vector<Eigen::Vector3d> const & vs, double m);

#endif //INTEGRATORS_COMPUTE_ENERGY_H
