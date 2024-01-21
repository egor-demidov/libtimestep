//
// Created by egor on 1/19/24.
//

#include "compute_energy.h"

double compute_kinetic_energy(std::vector<Eigen::Vector3d> const & vs, double m) {
    double val = 0.0;
    for (auto const & v : vs) {
        val += v.dot(v);
    }
    return m * val / 2.0;
}

double compute_linear_momentum(std::vector<Eigen::Vector3d> const & vs, double m) {
    Eigen::Vector3d val = Eigen::Vector3d::Zero();
    for (auto const & v : vs) {
        val += v;
    }
    return m * val.norm();
}
