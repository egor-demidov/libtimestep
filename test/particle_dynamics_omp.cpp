//
// Created by egor on 1/18/24.
//

#include <vector>
#include <algorithm>
#include <random>
#include <iostream>

// Will be using Eigen for linear algebra
// Any library that overloads arithmetic operators for vectors could work
#include <Eigen/Eigen>

#include <libtimestep/integrator/integrator.h>
#include <libtimestep/step_handler/step_handler.h>
#include <libtimestep/system/binary_system_omp.h>

#include "compute_energy.h"

// Implement a binary second-order granular system
// NOTE: this could have also been implemented in a unary system, but the binary system interface uses
// multithreading and will be faster for larger systems
// Here Velocity Verlet integration scheme because it leads to more stable particulate simulations
class GranularSystem : public binary_system_omp<Eigen::Vector3d, double, velocity_verlet_half, step_handler, GranularSystem, false> {
public:
    GranularSystem(double k, double m, double g, double gamma_c, double r_part,
                   std::vector<Eigen::Vector3d> x0, std::vector<Eigen::Vector3d> v0, double t0) :
            binary_system_omp<Eigen::Vector3d, double, velocity_verlet_half, step_handler, GranularSystem, false>(std::move(x0), std::move(v0),
                                                                                       t0, Eigen::Vector3d::Zero(), 0.0, *this, step_handler_instance),
                    k(k), m(m), g(g), gamma_c(gamma_c), r_part(r_part) {}

    // Compute the acceleration of particle i due to its interaction with particle j
    Eigen::Vector3d compute_acceleration(size_t i, size_t j,
                                         std::vector<Eigen::Vector3d> const & x [[maybe_unused]],
                                         std::vector<Eigen::Vector3d> const & v [[maybe_unused]],
                                         double t [[maybe_unused]]) {
        auto const & x_i = this->x[i];
        auto const & x_j = this->x[j];
        auto const & v_i = this->v[i];
        auto const & v_j = this->v[j];

        return (compute_elasticity(x_i, x_j, v_i, v_j) + compute_attraction(x_i, x_j)) / m;
    }

    Eigen::Vector3d compute_acceleration(long i [[maybe_unused]],
                                         std::vector<Eigen::Vector3d> const & x [[maybe_unused]],
                                         std::vector<Eigen::Vector3d> const & v [[maybe_unused]],
                                         double t [[maybe_unused]]) {

        return Eigen::Vector3d::Zero();
    }

private:
    // Elastic contact force - spring and dashpot
    [[nodiscard]] Eigen::Vector3d compute_elasticity(Eigen::Vector3d const & x1, Eigen::Vector3d const & x2,
                                       Eigen::Vector3d const & v1, Eigen::Vector3d const & v2) const {
        Eigen::Vector3d distance = x2 - x1;
        double distance_norm = distance.norm();
        double overlap = distance_norm - 2.0 * this->r_part;

        if (overlap >= 0.0)
            return Eigen::Vector3d::Zero();

        Eigen::Vector3d n = distance.normalized();
        double normalRelativeVelocity = (v2 - v1).dot(n);

        return (this->k * overlap + this->gamma_c * normalRelativeVelocity) * n;
    }

    // Linear attraction between particles
    [[nodiscard]] Eigen::Vector3d compute_attraction(Eigen::Vector3d const & x1, Eigen::Vector3d const & x2) const {
        return this->g * (x2 - x1).normalized();
    }

    const double k, m, g, gamma_c, r_part;

    step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;
};

bool do_overlap(std::vector<Eigen::Vector3d> const & particles, Eigen::Vector3d const & x, double r_part) {
    return std::any_of(particles.begin(), particles.end(), [&x, r_part] (auto const & particle) -> bool {
        return (particle - x).norm() <= 2.0 * r_part;
    });
}

int main() {
    const double dt = 0.001;                        // Integration time step
    const double t_tot = 50.0;                      // Duration of the simulation
    const auto n_steps = long(t_tot / dt);        // Number of time steps
    const double r_part = 0.1;                      // Radius of a particle
    const double k = 1000.0;                        // Elastic stiffness of aa particle
    const double m = 1.0;                           // Mass of a particle
    const double g = 0.2;                           // Attraction acceleration between particles
    const double gamma_c = 0.2;                    // Elastic (collision) damping coefficient

    const long seed = 0;                          // Deterministic seed for pRNG for reproducibility

    std::vector<Eigen::Vector3d> x0, v0;

    // Populate the vectors ...
    // Random particles will be generated in a box of size 1 around the origin
    std::mt19937_64 mt(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (long i = 0; i < 100; i ++) {
        Eigen::Vector3d x_part;
        do {
            x_part = {dist(mt), dist(mt), dist(mt)};
        } while (do_overlap(x0, x_part, r_part));
        x0.emplace_back(x_part);
    }

    v0.resize(x0.size(), Eigen::Vector3d::Zero());

    GranularSystem system(k, m, g, gamma_c, r_part, x0, v0, 0.0);

    for (long n = 0; n < n_steps; n ++) {
        system.do_step(dt);
    }

    const double target_linear_momentum = 1.39595e-13;
    const double tolerance = 5.0; // Percent
    std::cout << "Actual linear momentum: " << compute_linear_momentum(system.get_v(), m) << std::endl;

    if (compute_linear_momentum(system.get_v(), m) > target_linear_momentum * (1.0 + tolerance / 100.0))
        return EXIT_FAILURE;

    return 0;
}
