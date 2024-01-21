//
// Created by egor on 1/18/24.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <execution>
#include <random>
#include <chrono>

#include <Eigen/Eigen>
#include <libtimestep/integrator/integrator.h>
//#include <libtimestep/system/system.h>

#include "write_vtk.h"
#include "compute_energy.h"

class binary_system {
public:
    // k - stiffness constant
    // g - attraction constant
    // r - particle radius
    // m - particle mass
    // gamma_c - visco-elastic damping
    binary_system(double k, double g, double r, double m, double gamma_c) : k(k), g(g), r(r), m(m), gamma_c(gamma_c) {}

    binary_system(binary_system const &) = delete;

    void operator() (std::vector<Eigen::Vector3d>::const_iterator x_begin,
                     std::vector<Eigen::Vector3d>::const_iterator x_end,
                     std::vector<Eigen::Vector3d>::const_iterator v_begin [[maybe_unused]],
                     std::vector<Eigen::Vector3d>::iterator a_begin,
                     double t) const {

        std::vector<size_t> indices(x_end - x_begin);
        std::iota(indices.begin(), indices.end(), 0);

        std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [&indices, &x_begin, &v_begin, &a_begin, this] (size_t i) {
            auto const & x_i = *(x_begin + i);
            auto const & v_i = *(v_begin + i);
            auto & a_i = *(a_begin + i);

            a_i = Eigen::Vector3d::Zero();

            std::for_each(indices.begin(), indices.end(), [i, &x_i, &v_i, &a_i, &x_begin, &v_begin, this] (size_t j) {
                if (i == j)
                    return;

                auto const & x_j = *(x_begin + j);
                auto const & v_j = *(v_begin + j);

                a_i += this->compute_attraction(x_i, x_j) / m;
                a_i += this->compute_elasticity(x_i, x_j, v_i, v_j) / m;
            });
        });
    }

private:
    [[nodiscard]]
    Eigen::Vector3d compute_elasticity(Eigen::Vector3d const & x1, Eigen::Vector3d const & x2,
                                       Eigen::Vector3d const & v1, Eigen::Vector3d const & v2) const {
        Eigen::Vector3d distance = x2 - x1;
        double distance_norm = distance.norm();
        double overlap = distance_norm - 2.0 * this->r;

        if (overlap >= 0.0)
            return Eigen::Vector3d::Zero();

        Eigen::Vector3d n = distance.normalized();
        double normalRelativeVelocity = (v2 - v1).dot(n);

        return (this->k * overlap + this->gamma_c * normalRelativeVelocity) * n;
    }

    [[nodiscard]]
    Eigen::Vector3d compute_attraction(Eigen::Vector3d const & x1, Eigen::Vector3d const & x2) const {
        return this->g * (x2 - x1).normalized();
    }

    const double k, g, r, m, gamma_c;
};

bool do_overlap(std::vector<Eigen::Vector3d> const & particles, Eigen::Vector3d const & x, double r_part) {
    return std::any_of(particles.begin(), particles.end(), [&x, r_part] (auto const & particle) -> bool {
        return (particle - x).norm() <= 2.0 * r_part;
    });
}

int main() {
    const double dt = 0.001;                        // Integration time step
    const double t_tot = 50.0;                      // Duration of the simulation
    const auto n_steps = size_t(t_tot / dt);        // Number of time steps
    const size_t n_dumps = 300;                     // Total number of data dumps
    const size_t dump_period = n_steps / n_dumps;   // Number of time steps between dumps
    const double r_part = 0.1;                      // Radius of a particle
    const double k = 1000.0;                        // Elastic stiffness of aa particle
    const double m = 1.0;                           // Mass of a particle
    const double g = 0.2;                           // Attraction acceleration between particles
    const double gamma_c = 0.05;                    // Elastic (collision) damping coefficient

    std::vector<Eigen::Vector3d> x, v, a;

    // Populate the vectors ...
    // Random particles will be generated in a box of size 1 around the origin
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (size_t i = 0; i < 700; i ++) {
        Eigen::Vector3d x_part;
        do {
            x_part = {dist(mt), dist(mt), dist(mt)};
        } while (do_overlap(x, x_part, r_part));
        x.emplace_back(x_part);
    }

    // Print the number of particles in the system
    std::cout << x.size() << " particles have been created" << std::endl;
    v.resize(x.size(), Eigen::Vector3d::Zero());
    a.resize(x.size(), Eigen::Vector3d::Zero());

    binary_system system(k, g, r_part, m, gamma_c);
    velocity_verlet_half<std::vector<Eigen::Vector3d>, Eigen::Vector3d, double, binary_system>
            integrator(system, 0.0, x.begin(), x.end(), v.begin(), a.begin());

    std::vector<std::vector<Eigen::Vector3d>> solution_buffer;
    std::vector<double> time_buffer;

    // Start the execution timer
    auto start_time = std::chrono::high_resolution_clock::now();

    double t = 0.0;
    for (size_t n = 0; n < n_steps; n ++) {

        if (n % dump_period == 0) {
            write_vtk(x, r_part, n / dump_period, "data");

//            std::cout << "Kinetic energy: " << compute_kinetic_energy(v, m) << std::endl;
            std::cout << "Linear momentum: " << compute_linear_momentum(v, m) << std::endl;
        }

        time_buffer.emplace_back(t);
        solution_buffer.emplace_back(x);

        integrator.do_step(dt);

        t += dt;
    }

    // End the exceution times
    auto end_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    std::cout << elapsed_time << std::endl;

    return 0;
}
