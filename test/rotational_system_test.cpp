//
// Created by egor on 1/19/24.
//

#include <vector>
#include <iostream>
#include <random>
#include <chrono>

#include <Eigen/Eigen>

#include <libtimestep/rotational_integrator/rotational_integrator.h>
#include <libtimestep/rotational_system/rotational_system.h>
#include <libtimestep/rotational_step_handler/rotational_step_handler.h>

#include "../../soot-indentation/energy/compute_energy.h"
#include "write_vtk.h"

class RotationalGranularSystem : public rotational_binary_system<Eigen::Vector3d, double, rotational_velocity_verlet_half, rotational_step_handler, RotationalGranularSystem> {
public:
    RotationalGranularSystem(double k, double g, double r, double m, double gamma_c,
                             std::vector<Eigen::Vector3d> x0, std::vector<Eigen::Vector3d> v0,
                             std::vector<Eigen::Vector3d> theta0, std::vector<Eigen::Vector3d> omega0, double t0) :
            rotational_binary_system<Eigen::Vector3d, double, rotational_velocity_verlet_half, rotational_step_handler, RotationalGranularSystem>(std::move(x0), std::move(v0),
                                                                            std::move(theta0), std::move(omega0),
                                                                            t0, Eigen::Vector3d::Zero(), 0.0, *this, step_handler_instance),
                                                                            k(k), g(g), r(r), m(m), gamma_c(gamma_c) {}

    std::pair<Eigen::Vector3d, Eigen::Vector3d> compute_accelerations(size_t i, size_t j,
                                                                      std::vector<Eigen::Vector3d> const & x,
                                                                      std::vector<Eigen::Vector3d> const & v,
                                                                      std::vector<Eigen::Vector3d> const & theta [[maybe_unused]],
                                                                      std::vector<Eigen::Vector3d> const & omega [[maybe_unused]],
                                                                      double t [[maybe_unused]]) {
        auto const & xi = x[i];
        auto const & xj = x[j];
        auto const & vi = v[i];
        auto const & vj = v[j];

        auto [ai1, alphai1] = compute_elasticity(xi, xj, vi, vj);
        auto [ai2, alphai2] = compute_attraction(xi, xj);

        return std::make_pair(ai1+ai2, alphai1+alphai2);
    }

    std::pair<Eigen::Vector3d, Eigen::Vector3d> compute_accelerations(size_t i [[maybe_unused]],
                                         std::vector<Eigen::Vector3d> const & x [[maybe_unused]],
                                         std::vector<Eigen::Vector3d> const & v [[maybe_unused]],
                                         std::vector<Eigen::Vector3d> const & theta [[maybe_unused]],
                                         std::vector<Eigen::Vector3d> const & omega [[maybe_unused]],
                                         double t [[maybe_unused]]) {

        return std::make_pair(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
    }

private:
    [[nodiscard]]
    std::pair<Eigen::Vector3d, Eigen::Vector3d> compute_elasticity(Eigen::Vector3d const & x1, Eigen::Vector3d const & x2,
                                       Eigen::Vector3d const & v1, Eigen::Vector3d const & v2) const {
        Eigen::Vector3d distance = x2 - x1;
        double distance_norm = distance.norm();
        double overlap = distance_norm - 2.0 * this->r;

        if (overlap >= 0.0)
            return std::make_pair(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());

        Eigen::Vector3d n = distance.normalized();
        double normalRelativeVelocity = (v2 - v1).dot(n);

        return std::make_pair((this->k * overlap + this->gamma_c * normalRelativeVelocity) * n, Eigen::Vector3d::Zero());
    }

    [[nodiscard]]
    std::pair<Eigen::Vector3d, Eigen::Vector3d>  compute_attraction(Eigen::Vector3d const & x1, Eigen::Vector3d const & x2) const {
        return std::make_pair(this->g * (x2 - x1).normalized(), Eigen::Vector3d::Zero());
    }

    const double k, g, r, m, gamma_c;

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;
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

    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    // Populate the vectors ...
    // Random particles will be generated in a box of size 1 around the origin
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (size_t i = 0; i < 700; i ++) {
        Eigen::Vector3d x_part;
        do {
            x_part = {dist(mt), dist(mt), dist(mt)};
        } while (do_overlap(x0, x_part, r_part));
        x0.emplace_back(x_part);
    }

    // Print the number of particles in the system
    std::cout << x0.size() << " particles have been created" << std::endl;
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    RotationalGranularSystem system(k, g, r_part, m, gamma_c, x0, v0, theta0, omega0, 0.0);

    std::vector<std::vector<Eigen::Vector3d>> solution_buffer;
    std::vector<double> time_buffer;

    // Start the execution timer
    auto start_time = std::chrono::high_resolution_clock::now();

    double t = 0.0;
    for (size_t n = 0; n < n_steps; n ++) {

        if (n % dump_period == 0) {
            write_vtk(system.get_x(), r_part, n / dump_period, "data");

//            std::cout << "Kinetic energy: " << compute_kinetic_energy(system.get_v(), m) << std::endl;
            std::cout << "Linear momentum: " << compute_linear_momentum(system.get_v(), m) << std::endl;
        }

        time_buffer.emplace_back(t);
        solution_buffer.emplace_back(system.get_x());

        system.do_step(dt);

        t += dt;
    }

    // End the exceution times
    auto end_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    std::cout << elapsed_time << std::endl;

    return 0;
}

