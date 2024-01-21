#include <iostream>

#include <Eigen/Eigen>
#include <matplot/matplot.h>

#include <libtimestep/integrator.h>

// This is an example of a functor to be used with an integrator
// It can store the state and history of the simulation
// x and v containers will be populated with position and velocity values
// operator () then needs to compute acceleration for each position-velocity pair and store the computed values at
// the corresponding positions in the container

// This is a damped spring responding to an external force of magnitude g applied at time t=0
struct System {
    // Constructor of the system
    System(double m, double k, double gamma_d, double g) : m(m), k(k), gamma_d(gamma_d), g(g), omega_0(sqrt(k / m)) {}

    // Compute and return the acceleration
    void operator () (std::vector<double>::const_iterator x_begin,
                        std::vector<double>::const_iterator x_end,
                        std::vector<double>::const_iterator v_begin,
                        std::vector<double>::iterator a_begin, double t) const {

        std::transform(x_begin, x_end, v_begin, a_begin, [t, this] (double x, double v) -> double {
            return 1.0 / this->m * (this->g - this->k * x - this->gamma_d * v);
        });
    }

    [[nodiscard]]
    double x_exact(double t) const {
        return g / k * (1.0 - exp(-omega_0 * t) * (omega_0 * t + 1.0));
    }

    double m, k, gamma_d, g, omega_0;
};

double l2_norm_of_error(std::vector<double> const & x_exact, std::vector<double> const & x_numeric) {
    double result = 0.0;
    for (size_t n = 0; n < x_exact.size(); n ++) {
        result += pow(x_exact[n] - x_numeric[n], 2.0);
    }
    return sqrt(result);
}

int main() {
    std::vector<double> x = {0.0};
    std::vector<double> v = {0.0};
    std::vector<double> a = {0.0};
//    std::iterator<std::vector<double>> = vec.begin();

//    static_assert(std::is_same<decltype(*v.begin().base()), double &>::value);

    const double m = 1.0;
    const double k = 10.0;
    const double gamma_d = 2.0 * sqrt(m * k);
    const double g = 1.0;

    const double dt = 0.05;
    const double t_tot = 5.0;
    const auto n_steps = size_t(t_tot / dt);

    std::vector<double> x_euler;
    x_euler.reserve(n_steps);

    System system(m, k, gamma_d, g); // System object must exist for the storage duration of forward_euler
    forward_euler<std::vector<double>, double, double, System> integrator(system, 0.0, x.begin(), x.end(), v.begin(), a.begin());

    for (size_t n = 0; n < n_steps; n ++) {
        x_euler.emplace_back(x[0]);
        integrator.do_step(dt);
    }

    x[0] = 0;
    v[0] = 0;
    a[0] = 0;
    std::vector<double> a_copy = {0.0};
    std::vector<double> x_verlet;
    x_verlet.reserve(n_steps);

    velocity_verlet<std::vector<double>, double, double, System> velocityVerlet(system, 0.0, x.begin(), x.end(), v.begin(), a.begin(), a_copy.begin());

    for (size_t n = 0; n < n_steps; n ++) {
        x_verlet.emplace_back(x[0]);
        velocityVerlet.do_step(dt);
    }

    x[0] = 0;
    v[0] = 0;
    a[0] = 0;
    std::vector<double> x_verlet_half;
    x_verlet_half.reserve(n_steps);

    velocity_verlet_half<std::vector<double>, double, double, System> velocityVerletHalf(system, 0.0, x.begin(), x.end(), v.begin(), a.begin());

    for (size_t n = 0; n < n_steps; n ++) {
        x_verlet_half.emplace_back(x[0]);
        velocityVerletHalf.do_step(dt);
    }

    auto t_span = matplot::linspace(0.0, t_tot, n_steps);
    std::vector<double> x_exact(t_span.size());
    std::transform(t_span.begin(), t_span.end(), x_exact.begin(), [&system] (double t) {
        return system.x_exact(t);
    });

    std::cout << "L^2 norm euler " << l2_norm_of_error(x_exact, x_euler) << std::endl;
    std::cout << "L^2 norm verlet (wikipedia) " << l2_norm_of_error(x_exact, x_verlet) << std::endl;
    std::cout << "L^2 norm verlet (half) " << l2_norm_of_error(x_exact, x_verlet_half) << std::endl;

    auto fig = matplot::figure(false);
    auto ax = fig->current_axes();
    ax->hold(true);
    ax->plot(t_span, x_exact);
    ax->plot(t_span, x_euler);
    ax->plot(t_span, x_verlet);
    ax->legend({"exact", "Euler", "Verlet"});
    fig->show();

    return 0;
}
