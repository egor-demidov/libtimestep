//
// Created by egor on 1/21/24.
//

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <libtimestep/integrator/integrator.h>
#include <libtimestep/step_handler/step_handler.h>
#include <libtimestep/system/unary_system.h>

// Implement a unary second-order system
class OscillatorSystem : public unary_system<double, double, forward_euler, step_handler, OscillatorSystem> {
public:
    OscillatorSystem(double k, double m, double gamma_d,
                     std::vector<double> x0, std::vector<double> v0, double t0) :
            unary_system<double, double, forward_euler, step_handler, OscillatorSystem>(std::move(x0), std::move(v0), t0, 0.0, 0.0, *this, step_handler_instance),
            k(k), m(m), gamma_d(gamma_d) {}

    double compute_acceleration(size_t i,
                                std::vector<double> const & x,
                                std::vector<double> const & v,
                                double t [[maybe_unused]]) {
        auto const & x_i = this->get_x()[i];
        auto const & v_i = this->get_v()[i];

        return 1.0 / this->m * (1.0 - this->gamma_d * v_i - this->k * x_i);
    }

private:
    step_handler<std::vector<double>, double> step_handler_instance;
    const double k, m, gamma_d;
};

int main() {
    const double dt = 0.1; // Integration time step
    const double t_tot = 5.0; // Integration span
    const double k = 10.0; // Stiffness
    const double m = 1.0; // Mass
    const double gamma_d = 2.0 * sqrt(m * k); // Critically damped system
    const double omega_0 = sqrt(k / m); // Natural angular frequency
    const auto n_steps = size_t(t_tot / dt); // Number of integration time steps

    std::vector<double> x0 = {0.0};
    std::vector<double> v0 = {0.0};

    OscillatorSystem oscillator(k, m, gamma_d, x0, v0, 0.0);

    std::vector<double> t_span(n_steps + 1);
    std::vector<double> x_numerical(n_steps + 1);
    std::vector<double> x_exact(n_steps + 1);
    t_span[0] = 0.0;
    x_numerical[0] = x0[0];

    for (size_t n = 1; n < n_steps+1; n ++) {
        double t = dt * double(n);
        oscillator.do_step(dt);
        t_span[n] = t;
        x_numerical[n] = oscillator.get_x()[0];
    }

    std::transform(t_span.begin(), t_span.end(), x_exact.begin(), [omega_0, k] (double t) {
        return 1.0 / k * (1.0 - exp(-omega_0 * t) * (omega_0 * t + 1.0));
    });

    // Compute the L^2 norm of error
    double norm = 0.0;
    for (size_t i = 0; i < x_exact.size(); i ++) {
        norm += pow(x_exact[i] - x_numerical[i], 2.0);
    }
    norm = sqrt(norm);

    const double target_norm = 0.0152767;
    const double tolerance = 5.0; // Percent

    if (norm > target_norm * (1.0 + tolerance / 100.0))
        return EXIT_FAILURE;

    return 0;
}
