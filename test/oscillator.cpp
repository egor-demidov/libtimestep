//
// Created by egor on 1/21/24.
//

#include <vector>
#include <iostream>
#include <cmath>

#include <matplot/matplot.h>

#include <libtimestep/integrator/integrator.h>
#include <libtimestep/step_handler/step_handler.h>
#include <libtimestep/system/system.h>

class OscillatorSystem : public unary_system<double, double, forward_euler, step_handler> {
public:
    OscillatorSystem(double k, double m, double gamma_d,
                     std::vector<double> x0, std::vector<double> v0, double t0) :
            unary_system<double, double, forward_euler, step_handler>(std::move(x0), std::move(v0), t0, 0.0, 0.0),
            k(k), m(m), gamma_d(gamma_d) {}

    double compute_acceleration(size_t i, double t [[maybe_unused]]) override {
        auto const & x_i = get_x()[i];
        auto const & v_i = get_v()[i];

        return 1.0 / this->m * (1.0 - this->gamma_d * v_i - this->k * x_i);
    }

private:
    const double k, m, gamma_d;
};

int main() {
    const double dt = 0.01; // Integration time step
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
    t_span[0] = 0.0;
    x_numerical[0] = x0[0];

    for (size_t n = 0; n < n_steps; n ++) {
        double t = dt * double(n);
        oscillator.do_step(dt);
        t_span.emplace_back(t);
        x_numerical.emplace_back(oscillator.get_x()[0]);
    }

    auto fig = matplot::figure();
    auto ax = fig->current_axes();
    ax->plot(t_span, x_numerical);
    fig->show();

    return 0;
}
