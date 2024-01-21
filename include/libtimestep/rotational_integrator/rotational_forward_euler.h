//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_ROTATIONAL_FORWARD_EULER_H
#define INTEGRATORS_ROTATIONAL_FORWARD_EULER_H

// Integrator template that implements the Forward Euler integration scheme
template <
    typename field_container_t,
    typename field_value_t,
    typename real_t,
    typename functor_t,
    template <
        typename _field_container_t,
        typename _field_value_t>
    typename step_handler_t>
class rotational_forward_euler : public rotational_integrator<field_container_t, field_value_t, real_t, functor_t, step_handler_t> {
public:
    rotational_forward_euler(functor_t & acceleration_functor, real_t t0,
                  typename field_container_t::iterator x_begin,
                  typename field_container_t::iterator x_end,
                  typename field_container_t::iterator v_begin,
                  typename field_container_t::iterator a_begin,
                  typename field_container_t::iterator theta_begin,
                  typename field_container_t::iterator omega_begin,
                  typename field_container_t::iterator alpha_begin,
                  step_handler_t<field_container_t, field_value_t> step_handler) :
            rotational_integrator<field_container_t, field_value_t, real_t, functor_t, step_handler_t>(acceleration_functor, t0,
                                   x_begin, x_end, v_begin, a_begin, theta_begin, omega_begin, alpha_begin, std::move(step_handler)) {}

    void do_step(real_t dt) override {
        // Re-compute the acceleration
        this->update_acceleration();

        // Increment time
        this->t += dt;

        // Integrate position and velocity
        for (size_t n = 0; n < this->x_end_itr - this->x_begin_itr; n ++) {
            field_value_t const v = *(this->v_begin_itr + n);
            field_value_t const & a = *(this->a_begin_itr + n);
            field_value_t const & omega = *(this->omega_begin_itr + n);
            field_value_t const & alpha = *(this->alpha_begin_itr + n);

            this->step_handler.increment_x(n, v*dt + 0.5*a*dt*dt);
            this->step_handler.increment_v(n, a*dt);
            this->step_handler.increment_theta(n, omega*dt + 0.5*alpha*dt*dt);
            this->step_handler.increment_omega(n, alpha*dt);
        }
    }
};

#endif //INTEGRATORS_ROTATIONAL_FORWARD_EULER_H
