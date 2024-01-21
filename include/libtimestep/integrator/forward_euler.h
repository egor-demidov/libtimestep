//
// Created by egor on 1/18/24.
//

#ifndef INTEGRATORS_FORWARD_EULER_H
#define INTEGRATORS_FORWARD_EULER_H

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
class forward_euler : public integrator<field_container_t, field_value_t, real_t, functor_t, step_handler_t> {
public:
    forward_euler(functor_t & acceleration_functor, real_t t0,
                  typename field_container_t::iterator x_begin,
                  typename field_container_t::iterator x_end,
                  typename field_container_t::iterator v_begin,
                  typename field_container_t::iterator a_begin,
                  step_handler_t<field_container_t, field_value_t> step_handler) :
            integrator<field_container_t, field_value_t, real_t, functor_t, step_handler_t>(acceleration_functor, t0, x_begin, x_end, v_begin, a_begin, std::move(step_handler)) {}

    void do_step(real_t dt) override {
        // Re-compute the acceleration
        this->update_acceleration();

        // Increment time
        this->t += dt;

        // Integrate position and velocity
        for (size_t n = 0; n < this->x_end_itr - this->x_begin_itr; n ++) {
            field_value_t const & v = *(this->v_begin_itr + n);
            field_value_t const & a = *(this->a_begin_itr + n);

            this->step_handler.increment_x(n, v*dt + 0.5*a*dt*dt);
            this->step_handler.increment_v(n, a*dt);
        }
    }
};

#endif //INTEGRATORS_FORWARD_EULER_H
