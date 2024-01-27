//
// Created by egor on 1/18/24.
//

#ifndef INTEGRATORS_FORWARD_EULER_H
#define INTEGRATORS_FORWARD_EULER_H

// Integrator template that implements the Forward Euler integration scheme
// For position, uses a two-term Taylor expansion to reduce truncation error
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

    // Class constructor
    //
    // Notes:
    // Iterators passed to this constructor must remain valid for the duration of use of this integrator
    // x, v, and a buffers must be of the same size
    // The acceleration functor and the step handler must exist for the duration of use of this integrator
    forward_euler(functor_t & acceleration_functor,                                     // reference to a functor that computes acceleration
                  real_t t0,                                                            // time at which integration starts
                  typename field_container_t::iterator x_begin,                         // iterator pointing to the start of the x buffer
                  typename field_container_t::iterator x_end,                           // iterator pointing to the end of the x buffer
                  typename field_container_t::iterator v_begin,                         // iterator pointing to the start of the v buffer
                  typename field_container_t::iterator a_begin,                         // iterator pointing to the start of the a buffer
                  step_handler_t<field_container_t, field_value_t> & step_handler) :    // reference to an object that handles incrementing positions and velocities

          // Call the superclass constructor
          integrator<field_container_t, field_value_t, real_t, functor_t, step_handler_t>(
            acceleration_functor, t0, x_begin, x_end, v_begin, a_begin, step_handler) {}

    // Perform one time step
    void do_step(real_t dt /*time step*/) {
        // Re-compute the acceleration
        this->update_acceleration();

        // Increment time
        this->t += dt;

        // Integrate position and velocity
        for (size_t n = 0; n < this->x_end_itr - this->x_begin_itr; n ++) {
            field_value_t const & v = *(this->v_begin_itr + n);
            field_value_t const & a = *(this->a_begin_itr + n);

            this->step_handler.increment_x(n, v*dt + 0.5*a*dt*dt, this->x_begin_itr, this->v_begin_itr, this->a_begin_itr);
            this->step_handler.increment_v(n, a*dt, this->x_begin_itr, this->v_begin_itr, this->a_begin_itr);
        }
    }
};

#endif //INTEGRATORS_FORWARD_EULER_H
