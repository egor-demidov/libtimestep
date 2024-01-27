//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_ROTATIONAL_FORWARD_EULER_H
#define INTEGRATORS_ROTATIONAL_FORWARD_EULER_H

// Integrator template that implements the Forward Euler integration scheme (for rotating systems)
// using a two-term Taylor expansion for position and angle to reduce truncation error
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

    // Class constructor
    //
    // Notes:
    // Iterators passed to this constructor must remain valid for the duration of use of this integrator
    // x, v, a, theta, omega, and alpha buffers must be of the same size
    // The acceleration functor and the step handler must exist for the duration of use of this integrator
    rotational_forward_euler(functor_t & acceleration_functor,                                      // reference to a functor that computes translational and angular accelerations
                             real_t t0,                                                             // reference to a functor that computes acceleration
                              typename field_container_t::iterator x_begin,                         // iterator pointing to the start of the x buffer
                              typename field_container_t::iterator x_end,                           // iterator pointing to the end of the x buffer
                              typename field_container_t::iterator v_begin,                         // iterator pointing to the start of the v buffer
                              typename field_container_t::iterator a_begin,                         // iterator pointing to the start of the a buffer
                              typename field_container_t::iterator theta_begin,                     // iterator pointing to the start of the theta buffer
                              typename field_container_t::iterator omega_begin,                     // iterator pointing to the start of the omega buffer
                              typename field_container_t::iterator alpha_begin,                     // iterator pointing to the start of the alpha buffer
                              step_handler_t<field_container_t, field_value_t> & step_handler) :    // reference to an object that handles incrementing positions and velocities

          // Call the superclass constructor
          rotational_integrator<field_container_t, field_value_t, real_t, functor_t, step_handler_t>(acceleration_functor, t0,
               x_begin, x_end, v_begin, a_begin, theta_begin, omega_begin, alpha_begin, step_handler) {}

    // Perform one time step
    void do_step(real_t dt /*time step*/) {
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

            this->step_handler.increment_x(n, v*dt + 0.5*a*dt*dt, this->x_begin_itr, this->v_begin_itr, this->a_begin_itr, this->theta_begin_itr, this->omega_begin_itr, this->alpha_begin_itr);
            this->step_handler.increment_v(n, a*dt, this->x_begin_itr, this->v_begin_itr, this->a_begin_itr, this->theta_begin_itr, this->omega_begin_itr, this->alpha_begin_itr);
            this->step_handler.increment_theta(n, omega*dt + 0.5*alpha*dt*dt, this->x_begin_itr, this->v_begin_itr, this->a_begin_itr, this->theta_begin_itr, this->omega_begin_itr, this->alpha_begin_itr);
            this->step_handler.increment_omega(n, alpha*dt, this->x_begin_itr, this->v_begin_itr, this->a_begin_itr, this->theta_begin_itr, this->omega_begin_itr, this->alpha_begin_itr);
        }
    }
};

#endif //INTEGRATORS_ROTATIONAL_FORWARD_EULER_H
