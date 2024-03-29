//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_ROTATIONAL_UNARY_SYSTEM_H
#define INTEGRATORS_ROTATIONAL_UNARY_SYSTEM_H

#include "rotational_system.h"

// This is a base class for a simple rotational second order system
template <
    typename field_value_t,
    typename real_t,
    template <
        typename _field_container_t,
        typename _field_value_t,
        typename _real_t,
        typename _functor_t,
        template <
            typename __field_container_t,
            typename __field_value_t>
        typename _step_handler_t>
    typename integrator_t,
    template <
        typename _field_container_t,
        typename _field_value_t>
    typename step_handler_t,
    typename acceleration_handler_t>
class rotational_unary_system : public rotational_generic_system<field_value_t, real_t, integrator_t, step_handler_t,
        rotational_unary_system<field_value_t, real_t, integrator_t, step_handler_t, acceleration_handler_t>> {
public:
    typedef std::vector<field_value_t> field_container_t;
    typedef std::vector<long> index_container_t;

    // Class constructor
    //
    // Notes:
    // The acceleration handler and the step handler must exist for the duration of use of this second order system
    // This class itself acts as the acceleration functor for the integrator
    rotational_unary_system(field_container_t x0,                                               // container with initial positions
                            field_container_t v0,                                               // container with initial velocities
                             field_container_t theta0,                                          // container with initial angles
                             field_container_t omega0,                                          // container with initial angular velocities
                             real_t t0,                                                         // integration start time
                             field_value_t field_zero,                                          // zero value of the primary field type used
                             real_t real_zero,                                                  // zero value of the real number type used
                             acceleration_handler_t & acceleration_handler,                     // reference to the object that handles calculating accelerations FOR EACH FIELD
                             step_handler_t<field_container_t, field_value_t> & step_handler) : // reference to an object that handles incrementing positions and velocities

        // Call the superclass constructor
         rotational_generic_system<field_value_t, real_t, integrator_t, step_handler_t, rotational_unary_system>(std::move(x0),
            std::move(v0), std::move(theta0), std::move(omega0), t0, field_zero, real_zero, *this, step_handler),
            acceleration_handler(acceleration_handler) {}


    // This method is called by the integrator to compute accelerations
    void operator() (typename field_container_t::const_iterator x_begin [[maybe_unused]],
                     typename field_container_t::const_iterator x_end [[maybe_unused]],
                     typename field_container_t::const_iterator v_begin [[maybe_unused]],
                     typename field_container_t::iterator a_begin [[maybe_unused]],
                     typename field_container_t::const_iterator theta_begin [[maybe_unused]],
                     typename field_container_t::const_iterator omega_begin [[maybe_unused]],
                     typename field_container_t::iterator alpha_begin [[maybe_unused]],
                     real_t t) {

        this->reset_acceleration_buffers();

        std::for_each(this->indices.begin(), this->indices.end(), [t, this] (size_t i) {
            auto [a_i, alpha_i] = acceleration_handler.compute_acceleration(i, this->get_x(), this->get_v(), this->get_theta(), this->get_omega(), t);
            this->a[i] += a_i;
            this->alpha[i] += alpha_i;
        });
    }

private:
    acceleration_handler_t & acceleration_handler;
};

#endif //INTEGRATORS_ROTATIONAL_UNARY_SYSTEM_H
