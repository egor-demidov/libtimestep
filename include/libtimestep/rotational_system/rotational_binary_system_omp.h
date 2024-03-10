//
// Created by egor on 2/14/24.
//

#ifndef INTEGRATORS_ROTATIONAL_BINARY_SYSTEM_OMP_H
#define INTEGRATORS_ROTATIONAL_BINARY_SYSTEM_OMP_H

#include "rotational_system.h"

#include <omp.h>

// This is a base class for a second order rotational system where accelerations depend on binary
// interactions between fields
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
    typename acceleration_handler_t,
    bool have_unary_force>
class rotational_binary_system_omp : public rotational_generic_system<field_value_t, real_t, integrator_t, step_handler_t,
        rotational_binary_system_omp<field_value_t, real_t, integrator_t, step_handler_t, acceleration_handler_t, have_unary_force>> {
public:
    typedef std::vector<field_value_t> field_container_t;
    typedef std::vector<long> index_container_t;

    // Class constructor
    //
    // Notes:
    // The acceleration handler and the step handler must exist for the duration of use of this second order system
    // This class itself acts as the acceleration functor for the integrator
    rotational_binary_system_omp(field_container_t x0,                                              // container with initial positions
                             field_container_t v0,                                              // container with initial velocities
                             field_container_t theta0,                                          // container with initial angles
                             field_container_t omega0,                                          // container with initial angular velocities
                             real_t t0,                                                         // integration start time
                             field_value_t field_zero,                                          // zero value of the primary field type used
                             real_t real_zero,                                                  // zero value of the real number type used
                             acceleration_handler_t & acceleration_handler,                     // reference to the object that handles calculating accelerations FOR EACH BINARY INTERACTION
                             step_handler_t<field_container_t, field_value_t> & step_handler) : // reference to an object that handles incrementing positions and velocities

         // Call the superclass constructor
         rotational_generic_system<field_value_t, real_t, integrator_t, step_handler_t, rotational_binary_system_omp>(
            std::move(x0), std::move(v0), std::move(theta0), std::move(omega0), t0, field_zero, real_zero, *this, step_handler),
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

        #pragma omp parallel for default(none) shared(t)
        for (long i = 0; i < (long) this->indices.size(); i ++) {
            for (long j = 0; j < (long) this->indices.size(); j ++) {
                if (i == j) [[unlikely]]
                    continue;

                auto [a_i_new, alpha_i_new] = acceleration_handler.compute_accelerations(i, j, this->get_x(), this->get_v(), this->get_theta(), this->get_omega(), t);

                this->a[i] += a_i_new;
                this->alpha[i] += alpha_i_new;
            }

            // This is a compile-time conditional
            if constexpr (have_unary_force) {
                auto [a_i_new, alpha_i_new] = acceleration_handler.compute_accelerations(i, this->get_x(), this->get_v(), this->get_theta(), this->get_omega(), t);

                this->a[i] += a_i_new;
                this->alpha[i] += alpha_i_new;
            }
        }
    }

private:
    acceleration_handler_t & acceleration_handler;
};


#endif //INTEGRATORS_ROTATIONAL_BINARY_SYSTEM_OMP_H
