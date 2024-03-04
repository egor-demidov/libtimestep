//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_ROTATIONAL_STEP_HANDLER_H
#define INTEGRATORS_ROTATIONAL_STEP_HANDLER_H

#include <cstddef>

// Method that handles increments to position, velocity, angle, and angular velocity
// This is a generic implementation that merely adds the increments
// but when needed, more sophisticated step handlers can be implemented
// For example, in cases where increments are needed for additional computations
// custom step handlers will be needed
template <typename field_container_t, typename field_value_t>
struct rotational_step_handler {

    // This method increments the specified value in the x buffer
    void increment_x(long n,                                                                              // index of the value to increment
                     field_value_t const & dx,                                                              // value of the position increment
                     typename field_container_t::iterator x_begin_itr,                                      // iterator pointing to the start of the x buffer
                     typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {   // iterator pointing to the start of the alpha buffer

        *(x_begin_itr + n) += dx;
    }

    // This method increments the specified value in the v buffer
    void increment_v(long n,                                                                              // index of the value to increment
                     field_value_t const & dv,                                                              // value of the velocity increment
                     typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                     typename field_container_t::iterator v_begin_itr,                                      // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {   // iterator pointing to the start of the alpha buffer

        *(v_begin_itr + n) += dv;
    }

    // This method increments the specified value in the theta buffer
    void increment_theta(long n,                                                                              // index of the value to increment
                         field_value_t const & dtheta,                                                          // value of the angle increment
                         typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                         typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                         typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                         typename field_container_t::iterator theta_begin_itr,                                  // iterator pointing to the start of the theta buffer
                         typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                         typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {   // iterator pointing to the start of the alpha buffer

        *(theta_begin_itr + n) += dtheta;
    }

    // This method increments the specified value in the omega buffer
    void increment_omega(long n,                                                                              // index of the value to increment
                         field_value_t const & domega,                                                          // value of the angular velocity increment
                         typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                         typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                         typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                         typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                         typename field_container_t::iterator omega_begin_itr,                                  // iterator pointing to the start of the omega buffer
                         typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) const {   // iterator pointing to the start of the alpha buffer

        *(omega_begin_itr + n) += domega;
    }
};

#endif //INTEGRATORS_ROTATIONAL_STEP_HANDLER_H
