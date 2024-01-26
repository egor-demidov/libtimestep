//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_STEP_HANDLER_H
#define INTEGRATORS_STEP_HANDLER_H

#include <cstddef>

// Method that handles increments to position and velocity
// This is a generic implementation that merely adds the increments
// but when needed, more sophisticated step handlers can be implemented
// For example, in cases where increments are needed for additional computations
// custom step handlers will be needed
template <typename field_container_t, typename field_value_t>
struct step_handler {

    // This method increments the specified value in the x buffer
    void increment_x(size_t n,                                                                          // index of the value to increment
                    field_value_t const & dx,                                                           // value of the position increment
                    typename field_container_t::iterator x_begin_itr,                                   // iterator pointing to the start of the x buffer
                    typename field_container_t::const_iterator v_begin_itr [[maybe_unused]]) const {    // iterator pointing to the start of the v buffer

        *(x_begin_itr + n) += dx;
    }

    // This method increments the specified value in the v buffer
    void increment_v(size_t n,                                                                      // index of the value to increment
                     field_value_t const & dv,                                                      // value of the velocity increment
                     typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],       // iterator pointing to the start of the x buffer
                     typename field_container_t::iterator v_begin_itr) const {                      // iterator pointing to the start of the v buffer

        *(v_begin_itr + n) += dv;
    }
};

#endif //INTEGRATORS_STEP_HANDLER_H
