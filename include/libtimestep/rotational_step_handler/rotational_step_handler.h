//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_ROTATIONAL_STEP_HANDLER_H
#define INTEGRATORS_ROTATIONAL_STEP_HANDLER_H

#include <cstddef>

template <typename field_container_t, typename field_value_t>
struct rotational_step_handler {
    void increment_x(size_t n, field_value_t const & x, typename field_container_t::iterator x_begin_itr) const {
        *(x_begin_itr + n) += x;
    }

    void increment_v(size_t n, field_value_t const & v, typename field_container_t::iterator v_begin_itr) const {
        *(v_begin_itr + n) += v;
    }

    void increment_theta(size_t n, field_value_t const & theta, typename field_container_t::iterator theta_begin_itr) const {
        *(theta_begin_itr + n) += theta;
    }

    void increment_omega(size_t n, field_value_t const & omega, typename field_container_t::iterator omega_begin_itr) const {
        *(omega_begin_itr + n) += omega;
    }
};

#endif //INTEGRATORS_ROTATIONAL_STEP_HANDLER_H
