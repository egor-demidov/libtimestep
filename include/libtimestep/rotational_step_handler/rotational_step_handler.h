//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_ROTATIONAL_STEP_HANDLER_H
#define INTEGRATORS_ROTATIONAL_STEP_HANDLER_H

#include <cstddef>

template <typename field_container_t, typename field_value_t>
struct rotational_step_handler {
    void increment_x(size_t n,
                     field_value_t const & dx,
                     typename field_container_t::iterator x_begin_itr,
                     typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator omega_begin_itr [[maybe_unused]]) const {
        *(x_begin_itr + n) += dx;
    }

    void increment_v(size_t n,
                     field_value_t const & dv,
                     typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator v_begin_itr,
                     typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                     typename field_container_t::iterator omega_begin_itr [[maybe_unused]]) const {
        *(v_begin_itr + n) += dv;
    }

    void increment_theta(size_t n,
                         field_value_t const & dtheta,
                         typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator theta_begin_itr,
                         typename field_container_t::iterator omega_begin_itr [[maybe_unused]]) const {
        *(theta_begin_itr + n) += dtheta;
    }

    void increment_omega(size_t n,
                         field_value_t const & domega,
                         typename field_container_t::iterator x_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator v_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator theta_begin_itr [[maybe_unused]],
                         typename field_container_t::iterator omega_begin_itr) const {
        *(omega_begin_itr + n) += domega;
    }
};

#endif //INTEGRATORS_ROTATIONAL_STEP_HANDLER_H
