//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_STEP_HANDLER_H
#define INTEGRATORS_STEP_HANDLER_H

#include <cstddef>

template <typename field_container_t, typename field_value_t>
struct step_handler {
public:
    void increment_x(size_t n, field_value_t const & x, typename field_container_t::iterator x_begin_itr) const {
        *(x_begin_itr + n) += x;
    }

    void increment_v(size_t n, field_value_t const & v, typename field_container_t::iterator v_begin_itr) const {
        *(v_begin_itr + n) += v;
    }
};

#endif //INTEGRATORS_STEP_HANDLER_H
