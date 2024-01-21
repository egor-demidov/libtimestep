//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_STEP_HANDLER_H
#define INTEGRATORS_STEP_HANDLER_H

#include <cstddef>

template <typename field_container_t, typename field_value_t>
class step_handler {
public:
    step_handler(typename field_container_t::iterator x_begin,
                 typename field_container_t::iterator v_begin) : x_begin_itr(x_begin), v_begin_itr(v_begin) {}

    void increment_x(size_t n, field_value_t const & val) const {
        *(this->x_begin_itr + n) += val;
    }

    void increment_v(size_t n, field_value_t const & val) const {
        *(this->v_begin_itr + n) += val;
    }

protected:
    typename field_container_t::iterator x_begin_itr, v_begin_itr;
};

#endif //INTEGRATORS_STEP_HANDLER_H
