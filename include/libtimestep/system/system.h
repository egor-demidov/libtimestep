//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_SYSTEM_H
#define INTEGRATORS_SYSTEM_H

#include <vector>
#include <execution>

#include "../exception/exception.h"

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
        typename step_handler_t>
class generic_system {
public:
    typedef std::vector<field_value_t> field_container_t;
    typedef std::vector<size_t> index_container_t;

    // The containers cannot be resized after the integrator has been instantiated because that would
    // invalidate the iterators
    generic_system(field_container_t x0, field_container_t v0, real_t t0, field_value_t field_zero, real_t real_zero) :
            field_zero(std::move(field_zero)), real_zero(real_zero), x(std::move(x0)), v(std::move(v0)), a(x.size()),
            integrator(*this, t0, std::begin(this->x), std::end(this->x),
                       std::begin(this->v), std::begin(this->a), step_handler_t<field_container_t, field_value_t>(std::begin(this->x), std::begin(this->v))) {

        // Assert that x and v buffers are of the same size
        if (x.size() != v.size())
            throw SizeMismatchException("binary_system(field_container_t, field_container_t)");

        // Initialize index buffer
        indices.resize(x.size());

        // Generate a vector with particle indices
        std::iota(indices.begin(), indices.end(), 0);
    }

    // This method should set all entries in the acceleration buffer to zero
    void reset_acceleration_buffer() {
        std::fill(std::begin(a), std::end(a), field_zero);
    }

    // This method is called by the integrator to compute accelerations
    virtual void operator() (typename field_container_t::const_iterator x_begin [[maybe_unused]],
                     typename field_container_t::const_iterator x_end [[maybe_unused]],
                     typename field_container_t::const_iterator v_begin [[maybe_unused]],
                     typename field_container_t::iterator a_begin [[maybe_unused]],
                     real_t t) = 0;

    // This method is called from the driver program  to perform one time step of size dt
    void do_step(real_t dt) {
        this->integrator.do_step(dt);
    }

    // Getter for x buffer
    [[nodiscard]] field_container_t const & get_x() const {
        return this->x;
    }

    // Getter for v buffer
    [[nodiscard]] field_container_t const & get_v() const {
        return this->v;
    }

protected:
    const field_value_t field_zero;
    const real_t real_zero;

    index_container_t indices;
    field_container_t x, v, a;
    integrator_t<field_container_t, field_value_t, real_t, generic_system, step_handler_t> integrator;
};

#include "binary_system.h"
#include "unary_system.h"

#endif //INTEGRATORS_SYSTEM_H
