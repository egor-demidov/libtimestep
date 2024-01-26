//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_SYSTEM_H
#define INTEGRATORS_SYSTEM_H

#include <vector>
#include <execution>

#include "../exception/exception.h"

// This is a base class for a second order system
//
// Template parameters:
// primary field type
// real numer type
// integrator type
// step handler type
// step handler type
// acceleration functor type
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
        typename functor_t>
class generic_system {
public:
    typedef std::vector<field_value_t> field_container_t;
    typedef std::vector<size_t> index_container_t;

    // Class constructor
    //
    // Notes:
    // The acceleration functor and the step handler must exist for the duration of use of this second order system
    generic_system(field_container_t x0,                                                // container with initial positions
                   field_container_t v0,                                                // container with initial velocities
                   real_t t0,                                                           // integration start time
                   field_value_t field_zero,                                            // zero value of the primary field type used
                   real_t real_zero,                                                    // zero value of the real number type used
                   functor_t & functor_ref,                                             // reference to the object that handles calculating accelerations
                   step_handler_t<field_container_t, field_value_t> & step_handler) :   // reference to an object that handles incrementing positions and velocities

        // Initialize the member variables
        field_zero(std::move(field_zero)), real_zero(real_zero), x(std::move(x0)), v(std::move(v0)), a(x.size()),
        integrator(functor_ref, t0, std::begin(this->x), std::end(this->x),
                   std::begin(this->v), std::begin(this->a), step_handler) {

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

    // This method is called from the driver program to perform one time step of size dt
    void do_step(real_t dt /* time step */) {
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
    integrator_t<field_container_t, field_value_t, real_t, functor_t, step_handler_t> integrator;
};

#include "binary_system.h"
#include "unary_system.h"

#endif //INTEGRATORS_SYSTEM_H
