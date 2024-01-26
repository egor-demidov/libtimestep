//
// Created by egor on 1/17/24.
//

#ifndef INTEGRATORS_INTEGRATOR_H
#define INTEGRATORS_INTEGRATOR_H

#include <iterator>
#include <algorithm>
#include <type_traits>

// Base class for all integrators
template <
    typename field_container_t,
    typename field_value_t,
    typename real_t,
    typename functor_t,
    template <
        typename _field_container_t,
        typename _field_value_t>
    typename step_handler_t>
class integrator {
public:

    // Class constructor
    //
    // Notes:
    // Iterators passed to this constructor must remain valid for the duration of use of this integrator
    // x, v, and a buffers must be of the same size
    // The acceleration functor and the step handler must exist for the duration of use of this integrator
    integrator(functor_t & acceleration_functor,                                    // reference to a functor that computes acceleration
               real_t t0,                                                           // reference to a functor that computes acceleration
               typename field_container_t::iterator x_begin,                        // iterator pointing to the start of the x buffer
               typename field_container_t::iterator x_end,                          // iterator pointing to the end of the x buffer
               typename field_container_t::iterator v_begin,                        // iterator pointing to the start of the v buffer
               typename field_container_t::iterator a_begin,                        // iterator pointing to the start of the a buffer
               step_handler_t<field_container_t, field_value_t> & step_handler) :   // reference to an object that handles incrementing positions and velocities

        // Initialize all the member variables
        t(t0), x_begin_itr(x_begin), x_end_itr(x_end), v_begin_itr(v_begin), a_begin_itr(a_begin),
        acceleration_functor(acceleration_functor), step_handler(step_handler) {

        // Check the container type
        static_assert(std::is_same<decltype(*x_begin), field_value_t &>::value,
                "field_container_t must be a container of values of type field_value_t");
    }

protected:
    // This method should be called by any derived class when accelerations need to be recalculated
    void update_acceleration() const {
        // Create const iterators for field arrays that are not supposed to be modified by the acceleration functor
        typename field_container_t::const_iterator x_begin_const_itr = this->x_begin_itr;
        typename field_container_t::const_iterator x_end_const_itr = this->x_end_itr;
        typename field_container_t::const_iterator v_begin_const_itr = this->v_begin_itr;

        // Initialize the accelerations from the initial condition
        this->acceleration_functor(x_begin_const_itr, x_end_const_itr, v_begin_const_itr, this->a_begin_itr, this->t);
    }

    real_t t;
    typename field_container_t::iterator x_begin_itr, x_end_itr, v_begin_itr, a_begin_itr;
    functor_t & acceleration_functor;
    step_handler_t<field_container_t, field_value_t> & step_handler;
};

#include "forward_euler.h"
#include "velocity_verlet_half.h"

#endif //INTEGRATORS_INTEGRATOR_H
