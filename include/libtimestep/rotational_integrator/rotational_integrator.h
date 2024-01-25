//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_ROTATIONAL_INTEGRATOR_H
#define INTEGRATORS_ROTATIONAL_INTEGRATOR_H

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
class rotational_integrator {
public:
    rotational_integrator(functor_t & acceleration_functor, real_t t0,
               typename field_container_t::iterator x_begin,
               typename field_container_t::iterator x_end,
               typename field_container_t::iterator v_begin,
               typename field_container_t::iterator a_begin,
               typename field_container_t::iterator theta_begin,
               typename field_container_t::iterator omega_begin,
               typename field_container_t::iterator alpha_begin,
               step_handler_t<field_container_t, field_value_t> & step_handler) :
            t(t0), x_begin_itr(x_begin), x_end_itr(x_end), v_begin_itr(v_begin), a_begin_itr(a_begin),
            theta_begin_itr(theta_begin), omega_begin_itr(omega_begin),
            alpha_begin_itr(alpha_begin), acceleration_functor(acceleration_functor), step_handler(step_handler) {

        static_assert(std::is_same<decltype(*x_begin), field_value_t &>::value,
                      "field_container_t must be a container of values of type field_value_t");
    }

    virtual void do_step(real_t dt) = 0;

protected:
    void update_acceleration() const {
        // Create const iterators for field arrays that are not supposed to be modified by the acceleration functor
        typename field_container_t::const_iterator x_begin_const_itr = this->x_begin_itr;
        typename field_container_t::const_iterator x_end_const_itr = this->x_end_itr;
        typename field_container_t::const_iterator v_begin_const_itr = this->v_begin_itr;
        typename field_container_t::const_iterator theta_begin_const_itr = this->theta_begin_itr;
        typename field_container_t::const_iterator omega_begin_const_itr = this->omega_begin_itr;

        // Initialize the accelerations from the initial condition
        this->acceleration_functor(x_begin_const_itr, x_end_const_itr, v_begin_const_itr, this->a_begin_itr,
                                   theta_begin_const_itr, omega_begin_const_itr, this->alpha_begin_itr, this->t);
    }

    real_t t;
    typename field_container_t::iterator x_begin_itr, x_end_itr, v_begin_itr, a_begin_itr,
                            theta_begin_itr, omega_begin_itr, alpha_begin_itr;
    functor_t & acceleration_functor;
    step_handler_t<field_container_t, field_value_t> & step_handler;
};

#include "rotational_forward_euler.h"
#include "rotational_velocity_verlet_half.h"

#endif //INTEGRATORS_ROTATIONAL_INTEGRATOR_H
