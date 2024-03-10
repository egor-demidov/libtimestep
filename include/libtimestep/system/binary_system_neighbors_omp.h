//
// Created by egor on 3/4/24.
//

#ifndef INTEGRATORS_BINARY_SYSTEM_NEIGHBORS_OMP_H
#define INTEGRATORS_BINARY_SYSTEM_NEIGHBORS_OMP_H

#include <vector>

#include "system.h"

#include <omp.h>

// This is a base class for a second order system where accelerations depend on binary
// interactions between fields
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
        typename acceleration_handler_t,
        bool have_unary_force>
class binary_system_neighbors_omp : public generic_system<field_value_t, real_t, integrator_t, step_handler_t,
        binary_system_neighbors_omp<field_value_t, real_t, integrator_t, step_handler_t, acceleration_handler_t, have_unary_force>> {
public:
    typedef std::vector<field_value_t> field_container_t;
    typedef std::vector<long> index_container_t;

    binary_system_neighbors_omp(binary_system_neighbors_omp const &) = delete;

    // Class constructor
    //
    // Notes:
    // The acceleration handler and the step handler must exist for the duration of use of this second order system
    // This class itself acts as the acceleration functor for the integrator
    binary_system_neighbors_omp(long n_part,
                    real_t r_verlet,
                    field_container_t x0,                                                 // container with initial positions
                      field_container_t v0,                                                 // container with initial velocities
                      real_t t0,                                                            // integration start time
                      field_value_t field_zero,                                             // zero value of the primary field type used
                      real_t real_zero,                                                     // zero value of the real number type used
                      acceleration_handler_t & acceleration_handler,                        // reference to the object that handles calculating accelerations FOR EACH BINARY INTERACTION
                      step_handler_t<field_container_t, field_value_t> & step_handler) :    // reference to an object that handles incrementing positions and velocities

    // Call the superclass constructor
            generic_system<field_value_t, real_t, integrator_t, step_handler_t, binary_system_neighbors_omp>(std::move(x0),
                                                                                                   std::move(v0), t0, field_zero, real_zero, *this, step_handler),
                                                                                                   n_part(n_part),
                                                                                                    r_verlet(r_verlet),
                                                                                                   acceleration_handler(acceleration_handler),
                                                                                                   neighbor_list(n_part) {}

    // This method is called by the integrator to compute accelerations
    void operator() (typename field_container_t::const_iterator x_begin [[maybe_unused]],
                     typename field_container_t::const_iterator x_end [[maybe_unused]],
                     typename field_container_t::const_iterator v_begin [[maybe_unused]],
                     typename field_container_t::iterator a_begin [[maybe_unused]],
                     real_t t) {

        this->reset_acceleration_buffer();

#pragma omp parallel for default(none) shared(t)
        for (long i = 0; i < n_part; i ++) {
            for (long j : neighbor_list[i]) {
                if (i == j) [[unlikely]]
                    continue;

                this->a[i] += acceleration_handler.compute_acceleration(i, j, this->get_x(), this->get_v(), t);
            }

            if constexpr (have_unary_force) {
                this->a[i] += acceleration_handler.compute_acceleration(i, this->get_x(), this->get_v(), t);
            }
        }
    }

    // This method is called by the driver periodically to update the neighbor lists
    void update_neighbor_list() {
#pragma omp parallel for default(none)
        for (long i = 0; i < n_part; i ++) {
            neighbor_list[i].clear();

            for (long j = 0; j < n_part; j ++) {
                if (i == j)
                    continue;

                real_t distance = (this->get_x()[i] - this->get_x()[j]).norm();
                if (distance < r_verlet)
                    neighbor_list[i].emplace_back(j);
            }
        }
    }

private:
    const long n_part;
    const double r_verlet;
    acceleration_handler_t & acceleration_handler;
    std::vector<std::vector<long>> neighbor_list;
};

#endif //INTEGRATORS_BINARY_SYSTEM_NEIGHBORS_OMP_H
