//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_BINARY_SYSTEM_H
#define INTEGRATORS_BINARY_SYSTEM_H

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
    typename acceleration_handler_t>
class binary_system : public generic_system<field_value_t, real_t, integrator_t, step_handler_t,
        binary_system<field_value_t, real_t, integrator_t, step_handler_t, acceleration_handler_t>> {
public:
    typedef std::vector<field_value_t> field_container_t;
    typedef std::vector<size_t> index_container_t;

    // Class constructor
    //
    // Notes:
    // The acceleration handler and the step handler must exist for the duration of use of this second order system
    // This class itself acts as the acceleration functor for the integrator
    binary_system(field_container_t x0,                                                 // container with initial positions
                  field_container_t v0,                                                 // container with initial velocities
                  real_t t0,                                                            // integration start time
                  field_value_t field_zero,                                             // zero value of the primary field type used
                  real_t real_zero,                                                     // zero value of the real number type used
                  acceleration_handler_t & acceleration_handler,                        // reference to the object that handles calculating accelerations FOR EACH BINARY INTERACTION
                  step_handler_t<field_container_t, field_value_t> & step_handler) :    // reference to an object that handles incrementing positions and velocities

        // Call the superclass constructor
        generic_system<field_value_t, real_t, integrator_t, step_handler_t, binary_system>(std::move(x0),
            std::move(v0), t0, field_zero, real_zero, *this, step_handler), acceleration_handler(acceleration_handler) {}

    // This method is called by the integrator to compute accelerations
    void operator() (typename field_container_t::const_iterator x_begin [[maybe_unused]],
                     typename field_container_t::const_iterator x_end [[maybe_unused]],
                     typename field_container_t::const_iterator v_begin [[maybe_unused]],
                     typename field_container_t::iterator a_begin [[maybe_unused]],
                     real_t t) {

        this->reset_acceleration_buffer();

        std::for_each(std::execution::par_unseq, this->indices.begin(), this->indices.end(), [t, this] (size_t i) {
            std::for_each(this->indices.begin(), this->indices.end(), [t, i, this] (size_t j) {
                if (i == j)
                    return;

                this->a[i] += acceleration_handler.compute_acceleration(i, j, this->get_x(), this->get_v(), t);
            });
            this->a[i] += acceleration_handler.compute_acceleration(i, this->get_x(), this->get_v(), t);
        });
    }

private:
    acceleration_handler_t & acceleration_handler;
};

#endif //INTEGRATORS_BINARY_SYSTEM_H
