//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_BINARY_SYSTEM_H
#define INTEGRATORS_BINARY_SYSTEM_H

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

    // The containers cannot be resized after the integrator has been instantiated because that would
    // invalidate the iterators
    binary_system(field_container_t x0, field_container_t v0, real_t t0, field_value_t field_zero, real_t real_zero,
                  acceleration_handler_t & acceleration_handler, step_handler_t<field_container_t, field_value_t> & step_handler) :
            generic_system<field_value_t, real_t, integrator_t, step_handler_t, binary_system>(std::move(x0),
                                                                                std::move(v0),
                                                                                t0,
                                                                                field_zero,
                                                                                real_zero, *this, step_handler), acceleration_handler(acceleration_handler) {}

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
        });
    }

private:
    acceleration_handler_t & acceleration_handler;
};

#endif //INTEGRATORS_BINARY_SYSTEM_H
