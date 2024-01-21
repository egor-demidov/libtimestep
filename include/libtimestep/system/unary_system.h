//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_UNARY_SYSTEM_H
#define INTEGRATORS_UNARY_SYSTEM_H

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
class unary_system : public generic_system<field_value_t, real_t, integrator_t, step_handler_t> {
public:
    typedef std::vector<field_value_t> field_container_t;
    typedef std::vector<size_t> index_container_t;

    // The containers cannot be resized after the integrator has been instantiated because that would
    // invalidate the iterators
    unary_system(field_container_t x0, field_container_t v0, real_t t0, field_value_t field_zero, real_t real_zero) :
            generic_system<field_value_t, real_t, integrator_t, step_handler_t>(std::move(x0),
                                                                                std::move(v0),
                                                                                t0,
                                                                                field_zero,
                                                                                real_zero) {
    }

    // This method should compute the acceleration of particle i
    virtual field_value_t compute_acceleration(size_t i, real_t t) = 0;

    // This method is called by the integrator to compute accelerations
    void operator() (typename field_container_t::const_iterator x_begin [[maybe_unused]],
                     typename field_container_t::const_iterator x_end [[maybe_unused]],
                     typename field_container_t::const_iterator v_begin [[maybe_unused]],
                     typename field_container_t::iterator a_begin [[maybe_unused]],
                     real_t t) override {

        this->reset_acceleration_buffer();

        std::for_each(this->indices.begin(), this->indices.end(), [t, this] (size_t i) {
            this->a[i] += this->compute_acceleration(i, t);
        });
    }
};

#endif //INTEGRATORS_UNARY_SYSTEM_H
