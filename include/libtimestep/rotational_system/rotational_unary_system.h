//
// Created by egor on 1/21/24.
//

#ifndef INTEGRATORS_ROTATIONAL_UNARY_SYSTEM_H
#define INTEGRATORS_ROTATIONAL_UNARY_SYSTEM_H

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
class rotational_unary_system : public rotational_generic_system<field_value_t, real_t, integrator_t, step_handler_t, rotational_unary_system<field_value_t, real_t, integrator_t, step_handler_t>> {
public:
    typedef std::vector<field_value_t> field_container_t;
    typedef std::vector<size_t> index_container_t;

    rotational_unary_system(field_container_t x0, field_container_t v0,
                             field_container_t theta0, field_container_t omega0,
                             real_t t0, field_value_t field_zero, real_t real_zero) :
            rotational_generic_system<field_value_t, real_t, integrator_t, step_handler_t, rotational_unary_system>(std::move(x0),
                                                                                           std::move(v0),
                                                                                           std::move(theta0),
                                                                                           std::move(omega0),
                                                                                           t0, field_zero, real_zero, *this) {}

    // This method should compute the acceleration of particle i
    virtual std::pair<field_value_t, field_value_t> compute_accelerations(size_t i, real_t t) = 0;

    // This method is called by the integrator to compute accelerations
    void operator() (typename field_container_t::const_iterator x_begin [[maybe_unused]],
                     typename field_container_t::const_iterator x_end [[maybe_unused]],
                     typename field_container_t::const_iterator v_begin [[maybe_unused]],
                     typename field_container_t::iterator a_begin [[maybe_unused]],
                     typename field_container_t::const_iterator theta_begin [[maybe_unused]],
                     typename field_container_t::const_iterator omega_begin [[maybe_unused]],
                     typename field_container_t::iterator alpha_begin [[maybe_unused]],
                     real_t t) {

        this->reset_acceleration_buffers();

        std::for_each(this->indices.begin(), this->indices.end(), [t, this] (size_t i) {
            auto [a_i, alpha_i] = this->compute_acceleration(i, t);
            this->a[i] += a_i;
            this->alpha[i] += alpha_i;
        });
    }
};

#endif //INTEGRATORS_ROTATIONAL_UNARY_SYSTEM_H
