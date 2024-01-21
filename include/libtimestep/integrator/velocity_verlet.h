//
// Created by egor on 1/18/24.
//

#ifndef INTEGRATORS_VELOCITY_VERLET_H
#define INTEGRATORS_VELOCITY_VERLET_H

// Integrator template that implements the Velocity Verlet integration scheme
template <typename field_container_t, typename field_value_t, typename real_t, typename functor_t>
class velocity_verlet : public integrator<field_container_t, field_value_t, real_t, functor_t> {
public:
    velocity_verlet(functor_t & acceleration_functor, real_t t0,
                    typename field_container_t::iterator x_begin,
                    typename field_container_t::iterator x_end,
                    typename field_container_t::iterator v_begin,
                    typename field_container_t::iterator a_begin,
                    typename field_container_t::iterator a_copy_begin) :
            integrator<field_container_t, field_value_t, real_t, functor_t>(acceleration_functor, t0,
                                                                            x_begin,
                                                                            x_end,
                                                                            v_begin,
                                                                            a_begin), a_copy_begin_itr(a_copy_begin) {

        // Compute the accelerations
        this->update_acceleration();
    }

    void do_step(real_t dt) override {

        // Increment time
        this->t += dt;

        // Integrate position
        for (size_t n = 0; n < this->x_end_itr - this->x_begin_itr; n ++) {
            field_value_t & x = *(this->x_begin_itr + n);
            field_value_t const & v = *(this->v_begin_itr + n);
            field_value_t const & a = *(this->a_begin_itr + n);

            x += v*dt + 0.5*a*dt*dt;
        }

        this->swap_acceleration_buffers();
        this->update_acceleration();

        // Integrate velocity
        for (size_t n = 0; n < this->x_end_itr - this->x_begin_itr; n ++) {
            field_value_t & v = *(this->v_begin_itr + n);
            field_value_t const & a_new = *(this->a_begin_itr + n);
            field_value_t const & a_old = *(this->a_copy_begin_itr + n);

            v += 0.5*(a_old + a_new)*dt;
        }
    }

private:
    void swap_acceleration_buffers() {
        auto temp_iterator_holder = this->a_begin_itr;
        this->a_begin_itr = this->a_copy_begin_itr;
        this->a_copy_begin_itr = this->a_begin_itr;
    }

    typename field_container_t::iterator a_copy_begin_itr;
};

#endif //INTEGRATORS_VELOCITY_VERLET_H
