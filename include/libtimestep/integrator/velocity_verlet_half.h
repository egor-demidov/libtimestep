//
// Created by egor on 1/18/24.
//

#ifndef INTEGRATORS_VELOCITY_VERLET_HALF_H
#define INTEGRATORS_VELOCITY_VERLET_HALF_H

// Integrator template that implements the Velocity Verlet scheme from the review
template <
    typename field_container_t,
    typename field_value_t,
    typename real_t,
    typename functor_t,
    template <
        typename _field_container_t,
        typename _field_value_t>
    typename step_handler_t>
class velocity_verlet_half : public integrator<field_container_t,field_value_t, real_t, functor_t, step_handler_t> {
public:
    velocity_verlet_half(functor_t & acceleration_functor, real_t t0,
                    typename field_container_t::iterator x_begin,
                    typename field_container_t::iterator x_end,
                    typename field_container_t::iterator v_begin,
                    typename field_container_t::iterator a_begin,
                    step_handler_t<field_container_t, field_value_t> & step_handler) :
            integrator<field_container_t, field_value_t, real_t, functor_t, step_handler_t>(acceleration_functor, t0,
                                                                            x_begin,
                                                                            x_end,
                                                                            v_begin,
                                                                            a_begin,
                                                                            step_handler) {

        // Compute the accelerations
        this->update_acceleration();
    }

    void do_step(real_t dt) override {
        // If this is the first step, take half-a-timestep back in velocity
        if (!velocities_initialized){
            velocities_initialized = true;

            for (size_t n = 0; n < this->x_end_itr - this->x_begin_itr; n ++) {
                field_value_t const & a = *(this->a_begin_itr + n);

                this->step_handler.increment_v(n, -a * dt / 2.0, this->x_begin_itr, this->v_begin_itr);
            }

            this->update_acceleration();
        }

        // Integrate velocity and position
        for (size_t n = 0; n < this->x_end_itr - this->x_begin_itr; n ++) {
            field_value_t const & v = *(this->v_begin_itr + n);
            field_value_t const & a = *(this->a_begin_itr + n);

            this->step_handler.increment_v(n, a*dt, this->x_begin_itr, this->v_begin_itr);
            this->step_handler.increment_x(n, v*dt, this->x_begin_itr, this->v_begin_itr);
        }

        this->update_acceleration();
    }

private:
    bool velocities_initialized = false;
};

#endif //INTEGRATORS_VELOCITY_VERLET_HALF_H
