# libtimestep

This is a C++ template library that implements various time-integration schemes for
Ordinary Differential Equations (ODEs).

### Installation

This is a header-only library that does not require compiling an object an linking against it.
If you are using CMake to manage your C++ project, you can include libtimestep by adding
the following lines to your CMakeLists.txt file:

```cmake
include_directories(/path/to/libtimestep)
```

If your project is a git repository, libtimestep can be added as a submodule.
The following commands can be executed in your shell to create a `deps` folder in
your project and place the libtimestep submodule in it:

```shell
mkdir deps
git submodule add https://github.com/eg0000r-pub/libtimestep.git deps/libtimestep
```

Then, every time someone clones your submodule-containing project, they need to use:

```shell
git clone --recurse-submodules <your git repo>
```

### Implemented integration schemes

##### Forward Euler

The following recurrence relation is used to integrate a second-order system:
$$x_{t+\Delta t}=x_t+v_t\Delta t+\frac{1}{2}a(x_t,v_t)\Delta t^2$$
$$v_{t+\Delta t}=v_t + a(x_t,v_t)t\Delta t$$

##### Velocity Verlet

The following recurrence relation is used to integrate a second-order system:
$$v_{t+\Delta t/2}=v_{t-\Delta t/2}+a(x_t,v_{t-\Delta t / 2})\Delta t$$
$$x_{t+\Delta t}=x_t+v_{t+\Delta t/2}\Delta t$$

### Usage

This is an example where Forward Euler integration scheme to solve the damped
oscillator equation:
$$m\frac{d^2x}{dt^2}+\gamma\frac{dx}{dt}+kx=g(t)$$
where $m$ is mass, $k$ is stiffness, $\gamma$ is the damping coefficient, and $g(t)$ is the external
force applied to the system. We will be assuming that the system is initially at rest and
that the oscillator is critically damped, i.e. $\gamma=2\sqrt{mk}$. The external force is given by
$g(t)=u(t)$ where $u(t)$ is Heaviside step function. The analytical solution for this problem is then:
$$x(t)=\frac{1}{k}\left[1-e^{-\omega_0 t}\left(\omega_0t+1\right)\right]$$
where
$$\omega_o=\sqrt{\frac{k}{m}}$$

We will use `std::vector<double>` to store `x` and `v` fields. Any container that implements
`iterator` and `const_iterator` works. The reason why we are not using scalars to store fields
is because libtimestep interfaces can be used to solve systems of second-order ODEs. But
if you are not dealing with a system and need to solve a single ODE, then vectors with single values
can be used to represent that system.

The mathematical system is represented by a functor, i.e. an object that overloads 
`operator ()` and can behave as a function, but also has a constructor and member variables.
The functor that can represent our system can look like this:
```c++
class System {
public:
    System(double m, double k, double gamma_d) :
        m(m), k(k), gamma_d(gamma_d) {}
    
    void oprator (std::vector<double>::const_iterator x_begin,
                  std::vector<double>::const_iterator x_end,
                  std::vector<double>::const_iterator v_begin,
                  std::vector<double>::iterator a_begin,
                  double t) const {
        
        // This is where we compute the accelerations and place them in the
        //std:: vector represented by a_begin iterator
        
        // Since here we are solving a single equation. we do not need to loop over
        // std::vector elements
        
        // Create references for readability
        double const & x = *x_begin;
        double const & v = *v_begin;
        double & a = *a_begin;
        
        // Calculate the acceleration
        a = 1.0 / this->m * (1.0 - this->k * x - this->gamma_d * v);
    }
        
private;
    const double m, k, gamma_d;
};
```

If you are dealing with a multi-dimensional problem, i.e. 3D, and are using a linear
algebra library to represent vectors, i.e. libeigen's `Eigen::Vector3d`, then you could use
`std::vector<Eigen::Vector3d>` to represent fields of multi-dimensional systems.

But in this example we are solving a 1D problem. In the main function, let's initialize our fields:
```c++
std::vector<double> x = {0.0};
std::vector<double> v = {0.0};
std::vector<double> a = {0.0};
```
The values that we set our fields to will be the initial values for integration.
We also initialized a container for accelerations. While the value that we stored in 
it will not be used in the simulation, it is a buffer that will be used by the integrator
internally. Then, let's instantiate the system and the integrator:
```c++
System system(m, k, gamma_d); // System object must exist for the storage duration of forward_euler
forward_euler<std::vector<double>, double, double, System> integrator(system, 0.0, x.begin(), x.end(), v.begin(), a.begin());
```
Now we are ready to do the time-stepping and solve the system numerically. Since we want to store the solution at every
time step, let us create buffers to store `t` and `x` values.
```c++
const double dt = 0.01;
const double t_tot = 5.0;
const auto num_steps = size_t(t_tot / dt);

std::vector<double> t_buffer;
std::vector<double> x_buffer;

double t = 0.0;
for (size_t n = 0; n < num_steps; n ++) {
    t_buffer.emplace_back(t);
    x_buffer.emplace_back(x[0]);
    integrator.do_step(dt);
    t += dt;
}
```
