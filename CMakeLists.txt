cmake_minimum_required(VERSION 3.26)
project(integrators)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -flto=auto")

include_directories(deps/eigen)
add_subdirectory(deps/matplotplusplus)

find_package(TBB REQUIRED)

include_directories(include)

add_executable(integrators main.cpp)
add_executable(particle_dynamics test/particle_dynamics.cpp test/write_vtk.cpp test/compute_energy.cpp)
add_executable(system_test test/system_test.cpp test/write_vtk.cpp test/compute_energy.cpp)
add_executable(rotational_system_test test/rotational_system_test.cpp test/write_vtk.cpp test/compute_energy.cpp)

target_link_libraries(particle_dynamics PRIVATE TBB::tbb)
target_link_libraries(system_test PRIVATE TBB::tbb)
target_link_libraries(rotational_system_test PRIVATE TBB::tbb)
target_link_libraries(integrators matplot)
