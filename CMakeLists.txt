cmake_minimum_required(VERSION 3.16)
project(integrators)

include(CTest)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (${CMAKE_COMPILER_IS_GNUCXX})
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -flto=auto ${CMAKE_CXX_FLAGS}")
endif ()

include_directories(test/deps/eigen)

#add_subdirectory(test/deps/matplotplusplus)
set(CMAKE_FIND_DEBUG_MODE ON)

# GNU libstdc++ uses libtbb to implement C++ 17 parallel functionality
if (${CMAKE_COMPILER_IS_GNUCXX})
    find_package(TBB REQUIRED)
endif ()

include_directories(include)

#add_executable(integrators main.cpp)
#add_executable(particle_dynamics test/particle_dynamics.cpp test/write_vtk.cpp test/compute_energy.cpp)
#add_executable(system_test test/system_test.cpp test/write_vtk.cpp test/compute_energy.cpp)
add_executable(rotational_system_test test/rotational_system_test.cpp test/write_vtk.cpp ../soot-indentation/energy/compute_energy.cpp)
add_executable(oscillator_test test/oscillator.cpp)
add_executable(particle_dynamics_test test/particle_dynamics.cpp test/write_vtk.cpp ../soot-indentation/energy/compute_energy.cpp)

if (${CMAKE_COMPILER_IS_GNUCXX})
    #target_link_libraries(particle_dynamics PRIVATE TBB::tbb)
    #target_link_libraries(system_test PRIVATE TBB::tbb)
    target_link_libraries(rotational_system_test PRIVATE TBB::tbb)
    target_link_libraries(oscillator_test PRIVATE TBB::tbb)
    target_link_libraries(particle_dynamics_test PRIVATE TBB::tbb)
    #target_link_libraries(integrators matplot)
endif ()

add_test(NAME oscillator_test COMMAND ${CMAKE_BINARY_DIR}/oscillator_test)
add_test(NAME particle_dynamics_test COMMAND ${CMAKE_BINARY_DIR}/particle_dynamics_test)
