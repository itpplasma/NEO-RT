# Thick orbit example build configuration

add_executable(thick_orbit_example.x thick_orbit_example.f90)
target_link_libraries(thick_orbit_example.x neo_rt)

# Canonical frequency plotting example with fortplotlib
add_executable(plot_canonical_frequencies.x plot_canonical_frequencies.f90)
target_link_libraries(plot_canonical_frequencies.x neo_rt fortplotlib::fortplotlib)

# Simple frequency plotting example
add_executable(plot_frequencies_simple.x plot_frequencies_simple.f90)
target_link_libraries(plot_frequencies_simple.x neo_rt fortplotlib::fortplotlib)

# Demo plotting example with synthetic data
add_executable(plot_canonical_frequencies_demo.x plot_canonical_frequencies_demo.f90)
target_link_libraries(plot_canonical_frequencies_demo.x neo_rt fortplotlib::fortplotlib)

# Optional: add as a test that can be run to verify the example works
add_test(NAME thick_orbit_example
         COMMAND thick_orbit_example.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})