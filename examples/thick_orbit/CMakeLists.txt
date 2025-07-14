# Thick orbit example build configuration

add_executable(thick_orbit_example.x thick_orbit_example.f90)
target_link_libraries(thick_orbit_example.x neo_rt)

# Optional: add as a test that can be run to verify the example works
add_test(NAME thick_orbit_example
         COMMAND thick_orbit_example.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})