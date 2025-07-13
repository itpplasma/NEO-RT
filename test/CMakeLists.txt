file(GLOB TEST_SOURCES "test_*.f90")
foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name}.x ${test_src})
    target_link_libraries(${test_name}.x neo_rt util_for_test)
endforeach()

add_test(NAME test_timestep
         COMMAND test_timestep.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Test for POTATO frequencies (now always available with runtime dispatch)
add_test(NAME test_potato_frequencies
         COMMAND test_potato_frequencies.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})