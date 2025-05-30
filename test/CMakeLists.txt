file(GLOB TEST_SOURCES "test_*.f90")
foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name}.x ${test_src})
    target_link_libraries(${test_name}.x neo_rt util_for_test)
endforeach()

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/in_file
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${CMAKE_SOURCE_DIR}/examples/circ.bc
            ${CMAKE_CURRENT_BINARY_DIR}/in_file
    DEPENDS ${CMAKE_SOURCE_DIR}/examples/circ.bc
    COMMENT "Creating symlink to test.bc as in_file"
)

add_custom_target(link_test_bc ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/in_file
)

add_test(NAME test_timestep
         COMMAND test_timestep.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
