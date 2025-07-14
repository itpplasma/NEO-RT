file(GLOB TEST_SOURCES "test_*.f90")

# Exclude POTATO-specific tests when thick orbits are disabled
if(NOT USE_THICK_ORBITS)
    list(REMOVE_ITEM TEST_SOURCES 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_build.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_field_bridge.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_field_init.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_spline_real.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_find_bounce_real.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_physics_comparison.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_real_integration.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_comparison_plots.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test_potato_integration_fix.f90)
endif()

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name}.x ${test_src})
    target_link_libraries(${test_name}.x neo_rt util_for_test)
    
    # Add POTATO includes for POTATO-specific tests when thick orbits enabled
    if(USE_THICK_ORBITS)
        string(FIND ${test_src} "test_potato" test_is_potato)
        if(NOT ${test_is_potato} EQUAL -1)
            target_include_directories(${test_name}.x PRIVATE "${PROJECT_BINARY_DIR}/POTATO")
            add_dependencies(${test_name}.x potato_base)
        endif()
    endif()
endforeach()

add_test(NAME test_timestep
         COMMAND test_timestep.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Test for POTATO frequencies (now always available with runtime dispatch)
add_test(NAME test_potato_frequencies
         COMMAND test_potato_frequencies.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Test for POTATO library build integration (only when thick orbits enabled)
if(USE_THICK_ORBITS)
    add_test(NAME test_potato_build
             COMMAND test_potato_build.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    
    # Test for POTATO field bridge
    add_test(NAME test_potato_field_bridge
             COMMAND test_potato_field_bridge.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for POTATO field initialization
    add_test(NAME test_potato_field_init
             COMMAND test_potato_field_init.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for real POTATO spline implementation
    add_test(NAME test_potato_spline_real
             COMMAND test_potato_spline_real.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for real POTATO find_bounce implementation
    add_test(NAME test_potato_find_bounce_real
             COMMAND test_potato_find_bounce_real.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for real vs stub POTATO physics comparison
    add_test(NAME test_potato_physics_comparison
             COMMAND test_potato_physics_comparison.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for real POTATO integration
    add_test(NAME test_potato_real_integration
             COMMAND test_potato_real_integration.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for POTATO comparison plots generation
    add_test(NAME test_potato_comparison_plots
             COMMAND test_potato_comparison_plots.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for POTATO EFIT integration (Phase F.1.1)
    add_test(NAME test_potato_efit_integration
             COMMAND test_potato_efit_integration.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for EFIT-Boozer coordinate consistency (Phase F.1.2)
    add_test(NAME test_efit_boozer_consistency
             COMMAND test_efit_boozer_consistency.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for EFIT vs Boozer magnetic field comparison
    add_test(NAME test_efit_boozer_field_comparison
             COMMAND test_efit_boozer_field_comparison.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Simple test for EFIT vs Boozer equilibrium verification
    add_test(NAME test_efit_boozer_simple_comparison
             COMMAND test_efit_boozer_simple_comparison.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    
    # Test for fixing POTATO integration issues (Phase G.1.1)
    add_test(NAME test_potato_integration_fix
             COMMAND test_potato_integration_fix.x
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()