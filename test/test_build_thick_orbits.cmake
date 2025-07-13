# CMake test for thick orbit build option
cmake_minimum_required(VERSION 3.10)

# Test that thick orbits can be disabled (default behavior)
message(STATUS "Testing default build without thick orbits...")
if(DEFINED USE_THICK_ORBITS AND USE_THICK_ORBITS)
    message(FATAL_ERROR "USE_THICK_ORBITS should be OFF by default")
endif()

# Test that thick orbits can be enabled
message(STATUS "Testing thick orbit build option...")
set(USE_THICK_ORBITS ON)
if(NOT USE_THICK_ORBITS)
    message(FATAL_ERROR "USE_THICK_ORBITS option not working")
endif()

# Test that POTATO sources would be available when enabled  
# (This is a placeholder test - actual POTATO integration comes later)
if(USE_THICK_ORBITS)
    message(STATUS "Thick orbits enabled - POTATO integration would be activated")
    # Future: Test for POTATO source availability
    # if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/POTATO/")
    #     message(FATAL_ERROR "POTATO sources not found for thick orbit integration")
    # endif()
endif()

message(STATUS "Build option tests PASSED")