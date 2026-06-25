file(GLOB TEST_SOURCES "test_*.f90")

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name}.x ${test_src})
    target_link_libraries(${test_name}.x neo_rt "util_for_test_${PROJECT_NAME}")
endforeach()

add_test(NAME test_timestep
    COMMAND test_timestep.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_parallel
    COMMAND ${CMAKE_BINARY_DIR}/test_parallel.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

add_test(NAME test_neort_lib
    COMMAND ${CMAKE_BINARY_DIR}/test_neort_lib.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

# Generate a reference Boozer chartmap from circ.bc for test_chartmap_input.
#
# LIBNEO_EQDSK_BOOZ_PYTHON: path to the eqdsk-booz worktree python/ directory
#   (contains libneo/bc_to_booz_xform.py and libneo/booz_xform_to_boozer_chartmap.py).
#   Pass as -DLIBNEO_EQDSK_BOOZ_PYTHON=... at configure time, or set the
#   NEORT_LIBNEO_CONVERTERS env var at ctest time.
# LIBNEO_PYTHON: path to the libneo repo python/ directory (boozer_chartmap_writer etc.).
#   Defaults to the fetched libneo source tree python/ directory.
set(CIRC_CHARTMAP_NC "${CMAKE_BINARY_DIR}/circ_chartmap.nc")
set(CIRC_BC "${CMAKE_SOURCE_DIR}/examples/circ.bc")

if(NOT LIBNEO_EQDSK_BOOZ_PYTHON)
    set(LIBNEO_EQDSK_BOOZ_PYTHON "$ENV{NEORT_LIBNEO_CONVERTERS}")
endif()
if(NOT LIBNEO_PYTHON AND libneo_SOURCE_DIR)
    set(LIBNEO_PYTHON "${libneo_SOURCE_DIR}/python")
endif()
if(NOT LIBNEO_PYTHON)
    set(LIBNEO_PYTHON "$ENV{PYTHONPATH}")
endif()

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_test(NAME gen_circ_chartmap
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/gen_circ_chartmap.py"
            "${CIRC_CHARTMAP_NC}"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    if(LIBNEO_EQDSK_BOOZ_PYTHON)
        set_tests_properties(gen_circ_chartmap PROPERTIES
            FIXTURES_SETUP chartmap_fixture
            ENVIRONMENT
                "NEORT_LIBNEO_CONVERTERS=${LIBNEO_EQDSK_BOOZ_PYTHON};PYTHONPATH=${LIBNEO_PYTHON}:$ENV{PYTHONPATH}"
        )
    else()
        # LIBNEO_EQDSK_BOOZ_PYTHON not configured; rely on NEORT_LIBNEO_CONVERTERS
        # being set in the environment at ctest time (or set with -DLIBNEO_EQDSK_BOOZ_PYTHON=).
        set_tests_properties(gen_circ_chartmap PROPERTIES
            FIXTURES_SETUP chartmap_fixture
            ENVIRONMENT "PYTHONPATH=${LIBNEO_PYTHON}:$ENV{PYTHONPATH}"
        )
    endif()

    add_test(NAME test_chartmap_input
        COMMAND ${CMAKE_BINARY_DIR}/test_chartmap_input.x
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
    )
    set_tests_properties(test_chartmap_input PROPERTIES
        FIXTURES_REQUIRED chartmap_fixture
        ENVIRONMENT
            "CHARTMAP_FILE=${CIRC_CHARTMAP_NC};BC_FILE=${CIRC_BC}"
    )
else()
    add_test(NAME test_chartmap_input
        COMMAND ${CMAKE_BINARY_DIR}/test_chartmap_input.x
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
    )
endif()
