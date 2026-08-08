if(NOT DEFINED NEORT_SOURCE_DIR)
    message(FATAL_ERROR "NEORT_SOURCE_DIR is required")
endif()

if(DEFINED ENV{NEORT_TEST_TMPDIR})
    set(regeneration_dir "$ENV{NEORT_TEST_TMPDIR}/neort-circular-continuation-regeneration")
else()
    set(regeneration_dir "${NEORT_SOURCE_DIR}/.neort-circular-continuation-regeneration")
endif()

file(REMOVE_RECURSE "${regeneration_dir}")
file(MAKE_DIRECTORY "${regeneration_dir}")
execute_process(
    COMMAND fo exec gen_circular_flux_continuation
        "${regeneration_dir}/neort_circular_flux_continuation_symbolic.f90"
        "${regeneration_dir}/neort_circular_flux_continuation_limit_symbolic.f90"
        "${regeneration_dir}/circular_flux_continuation_generated.py"
    WORKING_DIRECTORY "${NEORT_SOURCE_DIR}/tools/gc_symbolics"
    RESULT_VARIABLE generation_status
)
if(NOT generation_status EQUAL 0)
    message(FATAL_ERROR "Fortsym circular continuation regeneration failed: ${generation_status}")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files
        "${regeneration_dir}/neort_circular_flux_continuation_symbolic.f90"
        "${NEORT_SOURCE_DIR}/src/generated/neort_circular_flux_continuation_symbolic.f90"
    RESULT_VARIABLE kernel_status
)
if(NOT kernel_status EQUAL 0)
    message(FATAL_ERROR "Generated circular continuation kernel differs")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files
        "${regeneration_dir}/neort_circular_flux_continuation_limit_symbolic.f90"
        "${NEORT_SOURCE_DIR}/src/generated/neort_circular_flux_continuation_limit_symbolic.f90"
    RESULT_VARIABLE limit_status
)
if(NOT limit_status EQUAL 0)
    message(FATAL_ERROR "Generated circular continuation limit kernel differs")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files
        "${regeneration_dir}/circular_flux_continuation_generated.py"
        "${NEORT_SOURCE_DIR}/POTATO/test/golden_record_resonance/circular_flux_continuation_generated.py"
    RESULT_VARIABLE table_status
)
if(NOT table_status EQUAL 0)
    message(FATAL_ERROR "Generated circular continuation table differs")
endif()
message(STATUS "Fortsym circular continuation outputs match")
