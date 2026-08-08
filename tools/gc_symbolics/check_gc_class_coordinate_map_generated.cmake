if(NOT DEFINED NEORT_SOURCE_DIR)
    message(FATAL_ERROR "NEORT_SOURCE_DIR is required")
endif()

if(DEFINED ENV{NEORT_TEST_TMPDIR})
    set(regeneration_dir "$ENV{NEORT_TEST_TMPDIR}/neort-gc-class-coordinate-map")
else()
    set(regeneration_dir "${NEORT_SOURCE_DIR}/.neort-gc-class-coordinate-map")
endif()

file(REMOVE_RECURSE "${regeneration_dir}")
file(MAKE_DIRECTORY "${regeneration_dir}")
execute_process(
    COMMAND fo exec gen_gc_class_coordinate_map
        "${regeneration_dir}/neort_gc_class_coordinate_map_symbolic.f90"
        "${regeneration_dir}/neort_gc_class_coordinate_map_inventory.txt"
    WORKING_DIRECTORY "${NEORT_SOURCE_DIR}/tools/gc_symbolics"
    RESULT_VARIABLE generation_status
)
if(NOT generation_status EQUAL 0)
    message(FATAL_ERROR "Fortsym class-map regeneration failed: ${generation_status}")
endif()

foreach(generated_name
        neort_gc_class_coordinate_map_symbolic.f90
        neort_gc_class_coordinate_map_inventory.txt)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E compare_files
            "${regeneration_dir}/${generated_name}"
            "${NEORT_SOURCE_DIR}/src/generated/${generated_name}"
        RESULT_VARIABLE compare_status
    )
    if(NOT compare_status EQUAL 0)
        message(FATAL_ERROR "Generated class-map file differs: ${generated_name}")
    endif()
endforeach()
message(STATUS "Fortsym POTATO class-coordinate outputs match")
