if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    cmake_minimum_required(VERSION 3.22)
    project(neort_symbolics LANGUAGES Fortran)
endif()

if(NOT TARGET fortsym)
    if(FORTSYM_SOURCE_DIR)
        add_subdirectory("${FORTSYM_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/fortsym")
    else()
        find_package(fortsym CONFIG REQUIRED)
    endif()
endif()

add_executable(gen_neort_separatrix_frequency
    gen_neort_separatrix_frequency.f90)
target_link_libraries(gen_neort_separatrix_frequency PRIVATE fortsym)
