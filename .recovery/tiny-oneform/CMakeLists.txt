cmake_minimum_required(VERSION 3.20)
project(neort_raw_oneform_verifier LANGUAGES Fortran)
enable_testing()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
set(CACHE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/cache)
add_library(raw_probe STATIC
    src/perturbation_field.f90
    src/perturbation_field_netcdf.f90
    src/raw_drive.f90)
target_include_directories(raw_probe PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY} ${CACHE_ROOT}/include
    ${CACHE_ROOT}/fortio-modules)
target_compile_options(raw_probe PRIVATE
    -g -O0 -fcheck=all -fbacktrace -Wall -Wextra -ffree-line-length-none)
target_link_libraries(raw_probe PUBLIC
    ${CACHE_ROOT}/libinterpolate.a ${CACHE_ROOT}/libfield_base.a
    ${CACHE_ROOT}/libneo_util.a ${CACHE_ROOT}/libfortio.a)
add_executable(test_raw_oneform.x src/test_raw_oneform.f90)
target_compile_options(test_raw_oneform.x PRIVATE -g -O0 -fcheck=all -fbacktrace)
target_link_libraries(test_raw_oneform.x PRIVATE raw_probe)
add_test(NAME test_raw_oneform COMMAND test_raw_oneform.x)
set_tests_properties(test_raw_oneform PROPERTIES TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
