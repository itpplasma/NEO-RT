cmake_minimum_required (VERSION 3.22)

project(POTATO)
enable_language(Fortran)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
add_compile_options(-cpp -g -fbacktrace -Wfatal-errors -std=f2008)
add_compile_options(-O3 -march=native -mtune=native)

include(Util)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS Fortran)

include(FetchContent)
FetchContent_Declare(
    fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG main
)
FetchContent_MakeAvailable(fortplot)

add_subdirectory(SRC)

add_executable (potato.x
    SRC/tt.f90
)
target_link_libraries(potato.x
potato
)


add_executable (test_profile.x
   SRC/test_profile.f90
)
target_link_libraries(test_profile.x
potato
)

add_executable (test_profile_fixed.x
   SRC/test_profile.f90
)
target_link_libraries(test_profile_fixed.x
potato_fixed
)

add_executable(check_neo2_profiles.x
   SRC/check_neo2_profiles.f90
)
target_link_libraries(check_neo2_profiles.x
   potato
)
