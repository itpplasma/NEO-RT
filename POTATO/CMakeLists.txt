cmake_minimum_required (VERSION 3.10)

project(POTATO)
enable_language(Fortran)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

# Override parent project's strict error flags for POTATO (legacy code)
set(CMAKE_Fortran_FLAGS "")
add_compile_options(-cpp -g -fbacktrace -std=f2008)
add_compile_options(-O3 -march=native -mtune=native)
# Allow warnings without stopping compilation
add_compile_options(-w)

include(Util)
find_package(BLAS)
find_package(LAPACK)

add_subdirectory(SRC)

# Temporarily disable potato.x due to missing dependencies
# add_executable (potato.x
#     SRC/tt.f90
# )
# target_link_libraries(potato.x
# potato
# ${BLAS_LIBRARIES}
# ${LAPACK_LIBRARIES}
# )


# Temporarily disable test_profile.x due to missing dependencies
# add_executable (test_profile.x
#    SRC/test_profile.f90
# )
# target_link_libraries(test_profile.x
# potato
# ${BLAS_LIBRARIES}
# ${LAPACK_LIBRARIES}
# )

# Temporarily disable test_profile_fixed.x due to missing dependencies
# add_executable (test_profile_fixed.x
#    SRC/test_profile.f90
# )
# target_link_libraries(test_profile_fixed.x
# potato_fixed
# ${BLAS_LIBRARIES}
# ${LAPACK_LIBRARIES}
# )
