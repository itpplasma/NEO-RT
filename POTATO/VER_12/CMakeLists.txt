cmake_minimum_required (VERSION 2.8.12)

project(POTATO)
enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
add_compile_options(-Og -cpp -g -fbacktrace -Wfatal-errors -std=f2008 -fopenmp)
   # -finit-real=nan -ffpe-trap=invalid,zero,overflow)
add_subdirectory(SRC)

find_package(BLAS)
find_package(LAPACK)

link_libraries(potato ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable (potato.x
    SRC/tt.f90
    )
