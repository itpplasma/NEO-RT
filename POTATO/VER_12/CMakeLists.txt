cmake_minimum_required (VERSION 2.8.12)

project(POTATO)
enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
add_compile_options(-Og -cpp -g -fbacktrace -Wfatal-errors -std=f2008 -fopenmp)
set(CMAKE_Fortran_FLAGS " -O3 -march=native -mtune=native -g -fbacktrace -finit-real=snan -fsignaling-nans -ffpe-trap=invalid,zero,overflow") #GNU
   # -finit-real=nan -ffpe-trap=invalid,zero,overflow)

add_subdirectory(SRC)
add_subdirectory("../../../contrib/vode" "vode")
add_dependencies(potato vode)

find_package(BLAS)
find_package(LAPACK)

link_libraries(potato
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   vode)

add_executable (potato.x
    SRC/tt.f90
    )
