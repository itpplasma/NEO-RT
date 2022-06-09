cmake_minimum_required (VERSION 2.8.12)

project(POTATO)
enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
add_compile_options(-cpp -g -fbacktrace -Wfatal-errors -std=f2008)
add_compile_options(-O3 -march=native -mtune=native) #GNU
# add_compile_options(-O3 -march=native -mtune=native -g -fbacktrace -finit-real=snan -fsignaling-nans -ffpe-trap=invalid,zero,overflow) #GNU
   # -finit-real=nan -ffpe-trap=invalid,zero,overflow)

add_subdirectory(SRC)
add_subdirectory("../../../contrib/vode" "vode")
add_dependencies(potato vode)

find_package(BLAS)
find_package(LAPACK)


add_executable (potato.x
    SRC/tt.f90
)
target_link_libraries(potato.x
potato
${BLAS_LIBRARIES}
${LAPACK_LIBRARIES}
vode)


add_executable (test_profile.x
   SRC/test_profile.f90
)
target_link_libraries(test_profile.x
potato
${BLAS_LIBRARIES}
${LAPACK_LIBRARIES}
vode)

add_executable (test_profile_fixed.x
   SRC/test_profile.f90
)
target_link_libraries(test_profile_fixed.x
potato_fixed
${BLAS_LIBRARIES}
${LAPACK_LIBRARIES}
vode)
