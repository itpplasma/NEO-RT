# The field, magfie, spline and kinds/constants sources live in libneo and are
# linked below (LIBNEO::magfie pulls in libneo's core). Only POTATO-specific
# sources are compiled here. binsrc/plag_coeff stay local: POTATO calls them as
# bare external subroutines, whereas libneo wraps them in *_sub modules.
add_library(potato_base
  logging_mod.f90
  potato_input_mod.f90
  alpha_lifetime_mod.f90
  velo.f90
  odeint_allroutines.f
  plag_coeff.f90
  find_all_roots.f90
  sample_matrix_mod.f90
  sample_matrix.f90
  sample_matrix_out.f90
  sub_potato.f90
  equimoments.f90
  sorting.f90
  binsrc.f90
  bmod_pert.f90
  resonant_int.f90
  eqmagprofs.f90
  box_counting.f90
  neo2_reader.f90
)
find_or_fetch(vode)
# vode (dvode_f90_m) uses labeled DO termination, deleted in f2018; build it with
# the permissive legacy standard rather than POTATO's f2018.
target_compile_options(vode PRIVATE -std=legacy)
target_link_libraries(potato_base PUBLIC
vode
LIBNEO::magfie
mc_efit
BLAS::BLAS
LAPACK::LAPACK
HDF5::HDF5
fortplot::fortplot
OpenMP::OpenMP_Fortran
)

add_library(potato
  profile_input.f90
)
target_link_libraries(potato PUBLIC potato_base)

add_library(potato_fixed
  profile_input_fixed.f90
)
target_link_libraries(potato_fixed PUBLIC potato_base)
