add_library(potato_base
  ${LIBNEO_SOURCE_DIR}/src/libneo_kinds.f90
  logging_mod.f90
  potato_input_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/math_constants.f90
  alpha_lifetime_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/MC/chamb_divB0.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/magfie_cyl.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/period_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/input_files.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/field_c_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/field_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/inthecore_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/field_eq_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/field_divB0.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/bdivfree_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/amn_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/theta_rz_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/extract_fluxcoord_mod.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/bdivfree.f90
  ${LIBNEO_SOURCE_DIR}/src/spl_three_to_five.f90
  ${LIBNEO_SOURCE_DIR}/src/magfie/spline5_RZ.f90
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
