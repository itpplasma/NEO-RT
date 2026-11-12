# spline does not export its Fortran module folder.
get_target_property(_potato_spline_modules spline Fortran_MODULE_DIRECTORY)
if(NOT _potato_spline_modules)
  get_target_property(_potato_spline_modules spline BINARY_DIR)
endif()
# Compile the same reader sources as NEO-RT, with OpenMP enabled for its caches.
add_library(potato_boozer_reader STATIC
  ../../src/logging.f90
  ../../src/util.f90
  ../../src/do_magfie_standalone.f90
)
target_link_libraries(potato_boozer_reader PUBLIC
  spline LIBNEO::boozer LIBNEO::coordinates OpenMP::OpenMP_Fortran
)
target_include_directories(potato_boozer_reader PUBLIC
  "${PROJECT_BINARY_DIR}" "${_potato_spline_modules}"
)

add_library(potato_base
  ${CMAKE_CURRENT_SOURCE_DIR}/../../src/kinetic_energy_selection.f90
  ${libneo_SOURCE_DIR}/src/util/libneo_kinds.f90
  logging_mod.f90
  wall_loss_mod.f90
  potato_input_mod.f90
  boozer_perturbation_mod.f90
  potato_invariants_mod.f90
  ${libneo_SOURCE_DIR}/src/util/math_constants.f90
  alpha_lifetime_mod.f90
  ${libneo_SOURCE_DIR}/src/MC/chamb_divB0.f90
  ${libneo_SOURCE_DIR}/src/magfie/magfie_cyl.f90
  ${libneo_SOURCE_DIR}/src/magfie/period_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/input_files.f90
  ${libneo_SOURCE_DIR}/src/magfie/field_c_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/field_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/inthecore_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/field_eq_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/field_divB0.f90
  ${libneo_SOURCE_DIR}/src/magfie/bdivfree_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/amn_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/theta_rz_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/extract_fluxcoord_mod.f90
  ${libneo_SOURCE_DIR}/src/magfie/bdivfree.f90
  ${libneo_SOURCE_DIR}/src/spl_three_to_five.f90
  ${libneo_SOURCE_DIR}/src/magfie/spline5_RZ.f90
  velo.f90
  odeint_allroutines.f
  plag_coeff.f90
  find_all_roots.f90
  sample_matrix_mod.f90
  sample_matrix.f90
  sample_matrix_out.f90
  poicut_tracing.f90
  sub_potato.f90
  profile_data_mod.f90
  equimoments.f90
  sorting.f90
  binsrc.f90
  bmod_pert.f90
  resonance_mode_bounds_mod.f90
  resonant_int.f90
  eqmagprofs.f90
  box_counting.f90
  neo2_reader.f90
)
# Keep the nonfinite coordinate guard effective under Fast optimization.
set_source_files_properties(bmod_pert.f90 PROPERTIES
  COMPILE_OPTIONS "-fno-finite-math-only")

itp_add_dependency(vode REF 5797785f0829334e758bcd5fe650b26abbe3db90
  TARGET vode)
# vode (dvode_f90_m) uses labeled DO termination, deleted in f2018; build it with
# the permissive legacy standard rather than POTATO's f2018.
target_compile_options(vode PRIVATE -std=legacy)
target_link_libraries(potato_base PUBLIC
potato_boozer_reader
spline
vode
BLAS::BLAS
LAPACK::LAPACK
LIBNEO::hdf5_tools
fortplot::fortplot
OpenMP::OpenMP_Fortran
)
target_include_directories(potato_base PUBLIC "${_potato_spline_modules}")

add_library(potato
  profile_input.f90
)
target_link_libraries(potato PUBLIC potato_base)

add_library(potato_fixed
  profile_input_fixed.f90
)
target_link_libraries(potato_fixed PUBLIC potato_base)
