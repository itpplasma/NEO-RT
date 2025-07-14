add_library(potato_base
  libneo_kinds.f90
  math_constants.f90
  alpha_lifetime_mod.f90
  chamb_divB0.f90
  magfie_cyl.f90
  period_mod.f90
  input_files.f90
  field_c_mod.f90
  field_mod.f90
  inthecore_mod.f90
  field_eq_mod.f90
  field_divB0.f90
  bdivfree_mod.f90
  amn_mod.f90
  theta_rz_mod.f90
  extract_fluxcoord_mod.f90
  bdivfree.f90
  spl_three_to_five.f90
  spline5_RZ.f90
  velo.f90
  odeint_allroutines.f
  plag_coeff.f90
  find_all_roots.f90
  sample_matrix.f90
  sample_matrix_out.f90
  sub_potato.f90
  equimoments.f90
  sorting.f90
  binsrc.f90
  bmod_pert.f90
  resonant_int.f90
  eqmagprofs.f90
  # box_counting.f90  # Disabled due to VODE conflict with dvode_f90_m
)
# Use VODE from parent project if available, otherwise fetch
if(NOT TARGET vode)
  find_or_fetch(vode)
endif()
target_link_libraries(potato_base PUBLIC vode)

add_library(potato
  profile_input.f90
)
target_link_libraries(potato PUBLIC potato_base)

add_library(potato_fixed
  profile_input_fixed.f90
)
target_link_libraries(potato_fixed PUBLIC potato_base)

# Remove circular dependency - we'll handle this differently
