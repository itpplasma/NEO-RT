gfortran -J OBJS --check=all -g -p -pg -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow SRC/alpha_lifetime_mod.f90 SRC/chamb_divB0.f90 SRC/magfie_cyl.f90 SRC/period_mod.f90 SRC/input_files.f90 SRC/field_c_mod.f90 SRC/field_mod.f90 SRC/inthecore_mod.f90 SRC/field_eq_mod.f90 SRC/field_divB0.f90 SRC/bdivfree_mod.f90 SRC/amn_mod.f90 SRC/theta_rz_mod.f90 SRC/extract_fluxcoord_mod.f90 SRC/bdivfree.f90 SRC/spline5_RZ.f90 SRC/velo.f90 SRC/odeint_allroutines.f SRC/plag_coeff.f90 SRC/find_all_roots.f90 SRC/sample_matrix.f90 SRC/sample_matrix_out.f90 SRC/sub_potato.f90 SRC/equimoments.f90 SRC/sorting.f90 SRC/binsrc.f90 SRC/profile_input.f90 SRC/bmod_pert.f90 SRC/resonant_int.f90 SRC/eqmagprofs.f90 SRC/tt.f90 -llapack
