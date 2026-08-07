program test_gc_full_fow_scales_boundary_cut
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_buchholz_boundary_symbolic, only: &
        evaluate_neort_buchholz_boundary_limits
    use neort_buchholz_cut_symbolic, only: evaluate_neort_buchholz_cut
    use neort_buchholz_cdot_symbolic, only: evaluate_neort_buchholz_cdot
    use neort_cylindrical_crossing_symbolic, only: &
        evaluate_neort_cylindrical_crossing_density
    use neort_full_fow_refinement_symbolic, only: &
        evaluate_neort_full_fow_refinement_scales
    use util_for_test, only: pass_test
    implicit none

    real(dp) :: scales(14), boundary(13), cut(5), cdot(3), crossing(5)
    real(dp), parameter :: hm_real = 0.6_dp, hm_imag = -0.8_dp
    real(dp), parameter :: distance = 0.04_dp
    real(dp), parameter :: distance_scale = 0.25_dp

    call evaluate_neort_full_fow_refinement_scales(2.0_dp, -3.0_dp, 0.5_dp, &
        -0.8_dp, 0.7_dp, 4.0_dp, 1.2_dp, -0.4_dp, hm_real, hm_imag, 0.5_dp, &
        -2.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp, &
        1.0_dp, 2.0_dp, 2.0_dp, 0.5_dp, 2.0_dp, 3.0_dp, 2.0_dp, 3.0_dp, &
        scales(1), scales(2), scales(3), scales(4), scales(5), scales(6), &
        scales(7), scales(8), scales(9), scales(10), scales(11), scales(12), &
        scales(13), scales(14))
    call require_close('scale R', scales(1), 2.0_dp)
    call require_close('scale Z', scales(2), 1.5_dp)
    call require_close('scale phi', scales(3), 0.5_dp)
    call require_close('scale p_parallel', scales(4), 0.4_dp)
    call require_close('scale mu', scales(5), 0.7_dp)
    call require_close('scale tau', scales(6), 2.0_dp)
    call require_close('scale omega_b', scales(7), 0.6_dp)
    call require_close('scale omega_phi', scales(8), 0.4_dp)
    call require_close('scale |Hm| real', scales(9), 0.6_dp)
    call require_close('scale |Hm| imag', scales(10), 0.4_dp)
    call require_close('scale |Hm|^2', scales(11), 0.5_dp)
    call require_close('scale root derivative', scales(12), 1.0_dp)
    call require_close('scale torque density', scales(13), 1.0_dp)
    call require_close('scale torque', scales(14), 1.0_dp)

    call evaluate_neort_buchholz_boundary_limits(distance, distance_scale, &
        3.0_dp, 2.0_dp, 2.0_dp, boundary(1), boundary(2), boundary(3), &
        boundary(4), boundary(5), boundary(6), boundary(7), boundary(8), &
        boundary(9), boundary(10), boundary(11), boundary(12), boundary(13))
    call require_close('C_tau', boundary(1), 0.5_dp)
    call require_close('regular finite form', boundary(2), 3.0_dp)
    call require_close('reflecting sqrt form', boundary(3), 0.8_dp)
    call require_close('reflecting derivative form', boundary(4), 10.0_dp)
    call require_close('separatrix log form', boundary(5), &
        -0.5_dp*log(distance/distance_scale))
    call require_close('separatrix positive tau form', boundary(6), &
        0.5_dp*log(distance_scale/distance))
    call require_close('X-point log form', boundary(7), &
        -log(distance/distance_scale))
    call require_close('X-point positive tau form', boundary(8), &
        log(distance_scale/distance))
    call require_close('regular class residual', boundary(9), 0.0_dp)
    call require_close('reflecting class residual', boundary(10), 0.0_dp)
    call require_close('separatrix class residual', boundary(11), 0.0_dp)
    call require_close('paired coefficient residual', boundary(12), 0.0_dp)
    call require_close('paired tau residual', boundary(13), 0.0_dp)

    call evaluate_neort_buchholz_cut(1.0_dp, 2.0_dp, 3.0_dp, 0.2_dp, 0.3_dp, &
        0.4_dp, 0.0_dp, 1.0_dp, 2.0_dp, cut(1), cut(2), cut(3), cut(4), &
        cut(5))
    call require_close('cut cross R', cut(1), -0.1_dp)
    call require_close('cut cross phi', cut(2), 0.2_dp)
    call require_close('cut cross Z', cut(3), -0.1_dp)
    call require_close('Buchholz cut C on section', cut(4), 0.0_dp)
    call require_close('cut identity residual', cut(5), 0.0_dp)

    call evaluate_neort_buchholz_cdot(1.0_dp, 2.0_dp, 3.0_dp, 0.3_dp, 0.4_dp, &
        0.5_dp, cdot(1), cdot(2), cdot(3))
    call require_close('Cdot at C=0', cdot(1), 2.6_dp)
    call require_close('Cdot transversality', cdot(2), 2.6_dp)
    call require_close('Cdot orientation', cdot(3), 2.6_dp)

    call evaluate_neort_cylindrical_crossing_density(-4.0_dp, cdot(1), &
        crossing(1), crossing(2), crossing(3), crossing(4), crossing(5))
    call require_close('signed R Bparallel_star', crossing(1), -4.0_dp)
    call require_close('signed crossing density', crossing(2), -10.4_dp)
    call require_close('positive phase-space crossing density', crossing(5), 10.4_dp)

    write (*, '(A)') 'test_gc_full_fow_scales_boundary_cut OK'
    call pass_test

contains

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual - expected) > 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'full-FOW scale/boundary/cut oracle failed'
        end if
    end subroutine require_close

end program test_gc_full_fow_scales_boundary_cut
