program test_harmonic_bounds
    use iso_fortran_env, only: dp => real64
    use neort, only: harmonic_bounds, configured_mth_max_abs => mth_max_abs
    use neort_config, only: config_t, set_config

    implicit none

    integer :: mth_min, mth_max
    type(config_t) :: config

    call harmonic_bounds(3, 1.2_dp, -1, mth_min, mth_max)
    if (mth_min /= -8 .or. mth_max /= 8) then
        error stop "automatic harmonic bounds changed"
    end if

    call harmonic_bounds(3, 1.2_dp, 0, mth_min, mth_max)
    if (mth_min /= 0 .or. mth_max /= 0) then
        error stop "mth_max_abs=0 must select only mth=0"
    end if

    call harmonic_bounds(3, 1.2_dp, 5, mth_min, mth_max)
    if (mth_min /= -5 .or. mth_max /= 5) then
        error stop "mth_max_abs=5 must select mth=-5:5"
    end if

    config%mth_max_abs = 5
    call set_config(config)
    if (configured_mth_max_abs /= 5) error stop "config_t did not set mth_max_abs"

    config = config_t()
    call set_config(config)
    if (configured_mth_max_abs /= -1) error stop "config_t default changed"
end program test_harmonic_bounds
