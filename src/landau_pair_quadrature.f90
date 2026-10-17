module neort_landau_pair_quadrature
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_gen_laguerre, gauss_legendre
    implicit none
    private

    real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)
    real(dp), parameter :: cm_weight_exponent = 0.5_dp
    real(dp), parameter :: relative_weight_exponent = 0.0_dp
    integer, parameter :: minimum_gaussian_order = 2, minimum_azimuth_order = 4

    ! Nodes and weights are read-only after initialization.
    type, public :: landau_quadrature_t
        real(dp), allocatable :: y(:), wy(:), t(:), wt(:)
        real(dp), allocatable :: mu(:), wmu(:), phi(:), wphi(:)
    end type landau_quadrature_t

    public :: initialize_landau_quadrature, validate_landau_quadrature

contains

    subroutine initialize_landau_quadrature(n_cm, n_rel, n_polar, n_azimuth, &
                                            quadrature)
        integer, intent(in) :: n_cm, n_rel, n_polar, n_azimuth
        type(landau_quadrature_t), intent(out) :: quadrature

        if (min(n_cm, n_rel, n_polar) < minimum_gaussian_order) &
            error stop 'Landau Gaussian quadrature orders must be at least two'
        if (n_azimuth < minimum_azimuth_order) &
            error stop 'Landau azimuth quadrature order must be at least four'
        allocate (quadrature%y(n_cm), quadrature%wy(n_cm))
        allocate (quadrature%t(n_rel), quadrature%wt(n_rel))
        allocate (quadrature%mu(n_polar), quadrature%wmu(n_polar))
        allocate (quadrature%phi(n_azimuth), quadrature%wphi(n_azimuth))
        call gauss_gen_laguerre(n_cm, cm_weight_exponent, &
                               quadrature%y, quadrature%wy)
        call gauss_gen_laguerre(n_rel, relative_weight_exponent, &
                               quadrature%t, quadrature%wt)
        call gauss_legendre(n_polar, quadrature%mu, quadrature%wmu)
        call initialize_azimuth(quadrature%phi, quadrature%wphi)
        call validate_landau_quadrature(quadrature)
    end subroutine initialize_landau_quadrature

    subroutine initialize_azimuth(phi, weights)
        real(dp), intent(out) :: phi(:), weights(:)
        real(dp) :: spacing
        integer :: i

        spacing = two_pi/real(size(phi), dp)
        do i = 1, size(phi)
            phi(i) = real(i - 1, dp)*spacing
        end do
        weights = spacing
    end subroutine initialize_azimuth

    subroutine validate_landau_quadrature(quadrature)
        type(landau_quadrature_t), intent(in) :: quadrature

        call validate_rule(quadrature%y, quadrature%wy, minimum_gaussian_order)
        call validate_rule(quadrature%t, quadrature%wt, minimum_gaussian_order)
        call validate_rule(quadrature%mu, quadrature%wmu, minimum_gaussian_order)
        call validate_rule(quadrature%phi, quadrature%wphi, minimum_azimuth_order)
        if (any(quadrature%y <= 0.0_dp)) &
            error stop 'Landau center-of-mass nodes must be positive'
        if (any(quadrature%t <= 0.0_dp)) &
            error stop 'Landau relative nodes must be positive'
        if (any(abs(quadrature%mu) > 1.0_dp)) &
            error stop 'Landau polar nodes must lie in [-1,1]'
        if (any(quadrature%phi < 0.0_dp) .or. any(quadrature%phi >= two_pi)) &
            error stop 'Landau azimuth nodes must lie in [0,2*pi)'
    end subroutine validate_landau_quadrature

    subroutine validate_rule(nodes, weights, minimum_order)
        real(dp), allocatable, intent(in) :: nodes(:), weights(:)
        integer, intent(in) :: minimum_order

        if (.not. allocated(nodes)) error stop 'Landau quadrature nodes are required'
        if (.not. allocated(weights)) &
            error stop 'Landau quadrature weights are required'
        if (size(nodes) < minimum_order) &
            error stop 'Landau quadrature order is too small'
        if (size(weights) /= size(nodes)) error stop 'Landau quadrature shapes differ'
        if (.not. all(ieee_is_finite(nodes))) &
            error stop 'Landau quadrature nodes must be finite'
        if (.not. all(ieee_is_finite(weights))) &
            error stop 'Landau quadrature weights must be finite'
        if (any(weights <= 0.0_dp)) &
            error stop 'Landau quadrature weights must be positive'
    end subroutine validate_rule

end module neort_landau_pair_quadrature
