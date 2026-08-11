module neort_separatrix
    ! Data-derived regularisation of the logarithmic thin-orbit limit.
    !
    ! For a smooth non-degenerate magnetic maximum the bounce and orbit
    ! integrals have the local form
    !
    !   F(x) = A log(1/x) + B + x [C log(1/x) + D] + O(x**2 log x),
    !
    ! with x = |eta-etatp|/etatp.  The coefficients are fitted from the
    ! smallest numerically resolved orbit distances.  The fit is performed in
    ! the scaled coordinate z=x/xscale, so the least-squares system remains
    ! well-conditioned while the physical distance and its derivative remain
    ! explicit at the call boundary.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_separatrix_frequency_symbolic, only: &
        neort_separatrix_frequency_kernel
    implicit none
    private

    integer, parameter, public :: SEPARATRIX_FIT_POINTS = 8

    type, public :: separatrix_model_t
        real(dp) :: tau(4) = 0.0_dp
        real(dp) :: drift(4) = 0.0_dp
        real(dp) :: xscale = 0.0_dp
        real(dp) :: tau_relative_rms = huge(1.0_dp)
        real(dp) :: drift_relative_rms = huge(1.0_dp)
        logical :: valid = .false.
    end type separatrix_model_t

    public :: fit_separatrix_model
    public :: evaluate_separatrix_model

contains

    subroutine fit_separatrix_model(x, tau_values, drift_values, model, ok)
        real(dp), intent(in) :: x(:), tau_values(:), drift_values(:)
        type(separatrix_model_t), intent(out) :: model
        logical, intent(out) :: ok

        real(dp) :: selected_x(SEPARATRIX_FIT_POINTS)
        real(dp) :: selected_tau(SEPARATRIX_FIT_POINTS)
        real(dp) :: selected_drift(SEPARATRIX_FIT_POINTS)
        real(dp) :: design(SEPARATRIX_FIT_POINTS, 4)
        real(dp) :: tau_coeff(4), drift_coeff(4)
        real(dp) :: tau_leading(2)
        real(dp) :: tau_rms, drift_rms
        integer :: selected, k

        model = separatrix_model_t()
        ok = .false.
        if (size(x) /= size(tau_values) .or. size(x) /= size(drift_values)) return
        if (size(x) < 4) return

        call select_nearest_boundary(x, tau_values, drift_values, selected_x, &
            selected_tau, selected_drift, selected)
        if (selected < 4) return

        model%xscale = maxval(selected_x(1:selected))
        if (.not. ieee_is_finite(model%xscale) .or. model%xscale <= 0.0_dp) return

        do k = 1, selected
            if (selected_x(k) <= 0.0_dp) return
            associate (z => selected_x(k)/model%xscale)
                design(k, 1) = log(1.0_dp/z)
                design(k, 2) = 1.0_dp
                design(k, 3) = z*design(k, 1)
                design(k, 4) = z
            end associate
        end do

        call least_squares(design, selected_tau, selected, tau_coeff, ok, tau_rms)
        if (.not. ok) return
        call least_squares(design, selected_drift, selected, drift_coeff, ok, drift_rms)
        if (.not. ok) return

        if (.not. ieee_is_finite(tau_coeff(1)) .or. tau_coeff(1) <= 0.0_dp) then
            ! The correction columns are below the numerical resolution of
            ! some orbit integrations.  In that case their unconstrained fit
            ! can steal the positive logarithmic coefficient.  Retain the
            ! same asymptotic boundary window, but use the resolved leading
            ! logarithmic submodel instead of accepting an unphysical fit.
            call least_squares(design(:, 1:2), selected_tau, selected, &
                tau_leading, ok, tau_rms)
            if (.not. ok .or. tau_leading(1) <= 0.0_dp) return
            tau_coeff = [tau_leading(1), tau_leading(2), 0.0_dp, 0.0_dp]
        end if
        if (.not. all(ieee_is_finite(tau_coeff))) then
            ok = .false.
            return
        end if
        if (.not. all(ieee_is_finite(drift_coeff))) then
            ok = .false.
            return
        end if

        model%tau = tau_coeff
        model%drift = drift_coeff
        model%tau_relative_rms = tau_rms
        model%drift_relative_rms = drift_rms
        model%valid = .true.
    end subroutine fit_separatrix_model

    subroutine evaluate_separatrix_model(model, x, tau, dtau_dx, omega, &
            domega_dx, drift, ddrift_dx, ok)
        type(separatrix_model_t), intent(in) :: model
        real(dp), intent(in) :: x
        real(dp), intent(out) :: tau, dtau_dx, omega, domega_dx, drift, ddrift_dx
        logical, intent(out) :: ok
        real(dp) :: z, dtau_dz, ddrift_dz, domega_dz

        tau = 0.0_dp
        dtau_dx = 0.0_dp
        omega = 0.0_dp
        domega_dx = 0.0_dp
        drift = 0.0_dp
        ddrift_dx = 0.0_dp
        ok = .false.
        if (.not. model%valid) return
        if (x <= 0.0_dp .or. x > model%xscale) return
        z = x/model%xscale
        call neort_separatrix_frequency_kernel(z, model%tau(1), model%tau(2), &
            model%tau(3), model%tau(4), model%drift(1), model%drift(2), &
            model%drift(3), model%drift(4), tau, dtau_dz, omega, domega_dz, &
            drift, ddrift_dz)
        dtau_dx = dtau_dz/model%xscale
        domega_dx = domega_dz/model%xscale
        ddrift_dx = ddrift_dz/model%xscale
        ok = ieee_is_finite(tau) .and. ieee_is_finite(dtau_dx) .and. &
            ieee_is_finite(omega) .and. ieee_is_finite(domega_dx) .and. &
            ieee_is_finite(drift) .and. ieee_is_finite(ddrift_dx)
    end subroutine evaluate_separatrix_model

    subroutine select_nearest_boundary(x, tau_values, drift_values, selected_x, &
            selected_tau, selected_drift, selected)
        real(dp), intent(in) :: x(:), tau_values(:), drift_values(:)
        real(dp), intent(out) :: selected_x(SEPARATRIX_FIT_POINTS)
        real(dp), intent(out) :: selected_tau(SEPARATRIX_FIT_POINTS)
        real(dp), intent(out) :: selected_drift(SEPARATRIX_FIT_POINTS)
        integer, intent(out) :: selected

        integer :: order(size(x)), i, j, candidate
        integer :: swap
        logical :: usable(size(x))

        selected_x = 0.0_dp
        selected_tau = 0.0_dp
        selected_drift = 0.0_dp
        selected = min(SEPARATRIX_FIT_POINTS, size(x))
        order = [(i, i=1, size(x))]
        usable = ieee_is_finite(x) .and. ieee_is_finite(tau_values) .and. &
            ieee_is_finite(drift_values)
        do i = 1, size(x) - 1
            do j = i + 1, size(x)
                if (usable(order(j))) then
                    if (.not. usable(order(i))) then
                        swap = order(i)
                        order(i) = order(j)
                        order(j) = swap
                    else if (x(order(j)) < x(order(i))) then
                        swap = order(i)
                        order(i) = order(j)
                        order(j) = swap
                    end if
                end if
            end do
        end do
        candidate = 0
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(order(i))) .or. &
                    .not. ieee_is_finite(tau_values(order(i))) .or. &
                    .not. ieee_is_finite(drift_values(order(i)))) cycle
            if (x(order(i)) <= 0.0_dp) cycle
            candidate = candidate + 1
            selected_x(candidate) = x(order(i))
            selected_tau(candidate) = tau_values(order(i))
            selected_drift(candidate) = drift_values(order(i))
            if (candidate == selected) exit
        end do
        selected = candidate
    end subroutine select_nearest_boundary

    subroutine least_squares(design, values, n, coefficients, ok, relative_rms)
        real(dp), intent(in) :: design(:, :), values(:)
        integer, intent(in) :: n
        real(dp), intent(out) :: coefficients(:), relative_rms
        logical, intent(out) :: ok
        integer :: ncolumns
        real(dp) :: q(n, size(coefficients)), r(size(coefficients), size(coefficients))
        real(dp) :: work(n, size(coefficients)), vector(n), rhs(size(coefficients)), residual(n)
        real(dp) :: column_scale(size(coefficients)), scale, correction
        integer :: i, j, k

        coefficients = 0.0_dp
        relative_rms = huge(1.0_dp)
        ok = .false.
        q = 0.0_dp
        r = 0.0_dp
        ncolumns = size(coefficients)
        if (size(design, 2) /= ncolumns .or. n < ncolumns) return
        work = design(1:n, 1:ncolumns)
        do j = 1, ncolumns
            column_scale(j) = sqrt(dot_product(work(:, j), work(:, j)))
            if (.not. ieee_is_finite(column_scale(j)) .or. &
                    column_scale(j) <= tiny(1.0_dp)) return
            work(:, j) = work(:, j)/column_scale(j)
        end do
        do j = 1, ncolumns
            vector = work(:, j)
            do k = 1, j - 1
                r(k, j) = dot_product(q(:, k), vector)
                vector = vector - r(k, j)*q(:, k)
            end do
            ! Reorthogonalisation prevents the close logarithmic columns from
            ! losing rank when the overlap contains only a narrow x interval.
            do k = 1, j - 1
                correction = dot_product(q(:, k), vector)
                r(k, j) = r(k, j) + correction
                vector = vector - correction*q(:, k)
            end do
            scale = sqrt(dot_product(vector, vector))
            if (.not. ieee_is_finite(scale) .or. &
                    scale <= sqrt(epsilon(1.0_dp))) return
            r(j, j) = scale
            q(:, j) = vector/scale
        end do

        rhs = 0.0_dp
        do j = 1, ncolumns
            rhs(j) = dot_product(q(:, j), values(1:n))
        end do
        do i = ncolumns, 1, -1
            coefficients(i) = rhs(i)
            do j = i + 1, ncolumns
                coefficients(i) = coefficients(i) - r(i, j)*coefficients(j)
            end do
            coefficients(i) = coefficients(i)/r(i, i)
        end do
        coefficients = coefficients/column_scale

        residual = matmul(design(1:n, :), coefficients) - values(1:n)
        scale = max(sqrt(dot_product(values(1:n), values(1:n))), tiny(1.0_dp))
        relative_rms = sqrt(dot_product(residual, residual))/scale
        ok = ieee_is_finite(relative_rms)
    end subroutine least_squares

end module neort_separatrix
