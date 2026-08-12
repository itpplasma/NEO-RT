module neort_separatrix
    ! Data-derived regularisation of the logarithmic thin-orbit limit.
    !
    ! Let x = |eta-etatp|/etatp be the pitch distance from the trapped-passing
    ! boundary.  Around a smooth non-degenerate magnetic maximum the orbit
    ! spends time (1/lambda) log(1/x) near the barrier, with lambda the
    ! linear instability exponent of that saddle, and the rest of the period
    ! is analytic in x.  Every orbit integral therefore shares one endpoint
    ! form,
    !
    !   I(x) = A log(1/x) + B + x [C log(1/x) + D] + O(x**2 log(1/x)),
    !
    ! and this module fits that form to two integrals at once: the bounce
    ! time tau, and the same integral weighted by the magnetic-drift
    ! integrand.  The bounce-averaged drift is the ratio of the two, which is
    ! how it is defined (Mackenbach et al., Phys. Plasmas 30, 093901 (2023),
    ! Eqs. 7/10/13).  Both leading coefficients come from the same barrier
    ! passage and differ only by the value of the drift integrand there, so
    ! the ratio has a finite separatrix limit num_a/tau_a.
    !
    ! The ratio is kept intact instead of being re-expanded.  Expanding it
    ! gives a series in 1/log(1/x) whose successive terms fall only by
    ! tau_b/(tau_a log(1/x)) - about 1.1 for the passing orbits of the shipped
    ! base example at its innermost resolved sample, so there the series does
    ! not even converge.  Where it does converge it does so far too slowly to
    ! be truncated after two terms and then extrapolated across the ten or
    ! more decades between the resolved samples and the boundary.  Fitting
    ! numerator and denominator separately resums it exactly.
    !
    ! The fit runs in the scaled coordinate z=x/xscale so the least-squares
    ! system stays well-conditioned, while the physical distance and its
    ! derivative remain explicit at the call boundary.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_separatrix_frequency_symbolic, only: &
        neort_separatrix_frequency_kernel
    implicit none
    private

    ! Minimum natural-log span log(x_max/x_min) of the fitted boundary window.
    ! The leading coefficient A is resolved only through the contrast of the
    ! log column across the window, so with relative sampling noise sigma and
    ! I/A of order ten the relative error of A is roughly 10*sigma/span.  The
    ! competing model error is the neglected O(x**2 log(1/x)) term, which is
    ! O(x_max**2) relative and stays below 1e-9 out to x_max = 1e-4, so the
    ! window can be widened freely.  A span of two e-folds keeps an order of
    ! margin over the measured sampling noise without approaching that limit.
    ! A fixed sample count cannot do this: it inherits whatever spacing the
    ! caller's pitch grid happens to have.
    real(dp), parameter, public :: SEPARATRIX_FIT_LOG_SPAN = 2.0_dp
    integer, parameter :: SEPARATRIX_FIT_TERMS = 4

    type, public :: separatrix_model_t
        real(dp) :: tau(4) = 0.0_dp
        real(dp) :: numerator(4) = 0.0_dp
        real(dp) :: xscale = 0.0_dp
        real(dp) :: tau_relative_rms = huge(1.0_dp)
        real(dp) :: numerator_relative_rms = huge(1.0_dp)
        logical :: valid = .false.
    end type separatrix_model_t

    public :: fit_separatrix_model
    public :: evaluate_separatrix_model

contains

    subroutine fit_separatrix_model(x, tau_values, drift_values, model, ok)
        ! drift_values are bounce averages, i.e. already divided by the bounce
        ! time.  The endpoint form holds for the un-normalised integral, so
        ! the numerator is reconstructed here and fitted on the same basis.
        real(dp), intent(in) :: x(:), tau_values(:), drift_values(:)
        type(separatrix_model_t), intent(out) :: model
        logical, intent(out) :: ok

        real(dp), allocatable :: selected_x(:), selected_tau(:), selected_drift(:)
        real(dp), allocatable :: selected_numerator(:)
        real(dp) :: design(size(x), SEPARATRIX_FIT_TERMS)
        real(dp) :: tau_coeff(SEPARATRIX_FIT_TERMS), num_coeff(SEPARATRIX_FIT_TERMS)
        real(dp) :: tau_leading(2), num_leading(2)
        real(dp) :: tau_rms, num_rms, xscale
        integer :: selected, k, nfit, first_window
        logical :: full_ok, tau_ok, num_ok

        model = separatrix_model_t()
        ok = .false.
        if (size(x) /= size(tau_values) .or. size(x) /= size(drift_values)) return
        if (size(x) < SEPARATRIX_FIT_TERMS) return

        allocate(selected_x(size(x)), selected_tau(size(x)), selected_drift(size(x)))
        allocate(selected_numerator(size(x)))
        call select_nearest_boundary(x, tau_values, drift_values, selected_x, &
            selected_tau, selected_drift, selected)
        if (selected < SEPARATRIX_FIT_TERMS) return
        selected_numerator = selected_drift*selected_tau

        first_window = resolved_window(selected_x, selected)

        ! The smallest samples can be under-resolved by the orbit quadrature
        ! and can reverse the physical positive logarithmic trend.  Increase
        ! the boundary window until the same asymptotic model is resolved;
        ! no fixed outer cutoff or hand-tuned coefficient is introduced.
        do nfit = first_window, selected
            xscale = maxval(selected_x(1:nfit))
            if (.not. ieee_is_finite(xscale) .or. xscale <= 0.0_dp) cycle

            do k = 1, nfit
                associate (z => selected_x(k)/xscale)
                    design(k, 1) = log(1.0_dp/z)
                    design(k, 2) = 1.0_dp
                    design(k, 3) = z*design(k, 1)
                    design(k, 4) = z
                end associate
            end do

            call least_squares(design, selected_tau, nfit, tau_coeff, tau_ok, tau_rms)
            call least_squares(design, selected_numerator, nfit, num_coeff, &
                num_ok, num_rms)
            full_ok = tau_ok .and. num_ok .and. all(ieee_is_finite(tau_coeff)) .and. &
                all(ieee_is_finite(num_coeff)) .and. tau_coeff(1) > 0.0_dp
            if (.not. full_ok) then
                ! The correction columns are below the numerical resolution
                ! of this window, or their unconstrained fit stole the
                ! positive logarithmic coefficient.  Retain the same
                ! boundary window and use the resolved leading model, which
                ! still gives the ratio its finite separatrix limit.
                call least_squares(design(:, 1:2), selected_tau, nfit, &
                    tau_leading, tau_ok, tau_rms)
                call least_squares(design(:, 1:2), selected_numerator, nfit, &
                    num_leading, num_ok, num_rms)
                if (.not. tau_ok .or. .not. num_ok .or. &
                        .not. all(ieee_is_finite(tau_leading)) .or. &
                        .not. all(ieee_is_finite(num_leading)) .or. &
                        tau_leading(1) <= 0.0_dp) cycle
                tau_coeff = [tau_leading(1), tau_leading(2), 0.0_dp, 0.0_dp]
                num_coeff = [num_leading(1), num_leading(2), 0.0_dp, 0.0_dp]
            end if

            model%xscale = xscale
            model%tau = tau_coeff
            model%numerator = num_coeff
            model%tau_relative_rms = tau_rms
            model%numerator_relative_rms = num_rms
            model%valid = .true.
            ok = .true.
            return
        end do
        ! No window produced a usable model.  Leave xscale at zero so callers
        ! that gate on x <= xscale fall back to their own representation
        ! instead of entering an invalid model.
        model = separatrix_model_t()
    end subroutine fit_separatrix_model

    pure function resolved_window(selected_x, selected) result(nfit)
        ! Smallest sample count that both spans SEPARATRIX_FIT_LOG_SPAN and
        ! leaves the fit over-determined, so that its residual is a residual
        ! and not an artefact of an exactly determined system.
        real(dp), intent(in) :: selected_x(:)
        integer, intent(in) :: selected
        integer :: nfit, k

        nfit = selected
        do k = min(SEPARATRIX_FIT_TERMS + 1, selected), selected
            if (log(selected_x(k)/selected_x(1)) >= SEPARATRIX_FIT_LOG_SPAN) then
                nfit = k
                return
            end if
        end do
    end function resolved_window

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
            model%tau(3), model%tau(4), model%numerator(1), model%numerator(2), &
            model%numerator(3), model%numerator(4), tau, dtau_dz, omega, domega_dz, &
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
        real(dp), intent(out) :: selected_x(:)
        real(dp), intent(out) :: selected_tau(:)
        real(dp), intent(out) :: selected_drift(:)
        integer, intent(out) :: selected

        integer :: order(size(x)), i, j, candidate
        integer :: swap
        logical :: usable(size(x))

        selected_x = 0.0_dp
        selected_tau = 0.0_dp
        selected_drift = 0.0_dp
        selected = 0
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
            if (candidate == size(selected_x)) exit
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
        real(dp) :: rank_tolerance
        integer :: i, j, k

        coefficients = 0.0_dp
        relative_rms = huge(1.0_dp)
        ok = .false.
        q = 0.0_dp
        r = 0.0_dp
        ncolumns = size(coefficients)
        if (size(design, 2) /= ncolumns .or. n < ncolumns) return
        ! The four-term correction basis is deliberately rejected when its
        ! extra columns are below square-root machine precision.  The
        ! two-term leading model has only one contrast direction and can
        ! still resolve a narrower physical window down to roundoff.
        rank_tolerance = epsilon(1.0_dp)
        if (ncolumns > 2) rank_tolerance = sqrt(epsilon(1.0_dp))
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
                    scale <= rank_tolerance) return
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
