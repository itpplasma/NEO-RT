module neort_conserving_native
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre, gauss_gen_laguerre
    use neort_conserving_basis, only: conserving_basis_t, &
                                  initialize_conserving_basis, evaluate_conserving_basis
    use neort_landau_pair_quadrature, only: landau_quadrature_t
    use neort_landau_projection, only: assemble_landau_projection
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: conserving_species_t
        real(dp) :: mass_kg = 0, charge_e = 0, density_m3 = 0
        real(dp) :: temperature_ev = 0, vth_m_s = 0, beta = 0, A1 = 0, A2 = 0
        real(dp) :: coulomb_log = 0
        type(conserving_basis_t) :: basis
        real(dp), allocatable :: mass(:, :), streaming(:, :), precession(:, :)
        complex(dp), allocatable :: hamiltonian(:, :), hamiltonian2(:, :)
        real(dp), allocatable :: theta(:), theta_weight(:), jacobian_si(:)
        real(dp), allocatable :: bmod_native(:)
    end type conserving_species_t
    public :: prepare_conserving_species, assemble_conserving_collisions
    public :: conserving_native_point

contains

    ! Caller has installed this selected species' complete native context.
    ! No process-global species state is installed from a parallel callback.
    subroutine prepare_conserving_species(speed_steps, pitch_steps, umax, &
                                          theta_order, velocity_order, cache, info)
        use neort_lorentz_native, only: validate_lorentz_geometry
        use neort_species, only: mixture_enabled, mixture_at_s, selected_species, &
                                 species_mixture
        use neort_profiles, only: vth, A1, A2
        use driftorbit, only: etatp, etadt, mph
        use do_magfie_mod, only: psi_pr, q, sign_theta
        use neort_collision_rates, only: collision_elementary_charge
        use util, only: qi, c
        integer, intent(in) :: speed_steps, pitch_steps, theta_order, velocity_order
        real(dp), intent(in) :: umax
        type(conserving_species_t), intent(out) :: cache
        integer, intent(out) :: info
        real(dp), allocatable :: nodes(:), weights(:), tail(:), wtail(:)
        real(dp), allocatable :: values(:), gradient(:, :), dvalues(:)
        real(dp), allocatable :: edges(:), mu(:), wm(:)
        real(dp) :: bmod, jac, ht, db, u, wu, xi, wxi, span, lo, hi
        real(dp) :: velocity(3), omega, measure, expected_speed
        complex(dp) :: hamiltonian
        integer :: n, it, iu, qu, edge, qx, branch, i, j

        info = -30
        if (.not. mixture_enabled) return
        if (.not. allocated(mixture_at_s)) return
        if (selected_species < 1 .or. selected_species > size(mixture_at_s)) return
        call validate_lorentz_geometry(info)
        if (info /= 0) return
        info = -31
        if (min(theta_order, velocity_order) < 2) return
        cache%coulomb_log = species_mixture%coulomb_log
        cache%mass_kg = mixture_at_s(selected_species)%mass_kg
        cache%charge_e = mixture_at_s(selected_species)%charge_e
        cache%density_m3 = mixture_at_s(selected_species)%density_m3
        cache%temperature_ev = mixture_at_s(selected_species)%temperature_ev
        cache%vth_m_s = 0.01_dp*vth
        expected_speed = sqrt(2.0_dp*collision_elementary_charge* &
                              cache%temperature_ev/cache%mass_kg)
        if (abs(cache%vth_m_s/expected_speed - 1.0_dp) > 2.0e-14_dp) return
        if (abs((qi/c)*10/(cache%charge_e*collision_elementary_charge) - &
                1.0_dp) > 2.0e-14_dp) return
        cache%A1 = A1
        cache%A2 = A2
        cache%beta = real(mph, dp)*q/((qi/c)*sign_theta*psi_pr)
        call initialize_conserving_basis(speed_steps, pitch_steps, umax, &
                                         etatp, etadt, cache%vth_m_s, cache%basis, info)
        if (info /= 0) return
        n = cache%basis%ndof
        allocate (cache%mass(n, n), cache%streaming(n, n), cache%precession(n, n))
        allocate (cache%hamiltonian(n, 3), cache%hamiltonian2(n, 3))
        allocate (cache%theta(theta_order), cache%theta_weight(theta_order))
        allocate (cache%jacobian_si(theta_order), cache%bmod_native(theta_order))
        allocate (nodes(velocity_order), weights(velocity_order))
        allocate (tail(velocity_order), wtail(velocity_order))
        allocate (values(n), gradient(3, n), dvalues(n))
        allocate (edges(4*pitch_steps + 3), mu(velocity_order), wm(velocity_order))
        call gauss_legendre(theta_order, cache%theta, cache%theta_weight)
        cache%theta = pi*(cache%theta + 1.0_dp)
        cache%theta_weight = pi*cache%theta_weight
        call gauss_legendre(velocity_order, nodes, weights)
        call gauss_gen_laguerre(velocity_order, 0.0_dp, tail, wtail)
        mu = nodes
        wm = weights
        cache%mass = 0.0_dp
        cache%streaming = 0.0_dp
        cache%precession = 0.0_dp
        cache%hamiltonian = (0.0_dp, 0.0_dp)
        cache%hamiltonian2 = (0.0_dp, 0.0_dp)
        do it = 1, theta_order
            call native_geometry(cache%theta(it), bmod, jac, ht, db)
            cache%jacobian_si(it) = jac
            cache%bmod_native(it) = bmod
            call pitch_edges(cache%basis, bmod, edges)
            do iu = 1, speed_steps + 1
                do qu = 1, velocity_order
                    if (iu <= speed_steps) then
                        span = umax/speed_steps
                        u = span*(real(iu - 1, dp) + (nodes(qu) + 1)/2)
                        wu = weights(qu)*span/2*u**2*exp(-u**2)
                    else
                        u = sqrt(umax**2 + tail(qu))
                        wu = wtail(qu)*u*exp(-umax**2)/2
                    end if
                    do edge = 1, size(edges) - 1
                        lo = edges(edge)
                        hi = edges(edge + 1)
                        if (hi <= lo) cycle
                        do qx = 1, velocity_order
                            xi = ((1 - mu(qx))*lo + (1 + mu(qx))*hi)/2
                            wxi = wm(qx)*(hi - lo)/2
                            velocity = cache%vth_m_s*u* &
                                       [sqrt(max(0.0_dp, 1 - xi**2)), 0.0_dp, xi]
                            call evaluate_conserving_basis(cache%basis, bmod, &
                                    velocity, values, gradient, branch, dvalues, ht, db)
                            call conserving_native_point(cache%theta(it), u, xi, &
                                                         omega, hamiltonian)
                            measure = 2*cache%density_m3/sqrt(pi)* &
                                      cache%theta_weight(it)*jac*wu*wxi
                            do j = 1, n
                                do i = 1, n
                                    cache%mass(i, j) = cache%mass(i, j) + &
                                                       measure*values(i)*values(j)
                                    cache%streaming(i, j) = cache%streaming(i, j) + &
                                                            measure*values(i)*dvalues(j)
                                    cache%precession(i, j) = cache%precession(i, j) + &
                                                       measure*values(i)*omega*values(j)
                                end do
                            end do
                            cache%hamiltonian(:, branch) = &
                               cache%hamiltonian(:, branch) + measure*values*hamiltonian
                            cache%hamiltonian2(:, branch) = &
                                cache%hamiltonian2(:, branch) + &
                                measure*values*u**2*hamiltonian
                        end do
                    end do
                end do
            end do
        end do
        info = 0
        if (.not. all(ieee_is_finite(cache%mass)) .or. &
            .not. all(ieee_is_finite(cache%streaming)) .or. &
            .not. all(ieee_is_finite(cache%precession)) .or. &
            .not. all(ieee_is_finite(real(cache%hamiltonian, dp))) .or. &
            .not. all(ieee_is_finite(aimag(cache%hamiltonian))) .or. &
            .not. all(ieee_is_finite(real(cache%hamiltonian2, dp))) .or. &
            .not. all(ieee_is_finite(aimag(cache%hamiltonian2)))) info = -32
    end subroutine prepare_conserving_species

    subroutine native_geometry(theta, bmod, jacobian_si, htheta_cm, dlogb_dtheta)
        use do_magfie_mod, only: do_magfie, s
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: bmod, jacobian_si, htheta_cm, dlogb_dtheta
        real(dp) :: sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        call do_magfie([s, 0.0_dp, theta], bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        jacobian_si = abs(sqrtg)*1.0e-6_dp
        htheta_cm = hctrvr(3)
        dlogb_dtheta = hder(3)
    end subroutine native_geometry

    subroutine conserving_native_point(theta, u, xi, omega, hamiltonian)
        use neort_transport, only: timestep_transport
        use neort_orbit, only: nvar, Om_tB_local_v2
        use neort_profiles, only: vth, Om_tE
        use driftorbit, only: mth, mph, sign_vpar, etatp, magdrift
        use do_magfie_mod, only: do_magfie, s, q
        use util, only: mi
        real(dp), intent(in) :: theta, u, xi
        real(dp), intent(out) :: omega
        complex(dp), intent(out) :: hamiltonian
        real(dp) :: y(nvar), derivative(nvar), bmod, sqrtg, hder(3)
        real(dp) :: hcovar(3), hctrvr(3), hcurl(3), eta, old_sign, speed
        integer :: old_mth
        call do_magfie([s, 0.0_dp, theta], bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        eta = (1.0_dp - xi**2)/bmod
        old_mth = mth
        old_sign = sign_vpar
        sign_vpar = sign(1.0_dp, xi)
        mth = -nint(mph*q)
        if (eta > etatp) mth = 0
        speed = u*vth
        y = 0.0_dp
        y(1) = theta
        y(2) = speed*xi
        call timestep_transport(speed, eta, nvar, 0.0_dp, y, derivative)
        hamiltonian = mi*speed**2/2*cmplx(derivative(3), derivative(4), dp)
        omega = real(mph, dp)*Om_tE
        if (magdrift) omega = omega + real(mph, dp)*speed**2* &
                              Om_tB_local_v2(eta, bmod, hder(1), hctrvr(3))
        sign_vpar = old_sign
        mth = old_mth
    end subroutine conserving_native_point

    subroutine pitch_edges(basis, bmod, edges)
        type(conserving_basis_t), intent(in) :: basis
        real(dp), intent(in) :: bmod
        real(dp), intent(out) :: edges(:)
        real(dp) :: eta, x, temporary
        integer :: k, n, j, positive_count
        n = 1
        edges(n) = 0.0_dp
        do k = 0, 2*basis%pitch_steps
            if (k <= basis%pitch_steps) then
                eta = basis%etatp*real(k, dp)/basis%pitch_steps
            else
                eta = basis%etatp + (basis%etadt - basis%etatp)* &
                      real(k - basis%pitch_steps, dp)/basis%pitch_steps
            end if
            x = sqrt(max(0.0_dp, 1.0_dp - eta*bmod))
            n = n + 1
            edges(n) = x
        end do
        positive_count = n
        do k = 2, n
            temporary = edges(k)
            j = k - 1
            do while (j >= 1)
                if (edges(j) <= temporary) exit
                edges(j + 1) = edges(j)
                j = j - 1
            end do
            edges(j + 1) = temporary
        end do
        edges(positive_count:2*positive_count - 1) = edges(1:positive_count)
        do k = 1, positive_count - 1
            edges(k) = -edges(2*positive_count - k)
        end do
    end subroutine pitch_edges

    subroutine assemble_conserving_collisions(caches, ntheta, quadrature, &
                                              matrix, info)
        type(conserving_species_t), intent(in) :: caches(:)
        integer, intent(in) :: ntheta
        type(landau_quadrature_t), intent(in) :: quadrature
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: mass(:), charge(:), density(:)
        integer, allocatable :: ndof(:)
        complex(dp), allocatable :: local_matrix(:, :)
        real(dp) :: local_b
        integer :: species, theta_index, n

        info = -33
        matrix = (0.0_dp, 0.0_dp)
        if (size(caches) == 0 .or. ntheta < 2) return
        n = 0
        allocate (mass(size(caches)), charge(size(caches)), density(size(caches)))
        allocate (ndof(size(caches)))
        do species = 1, size(caches)
            if (.not. allocated(caches(species)%theta)) return
            if (.not. allocated(caches(species)%theta_weight)) return
            if (.not. allocated(caches(species)%jacobian_si)) return
            if (.not. allocated(caches(species)%bmod_native)) return
            if (size(caches(species)%theta) /= ntheta) return
            if (size(caches(species)%theta_weight) /= ntheta) return
            if (size(caches(species)%jacobian_si) /= ntheta) return
            if (size(caches(species)%bmod_native) /= ntheta) return
            if (abs(caches(species)%coulomb_log/caches(1)%coulomb_log - &
                    1.0_dp) > 1.0e-13_dp) return
            if (maxval(abs(caches(species)%jacobian_si/ &
                           caches(1)%jacobian_si - 1.0_dp)) > 1.0e-13_dp) return
            if (maxval(abs(caches(species)%theta - caches(1)%theta)) > &
                1.0e-13_dp) return
            if (abs(caches(species)%temperature_ev/caches(1)%temperature_ev - &
                    1.0_dp) > 1.0e-13_dp) return
            if (maxval(abs(caches(species)%bmod_native/ &
                           caches(1)%bmod_native - 1.0_dp)) > 1.0e-13_dp) return
            mass(species) = caches(species)%mass_kg
            charge(species) = caches(species)%charge_e
            density(species) = caches(species)%density_m3
            ndof(species) = caches(species)%basis%ndof
            n = n + ndof(species)
        end do
        if (size(matrix, 1) /= n .or. size(matrix, 2) /= n) return
        allocate (local_matrix(n, n))
        do theta_index = 1, ntheta
            local_b = caches(1)%bmod_native(theta_index)
            call assemble_landau_projection(mass, charge, density, &
                    caches(1)%temperature_ev, caches(1)%coulomb_log, ndof, quadrature, &
                                            gradient, local_matrix)
            matrix = matrix + caches(1)%theta_weight(theta_index)* &
                     caches(1)%jacobian_si(theta_index)*local_matrix
        end do
        info = 0
    contains
        pure subroutine gradient(species_index, velocity, columns)
            integer, intent(in) :: species_index
            real(dp), intent(in) :: velocity(3)
            complex(dp), intent(out) :: columns(:, :)
            real(dp) :: values(size(columns, 2)), gradients(3, size(columns, 2))
            call evaluate_conserving_basis(caches(species_index)%basis, &
                                           local_b, velocity, values, gradients)
            columns = cmplx(gradients, 0.0_dp, dp)
        end subroutine gradient
    end subroutine assemble_conserving_collisions
end module neort_conserving_native
