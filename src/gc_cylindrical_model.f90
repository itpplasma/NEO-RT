module neort_gc_cylindrical_model
    !! Data-oriented real-space guiding-center contracts.
    !!
    !! Positions are ordered (R,Z,phi), while vectors use the orthonormal
    !! cylindrical basis (e_R,e_phi,e_Z).  R, Z, B, psi, charge, and mass use
    !! the CGS units already used by the direct GEQDSK/libneo path.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_perpendicular_invariant, only: &
        gc_buchholz_jk_from_mu_phys, gc_potato_jperp_from_mu_phys

    implicit none
    private

    real(dp), parameter, public :: GC_CYL_DEFAULT_C = 2.99792458e10_dp

    integer, parameter, public :: GC_CYL_SUCCESS = 0
    integer, parameter, public :: GC_CYL_INVALID_INPUT = 1
    integer, parameter, public :: GC_CYL_SINGULAR_BSTAR = 2
    integer, parameter, public :: GC_CYL_FIELD_ERROR = 10
    integer, parameter, public :: GC_CYL_POTENTIAL_ERROR = 11
    integer, parameter, public :: GC_CYL_STATE_ERROR = 12
    integer, parameter, public :: GC_CYL_START_ERROR = 13
    integer, parameter, public :: GC_CYL_INTEGRATOR_ERROR = 14
    integer, parameter, public :: GC_CYL_NO_RETURN = 15
    integer, parameter, public :: GC_CYL_INVARIANT_ERROR = 16
    integer, parameter, public :: GC_CYL_WALL_LOSS = 20
    integer, parameter, public :: GC_CYL_WALL_ERROR = 21
    integer, parameter, public :: GC_CYL_PERTURBATION_ERROR = 22
    integer, parameter, public :: GC_CYL_EQUILIBRIUM_DOMAIN = 30
    integer, parameter, public :: GC_CYL_MEASURE_UNAVAILABLE = 31

    integer, parameter, public :: GC_CYL_SECTION_PLANE = 1
    integer, parameter, public :: GC_CYL_SECTION_R = 2
    integer, parameter, public :: GC_CYL_SECTION_Z = 3
    integer, parameter, public :: GC_CYL_SECTION_PHI = 4

    type, public :: gc_cylindrical_state_t
        !! mu is the canonical physical magnetic moment mu_phys.  It is not
        !! Buchholz J_K and is not POTATO's normalized p^2(1-xi^2)/B.
        real(dp) :: R = 0.0_dp
        real(dp) :: Z = 0.0_dp
        real(dp) :: phi = 0.0_dp
        real(dp) :: p_parallel = 0.0_dp
        real(dp) :: mu = 0.0_dp
    end type gc_cylindrical_state_t

    type, public :: gc_cylindrical_field_sample_t
        real(dp) :: radius = 0.0_dp
        real(dp) :: b(3) = 0.0_dp
        real(dp) :: bhat(3) = 0.0_dp
        real(dp) :: bmod = 0.0_dp
        real(dp) :: db_dq(3, 3) = 0.0_dp
        real(dp) :: grad_b(3) = 0.0_dp
        real(dp) :: curl_bhat(3) = 0.0_dp
        real(dp) :: psi = 0.0_dp
        real(dp) :: grad_psi(3) = 0.0_dp
    end type gc_cylindrical_field_sample_t

    type, public :: gc_cylindrical_invariants_t
        !! magnetic_moment is mu_phys in the same units as state%mu.
        real(dp) :: energy = 0.0_dp
        real(dp) :: magnetic_moment = 0.0_dp
        real(dp) :: canonical_toroidal_momentum = 0.0_dp
    end type gc_cylindrical_invariants_t

    type, public :: gc_cylindrical_allowed_component_t
        integer :: component_id = 0
        integer :: sigma = 0
        real(dp) :: x_begin = 0.0_dp
        real(dp) :: x_end = 0.0_dp
        real(dp) :: canonical_begin = 0.0_dp
        real(dp) :: canonical_end = 0.0_dp
        real(dp) :: canonical_measure = 0.0_dp
        logical :: lower_root = .false.
        logical :: upper_root = .false.
    end type gc_cylindrical_allowed_component_t

    type, public :: gc_cylindrical_section_t
        integer :: kind = GC_CYL_SECTION_PLANE
        real(dp) :: point(3) = 0.0_dp
        real(dp) :: normal(2) = [0.0_dp, 1.0_dp]
        integer :: direction = 0
        integer :: winding = 0
        ! A physical R/Z cut is a Poincare section, not merely a zero of a
        ! local chart coordinate.  return_crossings records how many
        ! same-oriented cut crossings make one requested return.  It is kept
        ! explicit so a caller cannot silently use an opposite crossing or a
        ! synthetic t=0 root as a full field-line return.
        integer :: return_crossings = 1
        integer :: frequency_winding = 1
        real(dp) :: fieldline_q = 0.0_dp
    end type gc_cylindrical_section_t

    type, abstract, public :: gc_cylindrical_field_t
    contains
        procedure(gc_cylindrical_field_evaluate_i), deferred :: evaluate
    end type gc_cylindrical_field_t

    type, abstract, public :: gc_cylindrical_potential_t
    contains
        procedure(gc_cylindrical_potential_evaluate_i), deferred :: evaluate
    end type gc_cylindrical_potential_t

    type, abstract, public :: gc_cylindrical_wall_t
    contains
        procedure(gc_cylindrical_wall_evaluate_i), deferred :: evaluate
    end type gc_cylindrical_wall_t

    type, extends(gc_cylindrical_potential_t), public :: gc_cylindrical_zero_potential_t
    contains
        procedure :: evaluate => evaluate_zero_potential
    end type gc_cylindrical_zero_potential_t

    type, extends(gc_cylindrical_potential_t), public :: gc_cylindrical_linear_flux_potential_t
        real(dp) :: coefficient = 0.0_dp
        real(dp) :: psi_reference = 0.0_dp
    contains
        procedure :: evaluate => evaluate_linear_flux_potential
    end type gc_cylindrical_linear_flux_potential_t

    type, extends(gc_cylindrical_wall_t), public :: gc_cylindrical_polygon_wall_t
        real(dp), allocatable :: vertices(:, :)
        logical :: initialized = .false.
    contains
        procedure :: set_vertices => set_polygon_vertices
        procedure :: evaluate => evaluate_polygon_wall
    end type gc_cylindrical_polygon_wall_t

    abstract interface
        subroutine gc_cylindrical_allowed_value_i(x, vparallel_squared, &
                dvparallel_squared_dx, psi_star, dpsi_star_dx, status)
            import :: dp
            real(dp), intent(in) :: x
            real(dp), intent(out) :: vparallel_squared
            real(dp), intent(out) :: dvparallel_squared_dx
            real(dp), intent(out) :: psi_star
            real(dp), intent(out) :: dpsi_star_dx
            integer, intent(out) :: status
        end subroutine gc_cylindrical_allowed_value_i

        subroutine gc_cylindrical_field_evaluate_i(self, position, sample, status)
            import :: dp, gc_cylindrical_field_t, gc_cylindrical_field_sample_t
            class(gc_cylindrical_field_t), intent(in) :: self
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_field_sample_t), intent(out) :: sample
            integer, intent(out) :: status
        end subroutine gc_cylindrical_field_evaluate_i

        subroutine gc_cylindrical_potential_evaluate_i(self, position, field, &
                potential, gradient, status)
            import :: dp, gc_cylindrical_field_sample_t, &
                gc_cylindrical_potential_t
            class(gc_cylindrical_potential_t), intent(in) :: self
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_field_sample_t), intent(in) :: field
            real(dp), intent(out) :: potential, gradient(3)
            integer, intent(out) :: status
        end subroutine gc_cylindrical_potential_evaluate_i

        subroutine gc_cylindrical_wall_evaluate_i(self, position, distance, status)
            import :: dp, gc_cylindrical_wall_t
            class(gc_cylindrical_wall_t), intent(in) :: self
            real(dp), intent(in) :: position(3)
            real(dp), intent(out) :: distance
            integer, intent(out) :: status
        end subroutine gc_cylindrical_wall_evaluate_i
    end interface

    public :: make_gc_cylindrical_field_sample
    public :: gc_cylindrical_section_value
    public :: gc_cylindrical_section_rate
    public :: invariants_from_cylindrical_state
    public :: gc_cylindrical_invariant_residuals
    public :: state_from_cylindrical_invariants
    public :: gc_cylindrical_vparallel_squared
    public :: canonical_toroidal_momentum_from_state
    public :: canonical_flux_from_state
    public :: energy_from_state
    public :: gc_buchholz_jk_from_state
    public :: gc_potato_jperp_from_state
    public :: gc_cylindrical_allowed_value_i

contains

    pure subroutine make_gc_cylindrical_field_sample(radius, b, db_dq, psi, &
            grad_psi, sample, status)
        real(dp), intent(in) :: radius, b(3), db_dq(3, 3), psi, grad_psi(3)
        type(gc_cylindrical_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        real(dp) :: dbhat_dq(3, 3), curl_bhat(3)
        integer :: i, j

        sample = gc_cylindrical_field_sample_t()
        status = GC_CYL_INVALID_INPUT
        if (radius <= 0.0_dp) return
        if (.not. all(ieee_is_finite([radius, b, db_dq, psi, grad_psi]))) return
        sample%radius = radius
        sample%b = b
        sample%bmod = sqrt(dot_product(b, b))
        if (sample%bmod <= tiny(sample%bmod)) return
        sample%bhat = b/sample%bmod
        sample%db_dq = db_dq
        do j = 1, 3
            sample%grad_b(j) = dot_product(sample%bhat, db_dq(:, j))
            do i = 1, 3
                dbhat_dq(i, j) = (db_dq(i, j) &
                    -sample%bhat(i)*sample%grad_b(j))/sample%bmod
            end do
        end do
        curl_bhat(1) = dbhat_dq(3, 2) - dbhat_dq(2, 3)
        curl_bhat(2) = dbhat_dq(1, 3) - dbhat_dq(3, 1)
        curl_bhat(3) = dbhat_dq(2, 1) + sample%bhat(2)/radius &
            -dbhat_dq(1, 2)
        sample%curl_bhat = curl_bhat
        sample%psi = psi
        sample%grad_psi = grad_psi
        status = GC_CYL_SUCCESS
    end subroutine make_gc_cylindrical_field_sample

    pure function gc_cylindrical_section_value(section, state) result(value)
        type(gc_cylindrical_section_t), intent(in) :: section
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp) :: value

        select case (section%kind)
        case (GC_CYL_SECTION_PLANE)
            value = section%normal(1)*(state%R - section%point(1)) &
                +section%normal(2)*(state%Z - section%point(2))
        case (GC_CYL_SECTION_R)
            value = state%R - section%point(1)
        case (GC_CYL_SECTION_Z)
            value = state%Z - section%point(2)
        case (GC_CYL_SECTION_PHI)
            value = state%phi - section%point(3) &
                -2.0_dp*acos(-1.0_dp)*real(section%winding, dp)
        case default
            value = huge(1.0_dp)
        end select
    end function gc_cylindrical_section_value

    pure function gc_cylindrical_section_rate(section, derivative) result(rate)
        type(gc_cylindrical_section_t), intent(in) :: section
        real(dp), intent(in) :: derivative(5)
        real(dp) :: rate

        select case (section%kind)
        case (GC_CYL_SECTION_PLANE)
            rate = section%normal(1)*derivative(1) &
                +section%normal(2)*derivative(2)
        case (GC_CYL_SECTION_R)
            rate = derivative(1)
        case (GC_CYL_SECTION_Z)
            rate = derivative(2)
        case (GC_CYL_SECTION_PHI)
            rate = derivative(3)
        case default
            rate = 0.0_dp
        end select
    end function gc_cylindrical_section_rate

    pure subroutine invariants_from_cylindrical_state(state, field, potential, &
            mass, charge, c_light, invariants, status)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: potential, mass, charge, c_light
        type(gc_cylindrical_invariants_t), intent(out) :: invariants
        integer, intent(out) :: status

        invariants = gc_cylindrical_invariants_t()
        status = GC_CYL_INVALID_INPUT
        if (mass <= 0.0_dp .or. charge == 0.0_dp .or. c_light <= 0.0_dp) return
        if (state%mu < 0.0_dp .or. field%bmod <= 0.0_dp) return
        invariants%energy = energy_from_state(state, field, potential, charge, mass)
        invariants%magnetic_moment = state%mu
        invariants%canonical_toroidal_momentum = &
            canonical_toroidal_momentum_from_state(state, field, charge, c_light)
        if (.not. all(ieee_is_finite([invariants%energy, &
            invariants%magnetic_moment, invariants%canonical_toroidal_momentum]))) return
        status = GC_CYL_SUCCESS
    end subroutine invariants_from_cylindrical_state

    pure subroutine gc_cylindrical_invariant_residuals(state, field, potential, &
            mass, charge, c_light, invariants, energy_residual, &
            magnetic_moment_residual, canonical_momentum_residual, status)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: potential, mass, charge, c_light
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(out) :: energy_residual, magnetic_moment_residual
        real(dp), intent(out) :: canonical_momentum_residual
        integer, intent(out) :: status

        energy_residual = 0.0_dp
        magnetic_moment_residual = 0.0_dp
        canonical_momentum_residual = 0.0_dp
        status = GC_CYL_INVALID_INPUT
        if (mass <= 0.0_dp .or. charge == 0.0_dp .or. c_light <= 0.0_dp) return
        if (field%bmod <= 0.0_dp .or. state%mu < 0.0_dp) return
        energy_residual = energy_from_state(state, field, potential, charge, mass) &
            -invariants%energy
        magnetic_moment_residual = state%mu - invariants%magnetic_moment
        canonical_momentum_residual = &
            canonical_toroidal_momentum_from_state(state, field, charge, c_light) &
            -invariants%canonical_toroidal_momentum
        if (.not. all(ieee_is_finite([energy_residual, magnetic_moment_residual, &
            canonical_momentum_residual]))) return
        status = GC_CYL_SUCCESS
    end subroutine gc_cylindrical_invariant_residuals

    pure subroutine state_from_cylindrical_invariants(position, field, potential, &
            invariants, parallel_sign, mass, charge, c_light, state, status)
        real(dp), intent(in) :: position(3), potential
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        integer, intent(in) :: parallel_sign
        real(dp), intent(in) :: mass, charge, c_light
        type(gc_cylindrical_state_t), intent(out) :: state
        integer, intent(out) :: status

        real(dp) :: denominator, p_parallel, kinetic, p2_expected, mismatch

        state = gc_cylindrical_state_t()
        state%R = position(1)
        state%Z = position(2)
        state%phi = position(3)
        status = GC_CYL_START_ERROR
        if (parallel_sign == 0) return
        if (mass <= 0.0_dp .or. charge == 0.0_dp .or. c_light <= 0.0_dp) return
        if (field%bmod <= 0.0_dp .or. invariants%magnetic_moment < 0.0_dp) return
        denominator = state%R*field%bhat(2)
        if (abs(denominator) <= tiny(denominator)) return
        p_parallel = (invariants%canonical_toroidal_momentum &
            -charge*field%psi/c_light)/denominator
        ! A turning point is a valid launch.  Its sign is inherited from the
        ! requested branch only when the momentum is nonzero; rejecting zero
        ! here made the fixed-invariant return map unable to start at a
        ! p_parallel=0 turning point.
        if (p_parallel*real(parallel_sign, dp) < 0.0_dp) return
        kinetic = invariants%energy - invariants%magnetic_moment*field%bmod &
            -charge*potential
        if (kinetic < 0.0_dp) return
        p2_expected = 2.0_dp*mass*kinetic
        mismatch = abs(p_parallel**2 - p2_expected)
        if (mismatch > 1.0e-9_dp*max(abs(p2_expected), &
            tiny(p2_expected))) return
        state%p_parallel = p_parallel
        state%mu = invariants%magnetic_moment
        status = GC_CYL_SUCCESS
    end subroutine state_from_cylindrical_invariants

    pure function canonical_toroidal_momentum_from_state(state, field, charge, &
            c_light) result(value)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: charge, c_light
        real(dp) :: value

        value = charge*field%psi/c_light &
            +state%p_parallel*state%R*field%bhat(2)
    end function canonical_toroidal_momentum_from_state

    pure function canonical_flux_from_state(state, field, charge, c_light) &
            result(value)
        !! psi_star has the canonical sign P_phi = (q/c) psi_star.
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: charge, c_light
        real(dp) :: value

        value = field%psi + c_light/charge*state%p_parallel*state%R &
            *field%bhat(2)
    end function canonical_flux_from_state

    pure function gc_buchholz_jk_from_state(state, charge, c_light) result(value)
        !! Convert the canonical state magnetic moment to Buchholz J_K.
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp), intent(in) :: charge, c_light
        real(dp) :: value

        value = gc_buchholz_jk_from_mu_phys(state%mu, charge, c_light)
    end function gc_buchholz_jk_from_state

    pure function gc_potato_jperp_from_state(state, mass, bmod, v0) &
            result(value)
        !! Convert the canonical state magnetic moment to POTATO's
        !! normalized p^2(1-xi^2)/B.  v0 and bmod are explicit by design.
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp), intent(in) :: mass, bmod, v0
        real(dp) :: value

        value = gc_potato_jperp_from_mu_phys(state%mu, mass, bmod, v0)
    end function gc_potato_jperp_from_state

    pure function energy_from_state(state, field, potential, charge, mass) result(value)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: potential, charge, mass
        real(dp) :: value

        value = state%p_parallel**2/(2.0_dp*mass) &
            +state%mu*field%bmod + charge*potential
    end function energy_from_state

    pure function gc_cylindrical_vparallel_squared(energy, magnetic_moment, &
            potential, bmod, mass, charge) result(value)
        real(dp), intent(in) :: energy, magnetic_moment, potential, bmod
        real(dp), intent(in) :: mass, charge
        real(dp) :: value

        value = 2.0_dp*(energy - magnetic_moment*bmod - charge*potential)/mass
    end function gc_cylindrical_vparallel_squared

    pure subroutine evaluate_zero_potential(self, position, field, potential, &
            gradient, status)
        class(gc_cylindrical_zero_potential_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(out) :: potential, gradient(3)
        integer, intent(out) :: status

        associate (unused_self => self, unused_position => position, &
                unused_field => field)
        end associate
        potential = 0.0_dp
        gradient = 0.0_dp
        status = GC_CYL_SUCCESS
    end subroutine evaluate_zero_potential

    pure subroutine evaluate_linear_flux_potential(self, position, field, &
            potential, gradient, status)
        class(gc_cylindrical_linear_flux_potential_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(out) :: potential, gradient(3)
        integer, intent(out) :: status

        associate (unused_position => position)
        end associate
        potential = self%coefficient*(field%psi - self%psi_reference)
        gradient = self%coefficient*field%grad_psi
        status = GC_CYL_SUCCESS
    end subroutine evaluate_linear_flux_potential

    subroutine set_polygon_vertices(self, vertices, status)
        class(gc_cylindrical_polygon_wall_t), intent(inout) :: self
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(out) :: status

        self%initialized = .false.
        status = GC_CYL_INVALID_INPUT
        if (size(vertices, 1) /= 2 .or. size(vertices, 2) < 3) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (allocated(self%vertices)) deallocate(self%vertices)
        allocate(self%vertices(2, size(vertices, 2)))
        self%vertices = vertices
        self%initialized = .true.
        status = GC_CYL_SUCCESS
    end subroutine set_polygon_vertices

    subroutine evaluate_polygon_wall(self, position, distance, status)
        class(gc_cylindrical_polygon_wall_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: distance
        integer, intent(out) :: status

        integer :: i, j, n
        logical :: inside
        real(dp) :: edge_distance

        distance = -huge(1.0_dp)
        status = GC_CYL_WALL_ERROR
        if (.not. self%initialized) return
        n = size(self%vertices, 2)
        inside = .false.
        j = n
        distance = huge(1.0_dp)
        do i = 1, n
            edge_distance = point_segment_distance(position(1:2), &
                self%vertices(:, j), self%vertices(:, i))
            distance = min(distance, edge_distance)
            call toggle_polygon_crossing(position(1:2), self%vertices(:, j), &
                self%vertices(:, i), inside)
            j = i
        end do
        if (.not. inside) distance = -distance
        status = GC_CYL_SUCCESS
    end subroutine evaluate_polygon_wall

    pure subroutine toggle_polygon_crossing(point, first, second, inside)
        real(dp), intent(in) :: point(2), first(2), second(2)
        logical, intent(inout) :: inside

        real(dp) :: denominator, intersection

        if ((first(2) > point(2)) .eqv. (second(2) > point(2))) return
        denominator = second(2) - first(2)
        if (abs(denominator) <= tiny(denominator)) return
        intersection = first(1) + (point(2) - first(2)) &
            *(second(1) - first(1))/denominator
        if (point(1) < intersection) inside = .not. inside
    end subroutine toggle_polygon_crossing

    pure function point_segment_distance(point, first, second) result(distance)
        real(dp), intent(in) :: point(2), first(2), second(2)
        real(dp) :: distance
        real(dp) :: edge(2), offset(2), denominator, fraction

        edge = second - first
        offset = point - first
        denominator = dot_product(edge, edge)
        fraction = 0.0_dp
        if (denominator > 0.0_dp) then
            fraction = dot_product(offset, edge)/denominator
            fraction = max(0.0_dp, min(1.0_dp, fraction))
        end if
        distance = sqrt(dot_product(point - (first + fraction*edge), &
            point - (first + fraction*edge)))
    end function point_segment_distance

end module neort_gc_cylindrical_model
