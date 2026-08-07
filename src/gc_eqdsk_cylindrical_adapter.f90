module neort_gc_eqdsk_cylindrical_adapter
    !! Direct real-space view of the initialized GEQDSK/libneo field.
    !! Inputs use (R,Z,phi); field_eq itself is called as (R,phi,Z).
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, R0, a
    use geoflux_coordinates, only: geoflux_to_cyl
    use field_sub, only: dpsidr, dpsidz, field_eq, psif
    use neort_gc_cylindrical_model, only: GC_CYL_EQUILIBRIUM_DOMAIN, &
        GC_CYL_INVALID_INPUT, GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t, &
        gc_cylindrical_field_t, make_gc_cylindrical_field_sample

    implicit none
    private

    type, extends(gc_cylindrical_field_t), public :: eqdsk_cylindrical_field_t
        real(dp) :: field_scale = 1.0_dp
        real(dp) :: domain_R_min = 0.0_dp
        real(dp) :: domain_R_max = 0.0_dp
        logical :: domain_initialized = .false.
    contains
        procedure :: evaluate => evaluate_eqdsk_cylindrical_field
    end type eqdsk_cylindrical_field_t

    public :: initialize_eqdsk_cylindrical_field
    public :: configure_eqdsk_cylindrical_field
    public :: map_eqdsk_flux_position

contains

    subroutine initialize_eqdsk_cylindrical_field(path, field_scale, field, status)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: field_scale
        type(eqdsk_cylindrical_field_t), intent(out) :: field
        integer, intent(out) :: status

        field = eqdsk_cylindrical_field_t()
        status = GC_CYL_INVALID_INPUT
        if (len_trim(path) == 0 .or. field_scale <= 0.0_dp) return
        inp_swi = 11
        call read_boozer_file(trim(path))
        call configure_eqdsk_cylindrical_field(field_scale, field, status)
    end subroutine initialize_eqdsk_cylindrical_field

    subroutine configure_eqdsk_cylindrical_field(field_scale, field, status)
        real(dp), intent(in) :: field_scale
        type(eqdsk_cylindrical_field_t), intent(out) :: field
        integer, intent(out) :: status

        field = eqdsk_cylindrical_field_t()
        field%field_scale = field_scale
        status = GC_CYL_INVALID_INPUT
        if (inp_swi /= 11 .or. field_scale <= 0.0_dp) return
        if (.not. all(ieee_is_finite([R0, a])) .or. a <= 0.0_dp) return
        field%domain_R_min = R0 - a
        field%domain_R_max = R0 + a
        if (field%domain_R_min <= 0.0_dp) return
        field%domain_initialized = .true.
        status = GC_CYL_SUCCESS
    end subroutine configure_eqdsk_cylindrical_field

    subroutine map_eqdsk_flux_position(surface, theta, phi, position, status)
        real(dp), intent(in) :: surface, theta, phi
        real(dp), intent(out) :: position(3)
        integer, intent(out) :: status

        real(dp) :: xgeo(3), xcyl(3)

        position = 0.0_dp
        status = GC_CYL_INVALID_INPUT
        if (inp_swi /= 11 .or. surface < 0.0_dp .or. surface > 1.0_dp) return
        if (.not. all(ieee_is_finite([surface, theta, phi]))) return
        xgeo = [surface, theta, phi]
        call geoflux_to_cyl(xgeo, xcyl)
        if (.not. all(ieee_is_finite(xcyl)) .or. xcyl(1) <= 0.0_dp) return
        ! geoflux uses (R,phi,Z); the cylindrical backend uses (R,Z,phi).
        position = [xcyl(1), xcyl(3), xcyl(2)]
        status = GC_CYL_SUCCESS
    end subroutine map_eqdsk_flux_position

    subroutine evaluate_eqdsk_cylindrical_field(self, position, sample, status)
        class(eqdsk_cylindrical_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status
        real(dp) :: br, bphi, bz
        real(dp) :: dbr_dr, dbr_dphi, dbr_dz
        real(dp) :: dbphi_dr, dbphi_dphi, dbphi_dz
        real(dp) :: dbz_dr, dbz_dphi, dbz_dz
        real(dp) :: b(3), db_dq(3, 3), grad_psi(3), psi
        integer :: local_status

        sample = gc_cylindrical_field_sample_t()
        status = GC_CYL_INVALID_INPUT
        if (self%field_scale <= 0.0_dp .or. position(1) <= 0.0_dp) return
        if (self%domain_initialized) then
            if (position(1) < self%domain_R_min .or. &
                    position(1) > self%domain_R_max) then
                status = GC_CYL_EQUILIBRIUM_DOMAIN
                return
            end if
        end if
        call field_eq(position(1), position(3), position(2), br, bphi, bz, &
            dbr_dr, dbr_dphi, dbr_dz, dbphi_dr, dbphi_dphi, dbphi_dz, &
            dbz_dr, dbz_dphi, dbz_dz)
        b = self%field_scale*[br, bphi, bz]
        db_dq(:, 1) = self%field_scale*[dbr_dr, dbphi_dr, dbz_dr]
        db_dq(:, 2) = self%field_scale*[dbr_dphi, dbphi_dphi, dbz_dphi] &
            /position(1)
        db_dq(:, 3) = self%field_scale*[dbr_dz, dbphi_dz, dbz_dz]
        psi = self%field_scale*psif
        grad_psi = self%field_scale*[dpsidr, 0.0_dp, dpsidz]
        if (.not. all(ieee_is_finite([b, db_dq, psi, grad_psi]))) return
        call make_gc_cylindrical_field_sample(position(1), b, db_dq, psi, &
            grad_psi, sample, local_status)
        status = local_status
    end subroutine evaluate_eqdsk_cylindrical_field

end module neort_gc_eqdsk_cylindrical_adapter
