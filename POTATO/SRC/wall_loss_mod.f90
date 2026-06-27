module wall_loss_mod
    use libneo_kinds, only: dp
    implicit none
    private

    public :: load_wall, outside_wall, outside_wall_exact, wall_loaded
    public :: fastpath_resolves_circle, fastpath_resolves_ellipse

    logical :: wall_loaded = .false.
    real(dp), allocatable :: rwall(:), zwall(:)
    real(dp) :: rc, zc, sz
    real(dp) :: r_in2_c, r_out2_c, r_in2_e, r_out2_e

contains

    subroutine load_wall(filename)
        character(len=*), intent(in) :: filename
        integer :: iunit, ios, n
        real(dp) :: r, z
        logical :: exists

        inquire (file=filename, exist=exists)
        if (.not. exists) return
        open (newunit=iunit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) return

        n = 0
        do
            read (iunit, *, iostat=ios) r, z
            if (ios /= 0) exit
            n = n + 1
        end do
        if (n < 3) then
            close (iunit)
            return
        end if

        allocate (rwall(n), zwall(n))
        rewind (iunit)
        n = 0
        do
            read (iunit, *, iostat=ios) r, z
            if (ios /= 0) exit
            n = n + 1
            rwall(n) = r
            zwall(n) = z
        end do
        close (iunit)

        call set_bounding_shapes()
        wall_loaded = .true.
    end subroutine load_wall

    subroutine set_bounding_shapes()
        integer :: n
        real(dp) :: span_r, span_z

        n = size(rwall)
        rc = sum(rwall) / n
        zc = sum(zwall) / n

        span_r = maxval(rwall) - minval(rwall)
        span_z = maxval(zwall) - minval(zwall)
        if (span_z > 0.0_dp) then
            sz = span_r / span_z
        else
            sz = 1.0_dp
        end if

        call bounding_circles(1.0_dp, r_in2_c, r_out2_c)
        call bounding_circles(sz, r_in2_e, r_out2_e)
    end subroutine set_bounding_shapes

    subroutine bounding_circles(scale_z, r_in2, r_out2)
        real(dp), intent(in) :: scale_z
        real(dp), intent(out) :: r_in2, r_out2
        integer :: i, j, n
        real(dp) :: dx, dy, seg2, dist2

        n = size(rwall)
        r_in2 = huge(1.0_dp)
        r_out2 = 0.0_dp
        j = n
        do i = 1, n
            dist2 = (rwall(i) - rc)**2 + (scale_z*(zwall(i) - zc))**2
            if (dist2 > r_out2) r_out2 = dist2
            dx = rwall(i) - rwall(j)
            dy = scale_z*(zwall(i) - zwall(j))
            seg2 = dx*dx + dy*dy
            if (seg2 > 0.0_dp) then
                dist2 = (dx*scale_z*(zc - zwall(j)) - dy*(rc - rwall(j)))**2 / seg2
                if (dist2 < r_in2) r_in2 = dist2
            end if
            j = i
        end do
    end subroutine bounding_circles

    logical function outside_wall(r, z)
        real(dp), intent(in) :: r, z
        real(dp) :: d2

        outside_wall = .false.
        if (.not. wall_loaded) return

        d2 = (r - rc)**2 + (sz*(z - zc))**2
        if (d2 <= r_in2_e) return
        if (d2 >= r_out2_e) then
            outside_wall = .true.
            return
        end if
        outside_wall = outside_wall_exact(r, z)
    end function outside_wall

    logical function outside_wall_exact(r, z)
        real(dp), intent(in) :: r, z
        integer :: i, j, n
        logical :: inside

        outside_wall_exact = .false.
        if (.not. wall_loaded) return

        n = size(rwall)
        inside = .false.
        j = n
        do i = 1, n
            if ((zwall(i) > z) .neqv. (zwall(j) > z)) then
                if (r < (rwall(j) - rwall(i)) * (z - zwall(i)) &
                    / (zwall(j) - zwall(i)) + rwall(i)) then
                    inside = .not. inside
                end if
            end if
            j = i
        end do
        outside_wall_exact = .not. inside
    end function outside_wall_exact

    logical function fastpath_resolves_circle(r, z)
        real(dp), intent(in) :: r, z
        real(dp) :: d2

        d2 = (r - rc)**2 + (z - zc)**2
        fastpath_resolves_circle = d2 <= r_in2_c .or. d2 >= r_out2_c
    end function fastpath_resolves_circle

    logical function fastpath_resolves_ellipse(r, z)
        real(dp), intent(in) :: r, z
        real(dp) :: d2

        d2 = (r - rc)**2 + (sz*(z - zc))**2
        fastpath_resolves_ellipse = d2 <= r_in2_e .or. d2 >= r_out2_e
    end function fastpath_resolves_ellipse

end module wall_loss_mod
