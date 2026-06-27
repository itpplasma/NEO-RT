module wall_loss_mod
    use libneo_kinds, only: dp
    implicit none
    private

    public :: load_wall, outside_wall, wall_loaded

    logical :: wall_loaded = .false.
    real(dp), allocatable :: rwall(:), zwall(:)
    real(dp) :: rc, zc, r_in2, r_out2

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

        call set_bounding_circles()
        wall_loaded = .true.
    end subroutine load_wall

    subroutine set_bounding_circles()
        integer :: i, j, n
        real(dp) :: dx, dy, seg2, dist2

        n = size(rwall)
        rc = sum(rwall) / n
        zc = sum(zwall) / n
        r_in2 = huge(1.0_dp)
        r_out2 = 0.0_dp
        j = n
        do i = 1, n
            dist2 = (rwall(i) - rc)**2 + (zwall(i) - zc)**2
            if (dist2 > r_out2) r_out2 = dist2
            dx = rwall(i) - rwall(j)
            dy = zwall(i) - zwall(j)
            seg2 = dx*dx + dy*dy
            if (seg2 > 0.0_dp) then
                dist2 = (dx*(zc - zwall(j)) - dy*(rc - rwall(j)))**2 / seg2
                if (dist2 < r_in2) r_in2 = dist2
            end if
            j = i
        end do
    end subroutine set_bounding_circles

    logical function outside_wall(r, z)
        real(dp), intent(in) :: r, z
        integer :: i, j, n
        real(dp) :: d2
        logical :: inside

        outside_wall = .false.
        if (.not. wall_loaded) return

        d2 = (r - rc)**2 + (z - zc)**2
        if (d2 <= r_in2) return
        if (d2 >= r_out2) then
            outside_wall = .true.
            return
        end if

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
        outside_wall = .not. inside
    end function outside_wall

end module wall_loss_mod
