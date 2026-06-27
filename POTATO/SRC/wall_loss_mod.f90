module wall_loss_mod
    use libneo_kinds, only: dp
    implicit none
    private

    public :: load_wall, outside_wall, wall_loaded

    logical :: wall_loaded = .false.
    real(dp), allocatable :: rwall(:), zwall(:)

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
        wall_loaded = .true.
    end subroutine load_wall

    logical function outside_wall(r, z)
        real(dp), intent(in) :: r, z
        integer :: i, j, n
        logical :: inside

        outside_wall = .false.
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
        outside_wall = .not. inside
    end function outside_wall

end module wall_loss_mod
