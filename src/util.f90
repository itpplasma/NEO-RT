! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

module util
    implicit none
    save

    complex(8), parameter :: imun = (0d0, 1d0)

    real(8), parameter, public :: pi = 4*atan(1d0)

    real(8), parameter, public :: &
        qe = 4.803204d-10, & ! elementary charge
        !me  = 9.109382d-28,       & ! electron mass,
        mu = 1.660538d-24, & ! 1u
        c = 2.997925d+10, & ! speed of light
        !kb  = 1.381649d-16,       & ! Boltzmann constant
        eV = 1.602176d-12          ! 1 electron volt

    real(8), public :: qi = 1d0*qe, mi = 2.014d0*mu
contains

    subroutine disp(str, val)
        character(*) :: str
        real(8) :: val
        write (*, '("'//str//'" ES16.9E2)') val
    end subroutine disp

    function linspace(a, b, cnt)
        integer :: cnt
        real(8) :: linspace(cnt)
        integer :: i
        real(8) :: a, b, delta

        delta = (b - a)/(cnt - 1)
        linspace = a + delta*(/(i, i=0, cnt - 1)/)
    end function linspace

    ! From http://fortranwiki.org/fortran/show/newunit
    ! This is a simple function to search for an available unit.
    ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
    ! The UNIT value is returned by the function, and also by the optional
    ! argument. This allows the function to be used directly in an OPEN
    ! statement, and optionally save the result in a local variable.
    ! If no units are available, -1 is returned.
    integer function newunit(unit)
        integer, intent(out), optional :: unit
        ! local
        integer, parameter :: LUN_MIN = 100, LUN_MAX = 1000
        logical :: opened
        integer :: lun
        ! begin
        newunit = -1
        do lun = LUN_MIN, LUN_MAX
            inquire (unit=lun, opened=opened)
            if (.not. opened) then
                newunit = lun
                exit
            end if
        end do
        if (present(unit)) unit = newunit
    end function newunit

    subroutine readdata(filename, ncol, data)
        character(*), intent(in) :: filename
        integer, intent(in) :: ncol
        real(8), allocatable, intent(out) :: data(:, :)
        real(8) :: buf(ncol)
        integer :: lun, nrow, io, k

        open (unit=newunit(lun), file=filename, status="old")
        nrow = 0
        do
            read (lun, *, iostat=io) buf
            if (io /= 0) exit
            nrow = nrow + 1
        end do

        rewind (lun)

        allocate (data(nrow, 2))

        do k = 1, nrow
            read (lun, *) data(k, :)
        end do

        close (lun)
    end subroutine readdata

    subroutine clearfile(filename)
        character(*), intent(in) :: filename
        integer :: stat, lun

        open (unit=newunit(lun), iostat=stat, file=filename, status='old')
        if (stat == 0) close (lun, status='delete')
    end subroutine

    subroutine cached_spline(x, x_prev, spl_coeff, y_cached, y)
        use spline, only: spline_val_0

        real(8), intent(in) :: x, x_prev
        real(8), intent(in) :: spl_coeff(:, :, :)
        real(8), intent(in) :: y_cached(:, :)
        real(8), intent(out) :: y(:, :)

        real(8), parameter :: tol = 1.0d-12
        real(8) :: spl_val(3)
        integer :: j

        if (abs(x - x_prev) < tol) then
            y(:, :) = y_cached
            return
        end if

        do j = 1, size(spl_coeff, 3)
            spl_val(:) = spline_val_0(spl_coeff(:, :, j), x)
            y(:, j) = spl_val
        end do
    end subroutine

    subroutine cached_spline_value(x, x_prev, spl_coeff, val_cached, val)
        real(8), intent(in) :: x, x_prev
        real(8), intent(in) :: spl_coeff(:, :, :)
        real(8), intent(in) :: val_cached(:)
        real(8), intent(out) :: val(:)

        real(8) :: y_cached(3, size(val_cached, 1))
        real(8) :: y(3, size(val, 1))

        y_cached = 0d0
        y_cached(1, :) = val_cached
        call cached_spline(x, x_prev, spl_coeff, y_cached, y)
        val(:) = y(1, :)

    end subroutine cached_spline_value
end module util
