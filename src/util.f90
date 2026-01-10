! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

module util
    use iso_fortran_env, only: dp => real64

    implicit none
    save

    complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)

    real(dp), parameter, public :: pi = 4 * atan(1.0_dp)

    real(dp), parameter, public :: &
        qe = 4.803204e-10_dp, &  ! elementary charge
        !me  = 9.109382e-28_dp,       & ! electron mass,
        mu = 1.660538e-24_dp, &  ! 1u
        c = 2.997925e+10_dp, &  ! speed of light
        !kb  = 1.381649e-16_dp,       & ! Boltzmann constant
        eV = 1.602176e-12_dp  ! 1 electron volt

    real(dp), public :: qi = 1.0_dp * qe, mi = 2.014_dp * mu

    ! Flux-surface dependent (set from plasma data at each s)
    !$omp threadprivate (qi, mi)

contains

    subroutine disp(str, val)
        character(*) :: str
        real(dp) :: val
        write (*, '("'//str//'" ES16.9E2)') val
    end subroutine disp

    function linspace(a, b, cnt)
        integer :: cnt
        real(dp) :: linspace(cnt)
        integer :: i
        real(dp) :: a, b, delta

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
        real(dp), allocatable, intent(out) :: data(:, :)
        real(dp) :: buf(ncol)
        integer :: lun, nrow, io, k

        open (unit=newunit(lun), file=filename, status="old")
        nrow = 0
        do
            read (lun, *, iostat=io) buf
            if (io /= 0) exit
            nrow = nrow + 1
        end do

        rewind (lun)

        allocate (data(nrow, ncol))

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

    subroutine cached_spline(x, x_prev, spl_coeff, y)
        use spline, only: spline_val_0

        real(dp), intent(in) :: x, x_prev
        real(dp), intent(in) :: spl_coeff(:, :, :)
        real(dp), intent(inout) :: y(:, :)

        real(dp), parameter :: tol = 1.0e-12_dp
        real(dp) :: spl_val(3)
        integer :: j

        if (abs(x - x_prev) < tol) return

        do j = 1, size(spl_coeff, 3)
            spl_val(:) = spline_val_0(spl_coeff(:, :, j), x)
            y(:, j) = spl_val
        end do
    end subroutine

    subroutine check_file(path, required, error_msg)
        use logger, only: error

        character(len=*), intent(in) :: path
        logical, intent(in) :: required
        character(len=*), intent(in) :: error_msg
        logical :: exists

        inquire(file=path, exist=exists)
        if (.not. exists .and. required) then
            call error(trim(path) // trim(error_msg))
        end if
    end subroutine check_file

    function files_exist(plasma_file, profile_file) result(exist)
        character(len=*), intent(in) :: plasma_file, profile_file
        logical :: plasma_exists, profile_exists, exist

        inquire(file=plasma_file, exist=plasma_exists)
        inquire(file=profile_file, exist=profile_exists)
        exist = plasma_exists .and. profile_exists
    end function files_exist

end module util
