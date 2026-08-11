module neort_orbit_classes
    ! Periodic magnetic-well inventory used by the thin-orbit path.
    !
    ! A single global Bmin/Bmax pair is sufficient only for a one-well field.
    ! General Boozer surfaces can contain several wells separated by unequal
    ! maxima.  Each well therefore gets its own minimum and the lower of its
    ! two neighboring barriers.  Passing orbits continue to use the global
    ! maximum because they traverse the complete periodic field line.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie
    use neort_orbit, only: th0
    use driftorbit, only: Bmin, Bmax, etatp, etadt
    use util, only: pi
    implicit none
    private

    integer, parameter :: NTH_CLASS_SCAN = 1000
    integer, parameter :: MAX_ORBIT_CLASSES = NTH_CLASS_SCAN/2

    type, public :: orbit_class_t
        real(dp) :: theta_min = 0.0_dp
        real(dp) :: Bmin = 0.0_dp
        real(dp) :: Bbarrier = 0.0_dp
        real(dp) :: etatp = 0.0_dp
        real(dp) :: etadt = 0.0_dp
    end type orbit_class_t

    type(orbit_class_t), allocatable, save :: classes(:)
    real(dp), save :: global_theta = 0.0_dp
    real(dp), save :: global_Bmin = 0.0_dp
    real(dp), save :: global_Bmax = 0.0_dp
    real(dp), save :: global_etatp = 0.0_dp
    real(dp), save :: global_etadt = 0.0_dp

    !$omp threadprivate (classes, global_theta, global_Bmin, global_Bmax)
    !$omp threadprivate (global_etatp, global_etadt)

    public :: init_orbit_classes
    public :: orbit_class_count
    public :: select_orbit_class
    public :: select_global_passing
    public :: find_periodic_extrema

contains

    subroutine init_orbit_classes(s)
        real(dp), intent(in) :: s
        real(dp) :: theta(NTH_CLASS_SCAN), bmod(NTH_CLASS_SCAN)
        logical :: minima(NTH_CLASS_SCAN), maxima(NTH_CLASS_SCAN)
        integer :: minima_index(MAX_ORBIT_CLASSES)
        integer :: maxima_index(MAX_ORBIT_CLASSES)
        integer :: nmin, nmax, k, minimum_count, maximum_count
        real(dp) :: b_refined, theta_refined
        real(dp) :: maxima_value(NTH_CLASS_SCAN)
        integer :: previous_max, next_max

        do k = 1, NTH_CLASS_SCAN
            theta(k) = -pi + 2.0_dp*pi*real(k - 1, dp)/real(NTH_CLASS_SCAN, dp)
            bmod(k) = field_strength(s, theta(k))
        end do

        call find_periodic_extrema(bmod, minima, maxima, nmin, nmax)
        minimum_count = 0
        maximum_count = 0
        do k = 1, NTH_CLASS_SCAN
            if (minima(k)) then
                minimum_count = minimum_count + 1
                if (minimum_count <= MAX_ORBIT_CLASSES) minima_index(minimum_count) = k
            end if
            if (maxima(k)) then
                maximum_count = maximum_count + 1
                if (maximum_count <= MAX_ORBIT_CLASSES) maxima_index(maximum_count) = k
            end if
        end do

        if (nmin == 0 .or. nmax == 0) then
            error stop "magnetic field has no resolvable periodic well"
        end if
        if (nmin > MAX_ORBIT_CLASSES .or. nmax > MAX_ORBIT_CLASSES) then
            error stop "magnetic field has too many periodic wells"
        end if

        if (allocated(classes)) deallocate(classes)
        allocate(classes(nmin))
        global_Bmin = huge(1.0_dp)
        global_Bmax = -huge(1.0_dp)
        maxima_value = 0.0_dp
        do k = 1, nmax
            call refine_extremum(s, theta, bmod, maxima_index(k), &
                theta_refined, b_refined)
            maxima_value(maxima_index(k)) = b_refined
            global_Bmax = max(global_Bmax, b_refined)
        end do
        do k = 1, nmin
            call refine_extremum(s, theta, bmod, minima_index(k), &
                theta_refined, b_refined)
            previous_max = neighboring_maximum(minima_index(k), -1, maxima)
            next_max = neighboring_maximum(minima_index(k), 1, maxima)
            classes(k)%theta_min = theta_refined
            classes(k)%Bmin = b_refined
            classes(k)%Bbarrier = min(maxima_value(previous_max), maxima_value(next_max))
            if (.not. ieee_is_finite(classes(k)%Bbarrier) .or. &
                    classes(k)%Bbarrier <= classes(k)%Bmin) then
                error stop "magnetic well has an invalid separatrix"
            end if
            classes(k)%etatp = 1.0_dp/classes(k)%Bbarrier
            classes(k)%etadt = 1.0_dp/classes(k)%Bmin
            if (b_refined < global_Bmin) then
                global_Bmin = b_refined
                global_theta = theta_refined
            end if
        end do

        global_etatp = 1.0_dp/global_Bmax
        global_etadt = 1.0_dp/global_Bmin
        call select_global_passing()
    end subroutine init_orbit_classes

    pure subroutine find_periodic_extrema(bmod, minima, maxima, nmin, nmax)
        real(dp), intent(in) :: bmod(:)
        logical, intent(out) :: minima(:), maxima(:)
        integer, intent(out) :: nmin, nmax
        integer :: k, left, right, n

        n = size(bmod)
        minima = .false.
        maxima = .false.
        nmin = 0
        nmax = 0
        if (size(minima) /= n .or. size(maxima) /= n .or. n < 3) return
        do k = 1, n
            left = periodic_index(k - 1, n)
            right = periodic_index(k + 1, n)
            if (bmod(k) <= bmod(left) .and. bmod(k) < bmod(right)) then
                nmin = nmin + 1
                minima(k) = .true.
            end if
            if (bmod(k) >= bmod(left) .and. bmod(k) > bmod(right)) then
                nmax = nmax + 1
                maxima(k) = .true.
            end if
        end do
    end subroutine find_periodic_extrema

    integer function orbit_class_count()
        if (allocated(classes)) then
            orbit_class_count = size(classes)
        else
            orbit_class_count = 0
        end if
    end function orbit_class_count

    subroutine select_orbit_class(index)
        integer, intent(in) :: index

        if (.not. allocated(classes)) error stop "orbit classes are not initialized"
        if (index < 1 .or. index > size(classes)) error stop "invalid orbit class"
        th0 = classes(index)%theta_min
        Bmin = classes(index)%Bmin
        Bmax = classes(index)%Bbarrier
        etatp = classes(index)%etatp
        etadt = classes(index)%etadt
    end subroutine select_orbit_class

    subroutine select_global_passing()
        th0 = global_theta
        Bmin = global_Bmin
        Bmax = global_Bmax
        etatp = global_etatp
        etadt = global_etadt
    end subroutine select_global_passing

    real(dp) function field_strength(s, theta) result(bmod)
        real(dp), intent(in) :: s, theta
        real(dp) :: sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [s, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        if (.not. ieee_is_finite(bmod)) error stop "magnetic field strength is not finite"
        if (bmod <= 0.0_dp) error stop "magnetic field strength is not positive"
    end function field_strength

    subroutine refine_extremum(s, theta, bmod, index, theta_refined, b_refined)
        real(dp), intent(in) :: s, theta(:), bmod(:)
        integer, intent(in) :: index
        real(dp), intent(out) :: theta_refined, b_refined
        integer :: left, right
        real(dp) :: h, denominator, offset, scale

        left = periodic_index(index - 1, size(bmod))
        right = periodic_index(index + 1, size(bmod))
        h = theta(2) - theta(1)
        denominator = bmod(left) - 2.0_dp*bmod(index) + bmod(right)
        scale = max(1.0_dp, abs(bmod(left)), abs(bmod(index)), abs(bmod(right)))
        offset = 0.0_dp
        if (abs(denominator) > 64.0_dp*epsilon(1.0_dp)*scale) then
            offset = h*(bmod(right) - bmod(left))/(2.0_dp*denominator)
            offset = max(-0.5_dp*h, min(0.5_dp*h, offset))
        end if
        theta_refined = theta(index) + offset
        b_refined = field_strength(s, theta_refined)
    end subroutine refine_extremum

    integer function neighboring_maximum(index, direction, maxima)
        integer, intent(in) :: index, direction
        logical, intent(in) :: maxima(:)
        integer :: candidate, step

        candidate = index
        do step = 1, size(maxima)
            candidate = periodic_index(candidate + direction, size(maxima))
            if (maxima(candidate)) then
                neighboring_maximum = candidate
                return
            end if
        end do
        error stop "magnetic well has no neighboring barrier"
    end function neighboring_maximum

    pure integer function periodic_index(index, n)
        integer, intent(in) :: index, n
        periodic_index = modulo(index - 1, n) + 1
    end function periodic_index

end module neort_orbit_classes
