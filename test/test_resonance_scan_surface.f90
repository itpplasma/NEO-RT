program test_resonance_scan_surface
    use iso_fortran_env, only: real64
    use diag_resonance_scan, only: read_surface_table
    implicit none

    real(real64), allocatable :: surfaces(:)
    integer :: u

    open(newunit=u, file="surface_list_test.in", status="replace", action="write")
    write(u, *) 0.25_real64
    write(u, *) 0.50_real64
    write(u, *) 0.75_real64
    close(u)

    call read_surface_table("surface_list_test.in", surfaces)
    if (size(surfaces) /= 3) error stop "surface-list count mismatch"
    if (maxval(abs(surfaces - [0.25_real64, 0.50_real64, 0.75_real64])) > 1.0e-15_real64) then
        error stop "surface-list values mismatch"
    end if
    if (any(surfaces(2:) <= surfaces(:size(surfaces) - 1))) error stop "surface-list ordering mismatch"

    open(newunit=u, file="surface_list_test.in", status="old", action="read")
    close(u, status="delete")
end program test_resonance_scan_surface
