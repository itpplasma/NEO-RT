program test_poicut_midplane
   use field_eq_mod, only : allow_sol, nzet, zet, zmagaxis, rtf
   use field_sub, only : field_eq
   use field_mod, only : ampl, icall, iequil, ipert
   use input_files, only : convexfile
   use potato_input_mod, only : read_potato_input, npoicut, rho_pol_max
   use poicut_mod, only : npc, zpc_arr
   implicit none

   double precision, parameter :: cut_tolerance = 1.d-2
   character(len=1024) :: shifted_wall_file
   character(len=1024) :: argument
   double precision :: z_mid, theta
   double precision :: bmod, sqrtg
   double precision :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
   double precision :: br, bphi, bz
   double precision :: dbrdr, dbrdphi, dbrdz
   double precision :: dbphidr, dbphidphi, dbphidz
   double precision :: dbzdr, dbzdphi, dbzdz
   double precision :: x(3)
   integer :: output_unit, i
   double precision, parameter :: pi = 3.14159265358979323846d0
   integer, parameter :: nwall = 128
   double precision :: z_shift, wall_radius
   external :: find_poicut
   external :: magfie

   call get_command_argument(1, argument)
   shifted_wall_file = trim(argument)
   if (len_trim(shifted_wall_file) == 0) then
      error stop 'test_poicut_midplane requires a shifted wall output path'
   end if

   call read_potato_input('potato.in')
   call field_eq(200.d0, 0.d0, 0.d0, br, bphi, bz, dbrdr, dbrdphi, dbrdz, &
        dbphidr, dbphidphi, dbphidz, dbzdr, dbzdphi, dbzdz)
   z_shift = 0.75d0 * (zet(nzet) - zet(1))
   wall_radius = 2.d0 * max(rtf, max(abs(zet(1)), abs(zet(nzet))))
   open(newunit=output_unit, file=trim(shifted_wall_file), status='replace', &
        action='write')
   do i = 1, nwall
      theta = 2.d0 * pi * dble(i - 1) / dble(nwall)
      write(output_unit, *) rtf + wall_radius * cos(theta), &
           z_shift + wall_radius * sin(theta)
   end do
   close(output_unit)
   convexfile = trim(shifted_wall_file)
   icall = 1
   ipert = 0
   iequil = 1
   ampl = 1.d0

   x = [rtf, 0.d0, 0.d0]
   call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
   allow_sol = .false.
   zet = zet + z_shift
   zmagaxis = zmagaxis + z_shift
   z_mid = 0.5d0 * (zet(1) + zet(nzet))

   call find_poicut(rho_pol_max, npoicut)
   if (abs(zpc_arr(0) - z_mid) > cut_tolerance .or. &
       abs(zpc_arr(npc) - z_mid) > cut_tolerance) then
      error stop 'Poincare cut did not follow the shifted equilibrium midplane'
   end if

   print *, 'test_poicut_midplane PASSED'
end program test_poicut_midplane
