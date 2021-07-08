! 
  module signal_array_mod
!
    implicit none
!
    integer, dimension(:,:), allocatable :: signal_arr
    double precision, dimension(:,:), allocatable :: s_guess, theta_guess
    double precision, dimension(:,:), allocatable :: s_grid, theta_grid, R_mag_data, Z_mag_data
!
  end module signal_array_mod
!  
