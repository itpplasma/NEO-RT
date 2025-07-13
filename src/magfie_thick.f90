module magfie_thick_mod
    implicit none
    
    private
    public :: calculate_gyroradius, evaluate_field_with_gyroradius
    
contains

    subroutine calculate_gyroradius(v, B, qi, mi, rho_gyro)
        real(8), intent(in) :: v
        real(8), intent(in) :: B  
        real(8), intent(in) :: qi
        real(8), intent(in) :: mi
        real(8), intent(out) :: rho_gyro
        
        rho_gyro = mi * v / (qi * B)
    end subroutine calculate_gyroradius
    
    subroutine evaluate_field_with_gyroradius(x, rho_gyro, B_field, dB_dx)
        use do_magfie_mod, only: do_magfie
        
        real(8), intent(in) :: x(3)
        real(8), intent(in) :: rho_gyro
        real(8), intent(out) :: B_field(3)
        real(8), intent(out) :: dB_dx(3,3)
        real(8) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        
        B_field(1) = hctrvr(1) * bmod
        B_field(2) = hctrvr(2) * bmod  
        B_field(3) = hctrvr(3) * bmod
        
        dB_dx = 0.0d0
        dB_dx(1,1) = bder(1)
        dB_dx(2,2) = bder(2) 
        dB_dx(3,3) = bder(3)
    end subroutine evaluate_field_with_gyroradius
    
end module magfie_thick_mod