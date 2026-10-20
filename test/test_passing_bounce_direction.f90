program test_passing_bounce_direction
    use iso_fortran_env, only: dp=>real64
    use neort_orbit, only: bounce_integral, th0
    use driftorbit, only: etatp, sign_vpar
    use util, only: pi
    implicit none
    integer :: sigma,orientation
    real(dp), parameter :: angular_speed=1.7_dp
    real(dp) :: oracle_time,answer(3),initial(2)
    etatp=1.0_dp
    th0=0.3_dp
    initial=[th0,0.0_dp]
    oracle_time=2*pi/angular_speed
    do orientation=1,-1,-2
      do sigma=1,-1,-2
        sign_vpar=real(sigma,dp)
        answer=bounce_integral(1.0_dp,0.5_dp,2,initial,oracle_time/5,streaming)
        print *, 'sign,orientation,time,endpoint:',sigma,orientation,answer(1),answer(2)
        if(abs(answer(1)/oracle_time-1)>1e-10_dp) error stop 'one-turn time differs from analytic oracle'
        if(abs(answer(2)-(th0+sigma*orientation*2*pi))>1e-10_dp) error stop 'wrong signed one-turn endpoint'
      end do
    end do
contains
    subroutine streaming(v,eta,neq,t,y,ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: v,eta,t,y(neq)
        real(dp), intent(out) :: ydot(neq)
        ydot=0.0_dp
        ydot(1)=sign_vpar*orientation*angular_speed
    end subroutine streaming
end program
