module native_perturbation
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: native_perturbation_value
    logical, save :: initialized=.false., enabled=.false.
    integer :: ntor
    integer, allocatable :: gm(:), modes(:)
    real(dp), allocatable :: rho(:), s(:), gc(:,:,:,:), gs(:,:,:,:)
    real(dp), allocatable :: seed(:,:)
    real(dp) :: zero_disks(6)
    complex(dp), allocatable :: amplitude(:,:)
    real(dp), save :: last_rho=-1, last_theta=0
    !$omp threadprivate(last_rho,last_theta)

contains

    subroutine initialize_native()
        integer :: u, ng, nmg, ns, nm, nseed
        character(8) :: magic
        inquire(file='bmod_native.dat',exist=enabled)
        if (.not. enabled) return
        open(newunit=u,file='bmod_native.dat',form='unformatted', &
             access='stream',status='old',action='read')
        read(u) magic
        if (magic /= 'POTNAT02') error stop 'native perturbation: invalid format'
        read(u) ng,nmg,ns,nm,ntor,nseed
        if (min(ng,ns)<2 .or. min(nmg,nm,nseed)<1 .or. ntor==0) &
            error stop 'native perturbation: invalid dimensions'
        allocate(rho(ng),gm(nmg),s(ns),modes(nm),seed(4,nseed))
        allocate(gc(4,ng-1,nmg,3),gs(4,ng-1,nmg,3),amplitude(ns,nm))
        read(u) rho,gm,gc,gs,s,modes,amplitude,seed
        read(u) zero_disks
        close(u)
        if (any(rho(2:)<=rho(:ng-1)) .or. any(s(2:)<=s(:ns-1))) &
            error stop 'native perturbation: nonmonotone radial axes'
        if (any(.not.ieee_is_finite(gc)) .or. any(.not.ieee_is_finite(gs))) &
            error stop 'native perturbation: nonfinite geometry'
        if (any(.not.ieee_is_finite(real(amplitude,dp))) .or. &
            any(.not.ieee_is_finite(aimag(amplitude)))) &
            error stop 'native perturbation: nonfinite amplitudes'
        print *, 'bmod_native.dat read (native radial nodes, no R-Z resampling)'
    end subroutine initialize_native

    integer function interval(axis,x) result(j)
        real(dp), intent(in) :: axis(:),x
        integer :: lo,hi,mid
        lo=1; hi=size(axis)
        do while (hi-lo>1)
            mid=(lo+hi)/2
            if (axis(mid)>x) then
                hi=mid
            else
                lo=mid
            end if
        end do
        j=lo
    end function interval

    subroutine geometry(r,t,v,jac)
        real(dp), intent(in) :: r,t
        real(dp), intent(out) :: v(3),jac(2,2)
        real(dp) :: a,dr,dt,x
        integer :: j,k,c
        j=interval(rho,r); x=r-rho(j)
        v=0; jac=0
        do c=1,3
            do k=1,size(gm)
                call native_harmonic(x,t,real(gm(k),dp), &
                    gc(1,j,k,c),gc(2,j,k,c),gc(3,j,k,c),gc(4,j,k,c), &
                    gs(1,j,k,c),gs(2,j,k,c),gs(3,j,k,c),gs(4,j,k,c),a,dr,dt)
                v(c)=v(c)+a
                if (c<=2) then
                    jac(c,1)=jac(c,1)+dr
                    jac(c,2)=jac(c,2)+dt
                end if
            end do
        end do
    end subroutine geometry

    subroutine initial_seed(r,z,rh,t)
        real(dp), intent(in) :: r,z
        real(dp), intent(out) :: rh,t
        integer :: k
        k=minloc((seed(1,:)-r)**2+(seed(2,:)-z)**2,dim=1)
        rh=seed(3,k); t=seed(4,k)
    end subroutine initial_seed

    subroutine invert(r,z,rh,t,v,inside)
        real(dp), intent(in) :: r,z
        real(dp), intent(out) :: rh,t,v(3)
        logical, intent(out) :: inside
        real(dp) :: jac(2,2),delta(2),residual
        integer :: ipiv(2),info,iteration,attempt
        inside=.false.
        rh=last_rho; t=last_theta
        if (rh<rho(1)) call initial_seed(r,z,rh,t)
        do attempt=1,2
            do iteration=1,40
                call geometry(rh,t,v,jac)
                delta=[r-v(1),z-v(2)]
                residual=maxval(abs(delta))
                if (residual<2.e-12_dp) then
                    inside=.true.
                    last_rho=rh; last_theta=t
                    return
                end if
                call dgesv(2,1,jac,2,ipiv,delta,2,info)
                if (info/=0) exit
                ! A converged tangential coordinate with an outward Newton
                ! step at a geometry boundary is outside source support.
                if (abs(delta(2))<1.e-10_dp) then
                    if (rh<=rho(1)+1.e-13_dp .and. delta(1)<0) return
                    if (rh>=rho(size(rho))-1.e-13_dp .and. delta(1)>0) return
                end if
                rh=max(rho(1),min(rho(size(rho)),rh+ &
                    max(-.1_dp,min(.1_dp,delta(1)))))
                t=modulo(t+max(-.3_dp,min(.3_dp,delta(2))),8*atan(1._dp))
            end do
            call initial_seed(r,z,rh,t)
        end do
        write(error_unit,'(a,3es24.16)') &
            'native inverse failed R,Z,residual: ',r,z,residual
        error stop 'native perturbation: inverse map did not converge'
    end subroutine invert

    subroutine native_perturbation_value(r,z,value,handled)
        real(dp), intent(in) :: r,z
        complex(dp), intent(out) :: value
        logical, intent(out) :: handled
        real(dp) :: rh,t,v(3),sq,w
        complex(dp) :: a
        logical :: inside,ready
        integer :: j,k
        !$omp atomic read acquire
        ready=initialized
        if (.not.ready) then
            !$omp critical(native_perturbation_init)
            if (.not.initialized) then
                call initialize_native()
                !$omp atomic write release
                initialized=.true.
            end if
            !$omp end critical(native_perturbation_init)
        end if
        handled=enabled
        value=(0._dp,0._dp)
        if (.not.enabled) return
        if (.not.ieee_is_finite(r) .or. .not.ieee_is_finite(z)) &
            error stop 'native perturbation: nonfinite cylindrical position'
        ! Certified subsets of zero support avoid a singular inverse in the
        ! uncharted core. These disks come from the native Fourier contours,
        ! with an angular-derivative bound, not an assumed circular plasma.
        if (sum(([r,z]-zero_disks(1:2))**2)<zero_disks(3)) return
        if (sum(([r,z]-zero_disks(4:5))**2)>zero_disks(6)) return
        call invert(r,z,rh,t,v,inside)
        if (.not.inside) return
        sq=rh*rh
        if (sq<s(1) .or. sq>s(size(s))) return
        j=interval(s,sq); w=(sq-s(j))/(s(j+1)-s(j))
        do k=1,size(modes)
            a=(1-w)*amplitude(j,k)+w*amplitude(j+1,k)
            value=value+a*exp(cmplx(0._dp,modes(k)*t,dp))
        end do
        value=value*exp(cmplx(0._dp,-ntor*v(3),dp))
    end subroutine native_perturbation_value

    include 'generated/native_harmonic.inc'
end module native_perturbation
