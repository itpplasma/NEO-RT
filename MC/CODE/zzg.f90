  function zzg()
    real :: zzg, xi
    call random_number(xi)
    zzg=xi
  end function zzg

  function gauss() 
    ! Gaussian random variable by the Box-Muller method
    real :: gauss, r, th, xi(2)
    call random_number(xi)
    r = sqrt(-2.0*log(xi(1)))
    th = 2.0*pi*xi(2)
    gauss = r*sin(th)
  end function gauss
