from fffi import fortran_library, fortran_module

#libneo_rt = fortran_library('neo_rt', path='/home/calbert/build/NEO-RT/')
libneo_rt = fortran_library('neo_rt', path='/Users/ert/build/NEO-RT/')

magfie = fortran_module(libneo_rt, name='do_magfie_mod')
magfie.fdef("""
  integer :: m0b, n0b, nflux, nfp, nmode
  integer :: inp_swi
  double precision :: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds, q, dqds, iota, R0, a, eps, B0h, B00, bfac

  subroutine do_magfie_init
  end

  subroutine do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    double precision, dimension(:),       intent(in)         :: x
    double precision,                     intent(out)        :: bmod
    double precision,                     intent(out)        :: sqrtg
    double precision, dimension(:),       intent(out)        :: bder
    double precision, dimension(:),       intent(out)        :: hcovar
    double precision, dimension(:),       intent(out)        :: hctrvr
    double precision, dimension(:),       intent(out)        :: hcurl
  end
""")

magfie_pert = fortran_module(libneo_rt, name='do_magfie_pert_mod')
magfie_pert.fdef("""
  subroutine do_magfie_pert_init
  end

  subroutine do_magfie_pert_amp(x, bamp)
    real(8), dimension(:),       intent(in)         :: x
    complex(8),                  intent(out)        :: bamp
  end

  subroutine do_magfie_pert(x, bmod)
    real(8), dimension(:),       intent(in)         :: x
    complex(8),                  intent(out)        :: bmod
  end
""")

driftorb = fortran_module(libneo_rt, name='driftorbit')
#driftorb.fdef("""
#	
#""")

libneo_rt.compile(verbose=1)
magfie.load()
magfie_pert.load()
driftorb.load()


libneo_rt_mc = fortran_library('neo_rt_mc', path='/Users/ert/build/MC/')
libneo_rt_mc.fdef("""\
  subroutine magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    double precision, dimension(:),       intent(in)         :: x
    double precision,                     intent(out)        :: bmod
    double precision,                     intent(out)        :: sqrtg
    double precision, dimension(:),       intent(out)        :: bder
    double precision, dimension(:),       intent(out)        :: hcovar
    double precision, dimension(:),       intent(out)        :: hctrvr
    double precision, dimension(:),       intent(out)        :: hcurl
  end
""")
libneo_rt_mc.compile(verbose=1)
libneo_rt_mc.load()
