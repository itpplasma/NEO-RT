from os.path import expanduser
from fffi import FortranLibrary, FortranModule

libneo_rt = FortranLibrary('neo_rt', path=expanduser('~/build/NEO-RT/'))

magfie = FortranModule(libneo_rt, name='do_magfie_mod')
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

magfie_pert = FortranModule(libneo_rt, name='do_magfie_pert_mod')
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

driftorb = FortranModule(libneo_rt, name='driftorbit')
# driftorb.fdef("""
#
# """)

libneo_rt.compile(verbose=1)
magfie.load()
magfie_pert.load()
driftorb.load()


libneo_rt_mc = FortranLibrary('neo_rt_mc', path=expanduser('~/build/MC/'))
libneo_rt_mc.fdef("""\
  subroutine magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    double precision, dimension(3),       intent(in)         :: x
    double precision,                     intent(out)        :: bmod
    double precision,                     intent(out)        :: sqrtg
    double precision, dimension(3),       intent(out)        :: bder
    double precision, dimension(3),       intent(out)        :: hcovar
    double precision, dimension(3),       intent(out)        :: hctrvr
    double precision, dimension(3),       intent(out)        :: hcurl
  end

  subroutine cyl_coord(x, x_cyl)
    double precision, dimension(3),       intent(in)         :: x
    double precision, dimension(3),       intent(out)        :: x_cyl
  end

  subroutine velo(tau, z, vz)
    double precision              ,       intent(in)         :: tau
    double precision, dimension(3),       intent(in)         :: z
    double precision, dimension(3),       intent(out)        :: vz
  end
""")


parmot_mod = FortranModule(libneo_rt_mc, name='parmot_mod')
parmot_mod.fdef("""\
  double precision :: rmu,ro0
""")

libneo_rt_mc.compile(verbose=1)
libneo_rt_mc.load()
parmot_mod.load()

print(libneo_rt_mc.csource)
