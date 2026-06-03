# NEO-RT C port — progress

Goal: literal C port of the golden-path physics, deps as packages (DVODE -> GSL
msadams / SUNDIALS CVODE; no NetCDF on golden path; LAPACK dptsv kept for
spline). Gate: OMP_NUM_THREADS=1 golden record at rtol=1e-8 (see baseline doc).

Each module is cross-checked against the real Fortran at rtol 1e-12 via a
Fortran reference dumper linked to libneo_rt + a C test (pattern: port_c/test).

## Verified (bit-faithful to Fortran)
- [x] spline (spline_coeff/spline_val_0) -- 163 checks. Calls LAPACK dptsv.
- [x] do_magfie_mod field core: Boozer reader (inp_swi 8/9), profile/mode
      splines, do_magfie evaluation -- 297 checks (33 theta x 9 quantities).
- [x] do_magfie_pert_mod amplitude (do_magfie_pert_amp) -- ported, bn cross-check pending.

## Remaining (golden path, port order)
- [ ] profiles (neort_profiles): plasma.in/profile.in reading, thermodynamic forces.
- [ ] driftorbit: resonance condition, velocity-space setup.
- [ ] orbit (neort_orbit): bounce integrals -- **DVODE substitution, the 1e-8 risk**.
      Validate GSL msadams vs DVODE on one bounce integral BEFORE porting the rest.
- [ ] freq (neort_freq): canonical frequencies + frequency splines.
- [ ] magfie (neort_magfie): check_magfie sampling -> _magfie.out (bn, eps_exp columns).
- [ ] resonance, transport, nonlin (nonlin off in golden, must still link).
- [ ] neort driver + output writers (the 5 .out files), main.
- [ ] CMake/Makefile producing neo_rt_c.x; run_gate.sh against it.

## Key facts
- Golden cases: inp_swi=9, pertfile=T, vsteps=512, s=0.1..0.9.
- DVODE call: Adams (method_flag=10), rtol=1e-9, atol=1e-10, events nevents=2.
- Shared spline j_start is value-independent (x clamped in range), so a single
  static index is faithful.
- Cross-check harness builds: gfortran ref against build/libneo_rt.a +
  build/libspline.a + build/libvode.a -llapack -lblas.
