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
- [x] collis_alp (loacol_nbi/coleff/onseff) -- 21 checks (with profiles).
- [x] profiles (neort_profiles): plasma.in/profile.in reading, thermodynamic forces -- 21 checks (vth,M_t,Om_tE,A1,A2,collision state).
- [x] driftorbit: state module (declarations only; logic lives in freq/resonance/transport).
- [x] magfie (neort_magfie): flux-surface average (eps,B0,Bmin,Bmax,etatp,etadt,th0).
- [x] orbit (neort_orbit): bounce/bounce_time via CVODE+root-finding -- 18 checks vs DVODE
      on real trapped+passing orbits, all <1e-8 (worst 8.07e-9, event-location sensitivity).
- [ ] freq (neort_freq): canonical frequencies + frequency splines.
- [ ] magfie sampling (check_magfie) -> _magfie.out bn, eps_exp columns.
- [ ] resonance, transport, nonlin (nonlin off in golden, must still link).
- [ ] neort driver + output writers (the 5 .out files), main.
- [ ] CMake/Makefile producing neo_rt_c.x; run_gate.sh against it.

## Integrator substitution — de-risk result (DATA)

Packages available: GSL 2.8 (gsl_odeiv2_step_msadams), SUNDIALS CVODE 7 (CV_ADAMS).
Experiment: same smooth ODE (damped pendulum), same IC, same tolerances as the
code (Adams, rtol=1e-9, atol=1e-10), DVODE vs GSL msadams:
  GSL msadams vs DVODE: y1 rel=5.68e-9, y2 rel=5.49e-9  -> passes 1e-8, <2x margin.

Interpretation: the golden record IS DVODE's 1e-9-converged trajectory. A
substitute integrator reproduces it only to ~DVODE's truncation level (~5e-9),
so the 1e-8 gate passes with thin margin on smooth integrals. RISK: orbit bounce
integrals use root-finding events (nevents=2); event location depends on the
integrator interpolant, which can exceed the ~5e-9 gap and break 1e-8.
Mitigation if it fails: port DVODE's Adams (method_flag=10) path to C for
bit-reproducible trajectories, instead of a package. Decide after porting orbit
and running the real bounce integral through both.

UPDATE (after reading orbit.f90): bounce_integral uses DVODE ROOT-FINDING
(g_fcn=bounceroots, nevents=2) to locate the orbit turning point. GSL gsl_odeiv2
has NO event location -> GSL cannot reproduce this faithfully. SUNDIALS CVODE
HAS root-finding (CVodeRootInit) and in CV_ADAMS mode is the true package match
for DVODE (Adams + events). So the integrator is CVODE, not GSL.

RESOLVED: CVODE (CV_ADAMS + SUNNonlinSol_FixedPoint) vs DVODE on the proxy ODE at
rtol=1e-9/atol=1e-10 agree to 9e-15 (machine precision) -- they share the VODE
Adams algorithm lineage. This is the decision: use SUNDIALS CVODE. Margin under
the 1e-8 gate is ~6 orders. Build flags: -lsundials_cvode -lsundials_core -lmpi
(SUNDIALS here is MPI-enabled, SUN_COMM_NULL needs -lmpi). CV_ADAMS REQUIRES an
explicit fixed-point nonlinear solver (default Newton needs a linear solver and
segfaults). Events: CVodeRootInit(cv, 2, rootfn) mirrors g_fcn=bounceroots.
Still to confirm on the REAL bounce integral once orbit is ported, but the
proxy result removes the integrator as a project risk.

Orbit RHS (timestep): nvar=7, ydot(1)=vpar*hctrvr(3), ydot(2)=-0.5 v^2 eta
hctrvr(3) hder(3) bmod, ydot(3)=Om_tB/v^2 (magnetic drift), ydot(4:)=0 here
(perturbed-Hamiltonian integrands added by a transport-side timestep variant).
Events: G1=sign_vpar_htheta*(theta-th0), G2=sign_vpar_htheta*(2pi-(theta-th0)).

## Key facts
- Golden cases: inp_swi=9, pertfile=T, vsteps=512, s=0.1..0.9.
- DVODE call: Adams (method_flag=10), rtol=1e-9, atol=1e-10, events nevents=2.
- Shared spline j_start is value-independent (x clamped in range), so a single
  static index is faithful.
- Cross-check harness builds: gfortran ref against build/libneo_rt.a +
  build/libspline.a + build/libvode.a -llapack -lblas.

## CRITICAL: achievable cross-language tolerance vs the 1e-8 gate (DATA)

Measured CVODE-vs-DVODE agreement ladder (faithful C port, CVODE is the closest
package integrator to DVODE):
  proxy smooth ODE        9e-15
  real bounce integral    ~8e-9
  in-range frequencies    1.48e-8   <-- already just OVER the 1e-8 gate
  near-tpb extrapolation  ~2e-2 values / ~1.5e-1 deta-derivatives (tiny eta sliver)

Conclusion: the rtol=1e-8 golden gate effectively demands bit-reproduction of
DVODE's exact trajectory. Two DIFFERENT correct Adams integrators diverge ~8e-9
in bounce integrals, amplified to ~1.5e-8 in frequencies and far more in the
near-boundary extrapolation (slope built from 2 boundary bounce points over tiny
log-spaced eta gaps). The final transport (D11/D12/torque) integrates these, so
the C/Rust/Julia ports will likely match the golden to ~1e-7, NOT 1e-8.

DECISION NEEDED (affects all three ports): either
 (a) define the cross-language gate at ~1e-6..1e-7 (accept that different correct
     integrators differ at ~1e-8), or
 (b) port DVODE's exact Adams (method_flag=10) path for bit-reproducibility
     (contradicts "use a package", large per-language effort).
The orbit/field/profile layers are bit-faithful; only the integrator trajectory
divergence drives this.

DECISION (user): relax the cross-language gate to rtol ~1e-7. Different correct
Adams integrators differ at ~1e-8; 1e-7 gives clean margin and keeps the package
approach. Port gate = compare port .out files to golden.h5 at rtol 1e-7 (a
separate comparator; the Fortran's own test_golden_record stays at 1e-8).
