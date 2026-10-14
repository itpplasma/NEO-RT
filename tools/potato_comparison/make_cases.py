"""Materialize paired native inputs from one shared field/profile bundle."""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

AMU_G = 1.660538e-24
P_MASS_G = 1.6726e-24
EV_ERG = 1.602176e-12  # NEO-RT util.eV; differs from POTATO by < 4e-7
MASS = 2.014 * P_MASS_G / AMU_G
MODES = ("thin", "dph", "radial", "both")


def temperature(case, s):
    return case["T0_eV"] * (1 - 0.5 * (np.asarray(s) - case["s"]))


def write_profiles(case, destination):
    destination.mkdir(parents=True, exist_ok=True)
    grid = np.linspace(0.005, 0.98, 196)
    ti = temperature(case, grid)
    vt = np.sqrt(2 * ti * EV_ERG / (MASS * AMU_G))
    plasma = np.column_stack((grid, np.full_like(grid, 1e10),
                             np.zeros_like(grid), ti, np.ones_like(grid), ti))
    with (destination / "plasma.in").open("w") as f:
        f.write(f"% N am1 am2 Z1 Z2\n{len(grid)} {MASS:.16e} {MASS:.16e} 1 1\n")
        f.write("% s ni1 ni2 Ti1 Ti2 Te\n")
        np.savetxt(f, plasma, fmt="%.16e")
    np.savetxt(destination / "profile.in",
               np.column_stack((grid, np.zeros_like(grid), vt)), fmt="%.16e")
    density = [0.] * 9 + [1e10]
    poly = [0.] * 8 + [-case["T0_eV"] / 2,
                            case["T0_eV"] * (1 + case["s"] / 2)]
    with (destination / "profile_poly.in").open("w") as f:
        f.write("% Physical n [cm^-3], T [eV], Phi = 0; no energy scaling\n")
        f.write("% ni, Te, Ti, Phi; decreasing powers of s_pol\n")
        np.savetxt(f, [density, poly, poly, [0.] * 10], fmt="%.16e")


def materialize(root, case, destination):
    """Both solvers receive the same complete native bundle, with local aliases."""
    destination.mkdir(parents=True, exist_ok=True)
    for source in (root / "inputs/common", root / "inputs" / case["id"]):
        for path in source.iterdir():
            if path.is_file():
                shutil.copy2(path, destination / path.name)
    for name, target in (("in_file", "eq.bc"), ("in_file_pert", "pert.bc")):
        path = destination / name
        if path.is_symlink():
            path.unlink()
        path.symlink_to(target)


def potato_config(case):
    return f"""&potato_nml
    itest_type = 3
    edge_extension = .false.
    E_alpha = {case['T0_eV']:.16e}
    A_alpha = 2.014d0
    Z_alpha = 1d0
    rho_pol = 0.8d0
    rho_pol_max = 0.99d0
    scalfac_energy = 1d0
    scalfac_efield = 1d0
    Rmax_orbit = 200d0
    ntimstep = 30
    npoicut = 10000
    m_min = -3
    m_max = 3
    n_tor = -3
    nenerg = {case['nenerg']}
    thermen_max = 10d0
    nbox = 100
    unif_rho_pol = .true.
    adaptive_jperp = .true.
    npoi_init = 101
    nlagr_sampling = 3
    eps_sampling = 3d-3
    itermax_sampling = 5
    profile_file = 'profile_poly.in'
/
"""


def neort_config(case, rho, mode):
    vth = np.sqrt(2 * temperature(case, rho**2) * EV_ERG / (MASS * AMU_G))
    return f"""&params
    s = {rho**2:.16e}
    M_t = 0d0
    qs = 1d0
    ms = {MASS:.16e}
    vth = {vth:.16e}
    epsmn = 1d0
    pertfile_scale = 1d0
    m0 = 0
    mph = 3
    supban = .false.
    fow_radial = {'.true.' if mode in ('radial', 'both') else '.false.'}
    drift_phase = {'.true.' if mode in ('dph', 'both') else '.false.'}
    magdrift = .true.
    nopassing = .false.
    noshear = .false.
    pertfile = .true.
    nonlin = .false.
    comptorque = .true.
    bfac = 1d0
    efac = 1d0
    inp_swi = 9
    vsteps = 512
    log_level = 1
/
"""


def main():
    root, solver, cid, dest = sys.argv[1:5]
    root, dest = Path(root).resolve(), Path(dest).resolve()
    case = next(c for c in json.loads((root / "inputs/scan.json").read_text())
                ["cases"] if c["id"] == cid)
    materialize(root, case, dest)
    if solver == "potato":
        (dest / "potato.in").write_text(potato_config(case))
    else:
        rho, mode = float(sys.argv[5]), sys.argv[6]
        (dest / "driftorbit.in").write_text(neort_config(case, rho, mode))


if __name__ == "__main__":
    main()
