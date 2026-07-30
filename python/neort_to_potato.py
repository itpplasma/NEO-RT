#!/usr/bin/env python3
"""
Convert NEO-RT profile inputs to POTATO profile_poly.in.

Reads plasma.in and profile.in from a NEO-RT run and writes a polynomial
profile file compatible with POTATO's profile_poly.in format.

Density and temperature come from plasma.in.  The electrostatic potential
phi_e is derived from the ExB Mach number M_t in profile.in using:

    Er(s_pol) = sign_theta * dpsi_pol * vth(s) * M_t(s) / (c * R0)

    phi_e(s_pol) = -integral_0^s_pol Er ds_pol'   [statV, zero on axis]

RADIAL COORDINATE (this bit is easy to get wrong):
  NEO-RT's plasma.in / profile.in are tabulated against the normalized
  TOROIDAL flux s_tor.  POTATO's denstemp_of_psi / phielec_of_psi evaluate
  these polynomials at the normalized POLOIDAL flux
  s_pol = (psi - psi_axis)/(psi_sep - psi_axis).  Every profile must
  therefore be re-tabulated s_tor -> s_pol (libneo FluxConverter, from the
  EQDSK q-profile) BEFORE fitting.  Skipping this misplaces every profile
  radially -- it moved the Er=0 surface from rho_pol~0.90 to ~0.72.

FLUX GAUGE (likewise):
  dpsi_pol must be POTATO's own poloidal gauge,
  (PsiedgeVs - PsiaxisVs) * 1e8 Mx from the EQDSK header, because
  phielec_of_psi divides dPhi/ds_pol by exactly that.  NEO-RT's `psi_pr`
  in *_magfie_param.out is the TOROIDAL derivative dpsi_tor/ds (-3.59e7 Mx
  for AUG #30835) and using it here inflates Omega_ExB by x1.945 at every
  radius.

R0 is read from any *_magfie_param.out in the run directory (or given with
--R0).  sign_theta defaults to -1 (standard AUG convention).

Usage:
  python neort_to_potato.py plasma.in profile.in [options]
  python neort_to_potato.py --plot-only [--output profile_poly.in]

Examples:
  # Auto-read psi_pr and R0 from magfie_param files in same directory:
  python neort_to_potato.py plasma.in profile.in

  # Provide geometry parameters explicitly:
  python neort_to_potato.py plasma.in profile.in \\
      --psi-pr -35893832.0 --R0 171.47

  # Write profile and save diagnostic plots:
  python neort_to_potato.py plasma.in profile.in --plot

  # Zero phi_e (no electric field):
  python neort_to_potato.py plasma.in profile.in --phi-zero

  # Plot an existing profile_poly.in without input files:
  python neort_to_potato.py --plot-only
  python neort_to_potato.py --plot-only --output path/to/profile_poly.in
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

from libneo import FluxConverter, eqdsk_base

C_CGS = 2.9979e10    # cm/s


def _fortran_float(s):
    return float(s.replace('D', 'E').replace('d', 'e'))
EV_TO_ERG = 1.6e-12  # erg/eV
MU_CGS = 1.6726e-24  # g
POLY_ORDER = 9


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def read_plasma_in(path):
    """Return s [1], ni_1 [cm^-3], Ti_1 [eV], Te [eV], am1 [amu], Z1."""
    with open(path) as f:
        f.readline()
        parts = f.readline().split()
        n_pts = int(parts[0])
        am1 = _fortran_float(parts[1])
        Z1 = _fortran_float(parts[3])
        f.readline()
        rows = []
        for _ in range(n_pts):
            rows.append([_fortran_float(v) for v in f.readline().split()])
    data = np.array(rows)
    return data[:, 0], data[:, 1], data[:, 3], data[:, 5], am1, Z1


def read_profile_in(path):
    """Return s [1] and M_t (ExB Mach number) from first two columns."""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    rows.append([_fortran_float(parts[0]), _fortran_float(parts[1])])
                except ValueError:
                    pass
    data = np.array(rows)
    return data[:, 0], data[:, 1]


def read_profile_poly_in(path):
    """Return (poly_n, poly_Te, poly_Ti, poly_phi) from an existing profile_poly.in."""
    polys = []
    with open(path) as f:
        for line in f:
            if line.startswith(" #") or line.startswith("#"):
                continue
            coeffs = [float(v) for v in line.split()]
            if coeffs:
                polys.append(np.array(coeffs))
    if len(polys) < 4:
        raise ValueError(f"Expected 4 polynomial rows in {path}, got {len(polys)}")
    return polys[0], polys[1], polys[2], polys[3]


def read_magfie_param(path):
    """Parse a *_magfie_param.out file and return dict of named values."""
    values = {}
    pattern = re.compile(r'check_magfie:\s+(\S+)\s*=\s*([+-]?\d[\d.eE+\-]*)')
    with open(path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                values[m.group(1)] = float(m.group(2))
    return values


def find_magfie_params(directory):
    """Return parameter dict from the first *_magfie_param.out in directory."""
    files = glob.glob(os.path.join(directory, '*_magfie_param.out'))
    if not files:
        return {}
    return read_magfie_param(files[0])


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------

def read_eqdsk_geometry(eqdsk_path):
    """Return (stor2spol callable, dpsi_pol [Mx]) from an EQDSK g-file.

    dpsi_pol = psi_sep - psi_axis in POTATO's field_eq gauge; phielec_of_psi
    divides dPhi/ds_pol by exactly this, so Er must be built with it.
    """
    eq = eqdsk_base.read_eqdsk(eqdsk_path)
    converter = FluxConverter(eq["qprof"], eq["s_pol"])
    dpsi_pol = (eq["PsiedgeVs"] - eq["PsiaxisVs"]) * 1.0e8   # Wb -> Mx
    return (lambda s_tor: np.asarray(converter.stor2spol(s_tor))), dpsi_pol


def compute_vth(Ti_eV, am1):
    """Ion thermal velocity [cm/s]."""
    return np.sqrt(2.0 * Ti_eV * EV_TO_ERG / (am1 * MU_CGS))


def compute_er(s_plasma, Ti_eV, s_profile, M_t, am1, dpsi_pol, R0, sign_theta):
    """Er [statV per unit s_pol] from the ExB Mach number.

    Uses the POLOIDAL flux difference (see module docstring): with
    phi_e = -int Er ds_pol, POTATO recovers
    Omega_E = c dPhi/dpsi_pol = -sign_theta * vth * M_t/R0, matching NEO-RT's
    Om_tE = vth*M_t/R0 for the standard sign_theta = -1.

    Note M_t and vth are functions of the flux SURFACE; they are evaluated on
    the plasma.in surfaces (labelled by s_tor) and only the integration
    variable is s_pol.
    """
    vth = compute_vth(Ti_eV, am1)
    M_t_i = np.interp(s_plasma, s_profile, M_t)
    return sign_theta * dpsi_pol * vth * M_t_i / (C_CGS * R0)


def integrate_er_to_phi(s, er):
    """phi_e(s) [statV] via trapezoidal integration from axis outward.

    Matches import_neo2_profiles.f90 convention: phi_e = -integral_0^s Er ds',
    so phi_e(0)=0 and d(phi_e)/ds = -Er.
    """
    idx = np.argsort(s)
    s_d = s[idx]
    er_d = er[idx]
    phi_d = np.zeros(len(s))
    for i in range(1, len(s)):
        ds = s_d[i] - s_d[i - 1]
        phi_d[i] = phi_d[i - 1] - 0.5 * (er_d[i - 1] + er_d[i]) * ds
    phi = np.empty(len(s))
    phi[idx] = phi_d
    return phi


# ---------------------------------------------------------------------------
# Polynomial helpers
# ---------------------------------------------------------------------------

def polyfit_potato(s, y, order):
    """Polynomial fit; coefficients in descending power order (np.polyval)."""
    return np.polyfit(s, y, order)


def fmt_row(coeffs):
    return "  " + "  ".join(f"{c: .15e}" for c in coeffs)


def write_profile_poly(path, poly_n, poly_Te, poly_Ti, poly_phi, order, am1, Z1):
    """Write POTATO profile_poly.in."""
    with open(path, "w") as f:
        f.write(" # Generated from NEO-RT profiles by neort_to_potato.py\n")
        f.write(f" # Polynomial order: {order}  am1: {am1}  Z1: {Z1}\n")
        f.write(fmt_row(poly_n) + "\n")
        f.write(fmt_row(poly_Te) + "\n")
        f.write(fmt_row(poly_Ti) + "\n")
        f.write(fmt_row(poly_phi) + "\n")


# ---------------------------------------------------------------------------
# Plotting  (mirrors import_neo2_profiles.f90 diagnostic figures)
# ---------------------------------------------------------------------------

def _savefig(plt, fig, outdir, name):
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


def save_plots(outdir, s_plasma, ni_1, Ti_1, Te, s_profile, M_t,
               er, phi_e, poly_n, poly_Ti, poly_phi, order):
    """Save density_fit.png, temperature_fit.png, Er_fit.png, phi_e_fit.png."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available – skipping plots.", file=sys.stderr)
        return

    s_fine = np.linspace(0.0, 1.0, 300)

    fig, ax = plt.subplots()
    ax.plot(s_plasma, ni_1, "ro", ms=3, label="plasma.in data")
    ax.plot(s_fine, np.polyval(poly_n, s_fine), "b-", label="polynomial fit")
    ax.set_xlabel("Normalised poloidal flux")
    ax.set_ylabel("Ion density [cm⁻³]")
    ax.set_title("Ion density profile")
    ax.legend()
    _savefig(plt, fig, outdir, "density_fit.png")

    fig, ax = plt.subplots()
    ax.plot(s_plasma, Ti_1, "ro", ms=3, label="plasma.in data")
    ax.plot(s_fine, np.polyval(poly_Ti, s_fine), "b-", label="polynomial fit")
    ax.set_xlabel("Normalised poloidal flux")
    ax.set_ylabel("Ion temperature [eV]")
    ax.set_title("Ion temperature profile")
    ax.legend()
    _savefig(plt, fig, outdir, "temperature_fit.png")

    er_fine = np.polyval(polyfit_potato(s_plasma, er, order + 1), s_fine)
    fig, ax = plt.subplots()
    ax.plot(s_plasma, er, "ro", ms=3, label="from M_t")
    ax.plot(s_fine, er_fine, "b-", label="polynomial fit")
    ax.set_xlabel("Normalised poloidal flux")
    ax.set_ylabel("Er [statV per unit s]")
    ax.set_title("Radial electric field")
    ax.legend()
    _savefig(plt, fig, outdir, "Er_fit.png")

    fig, ax = plt.subplots()
    ax.plot(s_plasma, phi_e, "ro", ms=3, label="integrated Er")
    ax.plot(s_fine, np.polyval(poly_phi, s_fine), "b-", label="polynomial fit")
    ax.set_xlabel("Normalised poloidal flux")
    ax.set_ylabel("φₑ [statV]")
    ax.set_title("Electrostatic potential")
    ax.legend()
    _savefig(plt, fig, outdir, "phi_e_fit.png")


def save_plots_from_poly(outdir, poly_n, poly_Te, poly_Ti, poly_phi):
    """Plot polynomial profiles read from an existing profile_poly.in."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available – skipping plots.", file=sys.stderr)
        return

    s_fine = np.linspace(0.0, 1.0, 300)

    fig, ax = plt.subplots()
    ax.plot(s_fine, np.polyval(poly_n, s_fine), "b-")
    ax.set_xlabel("Normalised poloidal flux")
    ax.set_ylabel("Ion density [cm⁻³]")
    ax.set_title("Ion density profile")
    _savefig(plt, fig, outdir, "density_fit.png")

    fig, ax = plt.subplots()
    ax.plot(s_fine, np.polyval(poly_Ti, s_fine), "b-", label="Ti")
    ax.plot(s_fine, np.polyval(poly_Te, s_fine), "r-", label="Te")
    ax.set_xlabel("Normalised poloidal flux")
    ax.set_ylabel("Temperature [eV]")
    ax.set_title("Temperature profiles")
    ax.legend()
    _savefig(plt, fig, outdir, "temperature_fit.png")

    poly_er = -np.polyder(poly_phi)
    fig, ax = plt.subplots()
    ax.plot(s_fine, np.polyval(poly_er, s_fine), "b-")
    ax.set_xlabel("Normalised poloidal flux")
    ax.set_ylabel("Er [statV per unit s]")
    ax.set_title("Radial electric field")
    _savefig(plt, fig, outdir, "Er_fit.png")

    fig, ax = plt.subplots()
    ax.plot(s_fine, np.polyval(poly_phi, s_fine), "b-")
    ax.set_xlabel("Normalised poloidal flux")
    ax.set_ylabel("φₑ [statV]")
    ax.set_title("Electrostatic potential")
    _savefig(plt, fig, outdir, "phi_e_fit.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert NEO-RT plasma.in/profile.in to POTATO profile_poly.in",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("plasma_file", nargs="?", help="NEO-RT plasma.in")
    parser.add_argument("profile_file", nargs="?", help="NEO-RT profile.in")
    parser.add_argument("--output", default="profile_poly.in",
                        help="Output file (default: profile_poly.in)")
    parser.add_argument("--order", type=int, default=POLY_ORDER,
                        help="Polynomial order (default: 9)")
    parser.add_argument("--eqdsk", default="efit.eqdsk", metavar="FILE",
                        help="EQDSK g-file: supplies the s_tor->s_pol map and "
                             "the poloidal flux gauge (default: efit.eqdsk)")
    parser.add_argument("--dpsi-pol", type=float, default=None, metavar="MX",
                        help="Override psi_sep-psi_axis [Mx = G cm²]; normally "
                             "taken from the EQDSK header")
    parser.add_argument("--R0", type=float, default=None, metavar="CM",
                        help="Major radius [cm]")
    parser.add_argument("--sign-theta", type=float, default=-1.0,
                        help="Sign of poloidal field (default: -1)")
    parser.add_argument("--phi-zero", action="store_true",
                        help="Set phi_e = 0")
    parser.add_argument("--plot", action="store_true",
                        help="Save diagnostic plots alongside the output file")
    parser.add_argument("--plot-only", action="store_true",
                        help="Save diagnostic plots without writing profile_poly.in; "
                             "reads existing --output file if no input files given")
    args = parser.parse_args()

    have_inputs = args.plasma_file is not None and args.profile_file is not None

    if not have_inputs:
        if not args.plot_only:
            parser.error("plasma_file and profile_file are required unless --plot-only is used")
        poly_n, poly_Te, poly_Ti, poly_phi = read_profile_poly_in(args.output)
        outdir = os.path.dirname(os.path.abspath(args.output))
        save_plots_from_poly(outdir, poly_n, poly_Te, poly_Ti, poly_phi)
        return

    s_plasma, ni_1, Ti_1, Te, am1, Z1 = read_plasma_in(args.plasma_file)
    s_profile, M_t = read_profile_in(args.profile_file)

    # plasma.in/profile.in are tabulated in s_tor; POTATO evaluates every
    # polynomial at s_pol.  Re-tabulate BEFORE fitting (see module docstring).
    if not os.path.exists(args.eqdsk):
        parser.error(
            f"EQDSK '{args.eqdsk}' not found. It is required for the "
            "s_tor->s_pol map and the poloidal flux gauge; pass --eqdsk."
        )
    stor2spol, dpsi_pol_eqdsk = read_eqdsk_geometry(args.eqdsk)
    spol_plasma = stor2spol(s_plasma)
    print(f"s_tor -> s_pol via {args.eqdsk}: "
          f"s_tor[{s_plasma.min():.4f},{s_plasma.max():.4f}] -> "
          f"s_pol[{spol_plasma.min():.4f},{spol_plasma.max():.4f}]")

    poly_n  = polyfit_potato(spol_plasma, ni_1, args.order)
    poly_Ti = polyfit_potato(spol_plasma, Ti_1, args.order)
    poly_Te = polyfit_potato(spol_plasma, Te,   args.order)

    er    = None
    phi_e = None

    if args.phi_zero:
        poly_phi = np.zeros(args.order + 1)
        print("phi_e set to zero (--phi-zero)")
    else:
        # POLOIDAL gauge from the EQDSK -- NOT NEO-RT's toroidal psi_pr
        dpsi_pol = args.dpsi_pol if args.dpsi_pol is not None else dpsi_pol_eqdsk
        R0 = args.R0
        if R0 is None:
            run_dir = os.path.dirname(os.path.abspath(args.plasma_file))
            R0 = find_magfie_params(run_dir).get("R0")

        if R0 is None:
            print(
                "Warning: R0 not found – setting phi_e = 0.\n"
                "  Provide --R0, or place a *_magfie_param.out file next to "
                "plasma.in.",
                file=sys.stderr,
            )
            poly_phi = np.zeros(args.order + 1)
        else:
            er    = compute_er(s_plasma, Ti_1, s_profile, M_t, am1,
                               dpsi_pol, R0, args.sign_theta)
            # integrate over s_pol, the variable POTATO differentiates in
            phi_e = integrate_er_to_phi(spol_plasma, er)
            poly_phi = polyfit_potato(spol_plasma, phi_e, args.order)
            print(
                f"phi_e from M_t  "
                f"(dpsi_pol={dpsi_pol:.4e} Mx [EQDSK poloidal gauge], "
                f"R0={R0:.2f} cm, sign_theta={args.sign_theta:+.0f})"
            )

    if not args.plot_only:
        write_profile_poly(
            args.output, poly_n, poly_Te, poly_Ti, poly_phi,
            args.order, am1, Z1,
        )
        print(f"Written {args.output}")

    if args.plot or args.plot_only:
        if er is None:
            er    = np.zeros(len(s_plasma))
            phi_e = np.zeros(len(s_plasma))
        outdir = os.path.dirname(os.path.abspath(args.output))
        save_plots(outdir, s_plasma, ni_1, Ti_1, Te, s_profile, M_t,
                   er, phi_e, poly_n, poly_Ti, poly_phi, args.order)


if __name__ == "__main__":
    main()
