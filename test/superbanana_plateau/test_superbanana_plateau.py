#!/usr/bin/env python3
"""
Shaing superbanana-plateau (SBP) regression test for NEO-RT.

Runs NEO-RT with ``supban = .true.`` on a circular tokamak and compares the
numeric SBP diffusion coefficient D11 against an independent analytic
superbanana-plateau prediction built from the g(kappa) drift resonance of
Shaing 2009 PPCF 51 035009 Eq. (8) and the "Superbanana plateau" subsection of
doc/driftorbit.lyx.

The analytic model shares the field geometry with the code (parsed from the
generated ``*_magfie_param.out``) and evaluates the same resonance-line D11
integral, so a match confirms the reinstated SBP branch reproduces the
large-aspect-ratio analytic flux (historical marker commit ec133c3,
"correct superbanana plateau fluxes").
"""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

TEST_DIR = Path(__file__).resolve().parent
CASE_NAME = "driftorbit_test"

# Perturbation / species parameters fixed by driftorbit_test.in.
VTH = 1.0e8
M_T = 0.001
MPH = 18
M0 = 0
EPSMN = 1.0e-3
S = 3.9063e-2

# cgs constants (match src/util.f90).
C = 2.997925e10
QE = 4.803204e-10
AMU = 1.660538e-24
QI = 1.0 * QE
MI = 2.014 * AMU
SIGN_THETA = -1.0  # left-handed Boozer angles (src/do_magfie_standalone.f90)

# Analytic-model comparison tolerance.  The residual is the O(eps) difference
# between the large-aspect-ratio analytic bounce quantities and the numeric
# orbit integration (~3% at this flux surface, eps ~ 0.05).
D11_TOLERANCE = 0.06


def find_executable() -> Path:
    candidates = [
        TEST_DIR.parent.parent / "build" / "neo_rt.x",
        Path(os.environ.get("NEO_RT_EXE", "")),
    ]
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    result = shutil.which("neo_rt.x")
    if result:
        return Path(result)
    pytest.skip("neo_rt.x executable not found")


@pytest.fixture(scope="module")
def executable() -> Path:
    return find_executable()


def setup_work_dir(work_dir: Path) -> None:
    input_file = TEST_DIR / f"{CASE_NAME}.in"
    if not input_file.exists():
        pytest.fail(f"Required input file missing: {input_file}")
    (work_dir / f"{CASE_NAME}.in").symlink_to(input_file)

    in_file_src = TEST_DIR / "in_file"
    if not in_file_src.exists():
        pytest.fail(f"Required input file missing: {in_file_src}")
    (work_dir / "in_file").symlink_to(in_file_src.resolve())


def run_neort(executable: Path, work_dir: Path) -> None:
    result = subprocess.run(
        [str(executable), CASE_NAME],
        cwd=work_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"neo_rt.x failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )


def parse_field_params(work_dir: Path) -> dict:
    """Read flux-surface geometry from the generated magfie_param output."""
    param_file = work_dir / f"{CASE_NAME}_magfie_param.out"
    wanted = ("eps", "R0", "a", "psi_pr", "B0", "q", "iota", "dVds")
    params: dict = {}
    with open(param_file, "r") as f:
        for line in f:
            fields = line.split()
            # Lines read: "check_magfie: <name> = <value>"
            if len(fields) >= 4 and fields[0] == "check_magfie:" and fields[2] == "=":
                name = fields[1]
                if name in wanted:
                    params[name] = float(fields[3])
    missing = [k for k in wanted if k not in params]
    if missing:
        pytest.fail(f"Missing field params {missing} in {param_file}")
    return params


# --- analytic superbanana-plateau model ------------------------------------


def comelp(k2: float):
    """Complete elliptic integrals K(k), E(k) with k^2 = k2 (matches Fortran)."""
    pk = 1.0 - k2
    ak = (((0.01451196212 * pk + 0.03742563713) * pk + 0.03590092383) * pk
          + 0.09666344259) * pk + 1.38629436112
    bk = (((0.00441787012 * pk + 0.03328355346) * pk + 0.06880248576) * pk
          + 0.12498593597) * pk + 0.5
    K = ak - bk * np.log(pk)
    ae = (((0.01736506451 * pk + 0.04757383546) * pk + 0.0626060122) * pk
          + 0.44325141463) * pk + 1.0
    be = (((0.00526449639 * pk + 0.04069697526) * pk + 0.09200180037) * pk
          + 0.2499836831) * pk
    E = ae - be * np.log(pk)
    return K, E


def d11_analytic(p: dict, neta: int = 3000, npsi: int = 6001) -> float:
    """Analytic SBP D11 (normalized by the ripple-plateau Dp), same field."""
    eps, R0, a = p["eps"], p["R0"], p["a"]
    B0, psi_pr, q = p["B0"], p["psi_pr"], p["q"]
    iota, dVds = p["iota"], p["dVds"]

    om_te = VTH * M_T / R0
    chi = SIGN_THETA * psi_pr
    depsdr = eps / (a * np.sqrt(S))     # = 1/R0 for a circular flux surface
    hth = iota / R0                     # contravariant poloidal, leading order
    meff = M0 + q * MPH

    Dp = np.pi * VTH**3 / (16.0 * R0 * iota * (QI * B0 / (MI * C))**2)
    dsdreff = 2.0 / a * np.sqrt(S)
    vmin, vmax = 0.01 * VTH, 5.0 * VTH

    etatp = 1.0 / (B0 * (1.0 + eps))
    etadt = 1.0 / (B0 * (1.0 - eps))
    etas = np.linspace(etatp * (1.0 + 1e-8), etadt * (1.0 - 1e-8), neta)

    psi = np.linspace(-np.pi / 2, np.pi / 2, npsi)
    spsi = np.sin(psi)

    integrand = np.zeros(neta)
    for i, eta in enumerate(etas):
        k2 = (1.0 - eta * B0 * (1.0 - eps)) / (2.0 * eta * B0 * eps)
        if k2 <= 1e-12 or k2 >= 1.0 - 1e-12:
            continue
        K, E = comelp(k2)
        drift_shape = 2.0 * E / K - 1.0
        if drift_shape <= 0.0:  # no ell=0 resonance for these pitch values
            continue

        # Resonance Om_tE + <Om_tB> = 0 with <Om_tB> = -C0*drift_shape*v^2
        # gives v^2 explicitly (removes the velocity-space singularity).
        c0 = (C * 0.5 * MI * eta * B0 / (QI * chi)) * depsdr
        v2 = om_te / (c0 * drift_shape)
        if v2 <= 0.0:
            continue
        v = np.sqrt(v2)
        if v < vmin or v > vmax:
            continue
        ux = v / VTH

        # Bounce time and bounce-averaged perturbation for B = B0(1 - eps cos th),
        # h^theta = iota/R0 (large-aspect-ratio circular model field).
        kap = np.sqrt(k2)
        th = 2.0 * np.arcsin(kap * spsi)
        weight = 1.0 / np.sqrt(1.0 - k2 * spsi**2)
        Bth = B0 * (1.0 - eps * np.cos(th))
        Hb = EPSMN * np.trapezoid((2.0 - eta * Bth) * np.cos(meff * th) * weight,
                                  psi) / (2.0 * K)
        taub = 8.0 * K / (v * np.sqrt(2.0 * eta * B0 * eps) * hth)

        Hmn2 = Hb**2 * (MI * v2 / 2.0)**2
        d11int = (np.pi**1.5 * MPH**2 * C**2 * q * VTH
                  / (QI**2 * dVds * abs(psi_pr))
                  * ux**3 * np.exp(-ux**2) * taub * Hmn2)
        # Resonance-line Jacobian: |d Omph / d ux| = 2 Om_tE / ux at resonance.
        integrand[i] = d11int * ux / (2.0 * MPH * om_te)

    acc = np.trapezoid(integrand, etas)
    return dsdreff**(-2) * acc / Dp


def load_d11(work_dir: Path) -> float:
    data = np.loadtxt(work_dir / f"{CASE_NAME}.out")
    return float(data[4])  # column: total D11


def test_superbanana_plateau(executable: Path, tmp_path: Path) -> None:
    """Numeric SBP D11 matches the analytic superbanana-plateau prediction."""
    work_dir = tmp_path / "superbanana_plateau"
    work_dir.mkdir()
    setup_work_dir(work_dir)

    run_neort(executable, work_dir)

    d11_code = load_d11(work_dir)
    assert d11_code > 0.0, "SBP D11 is zero: ell=0 resonance not found"

    params = parse_field_params(work_dir)
    d11_pred = d11_analytic(params)

    rel_err = abs(d11_code - d11_pred) / d11_pred
    assert rel_err < D11_TOLERANCE, (
        f"SBP D11 = {d11_code:.4e} deviates from analytic prediction "
        f"{d11_pred:.4e} by {rel_err * 100:.1f}% "
        f"(tolerance: {D11_TOLERANCE * 100:.0f}%)"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
