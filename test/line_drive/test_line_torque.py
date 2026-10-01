"""End-to-end collapse of drive_form = 'line' onto the Boozer scalar route.

Runs neo_rt.x twice on the golden-record equilibrium with the ideal_helical
perturbation, once per drive form, and compares torque density and D11 per
orbit class.  The perturbation is built so that the Boozer scalar applies
exactly; the residual difference comes from the bounce-integral tolerance and
rotation terms of order (Om_E R / v)^2 that both routes treat differently.
"""
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "test" / "golden_record" / "input"
NAMELIST = """&params
    s = 0.5, m_t = 0.036, qs = 1.0, ms = 2.014, epsmn = 1.0e-3, m0 = 2, mph = 1,
    magdrift = .true., nopassing = .false., noshear = .false., pertfile = .false.,
    nonlin = .false., bfac = 1.0, efac = 1.0, inp_swi = 9, vsteps = 256,
    vmax_over_vth = 3.0, log_level = -1, comptorque = .true.,
    drive_form = '{form}', pert_model = 'ideal_helical'
/
"""
RTOL = 1.0e-3


def executable():
    for path in (ROOT / "build" / "neo_rt.x", Path(os.environ.get("NEO_RT_EXE", ""))):
        if path.is_file():
            return path
    found = shutil.which("neo_rt.x")
    if found is None:
        pytest.skip("neo_rt.x not found")
    return Path(found)


def run(form, work):
    work.mkdir()
    for name in ("in_file", "plasma.in", "profile.in"):
        (work / name).symlink_to(INPUT / name)
    (work / "run.in").write_text(NAMELIST.format(form=form))
    subprocess.run([str(executable()), "run"], cwd=work, check=True,
                   capture_output=True)
    torque = np.loadtxt(work / "run_torque.out", skiprows=1)[3:6]
    d11 = np.loadtxt(work / "run.out", skiprows=1)[1:4]
    return torque, d11


def test_line_matches_boozer_torque(tmp_path):
    torque_b, d11_b = run("boozer", tmp_path / "boozer")
    torque_l, d11_l = run("line", tmp_path / "line")
    assert np.all(np.abs(torque_b) > 0.0), "vacuous: zero Boozer torque"
    np.testing.assert_allclose(torque_l, torque_b, rtol=RTOL)
    np.testing.assert_allclose(d11_l, d11_b, rtol=RTOL)
