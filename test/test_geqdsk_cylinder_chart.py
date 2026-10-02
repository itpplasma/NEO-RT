"""Independent finite-aspect circular flux/straight-angle reference."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"python"))
SPEC = importlib.util.spec_from_file_location("geqdsk_cylinder_chart", ROOT/"python/geqdsk_cylinder_chart.py")
chart = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(chart)


def test_circular_flux_and_signed_q_against_exact_finite_aspect_formula():
    R0, c, F, radius = 165., 800., -2960000., 27.
    R = np.linspace(R0-65, R0+65, 65)
    Z = np.linspace(-65., 65., 65)
    RR, ZZ = np.meshgrid(R, Z)
    data = dict(R=R, Z=Z, psi=-c*((RR-R0)**2+ZZ**2)/2,
                F=np.full(65, F), psi_axis=0., psi_edge=-c*40**2/2,
                Raxis=R0, Zaxis=0.)
    eq = chart.Equilibrium(data)
    theta, Rc, Zc, q, chi, rho = eq.contour(-c*radius**2/2, [R0+radius, 0.], 512)
    expected_q = F/(c*np.sqrt(R0**2-radius**2))
    expected_flux = abs(F)*(R0-np.sqrt(R0**2-radius**2))
    np.testing.assert_allclose(q, expected_q, rtol=1e-11)
    np.testing.assert_allclose(eq.toroidal_flux_per_radian(chi, rho), expected_flux, rtol=1e-11)
    np.testing.assert_allclose((Rc-R0)**2+Zc**2, radius**2, rtol=1e-11)
    assert np.all(np.diff(theta) > 0)
    # Closed-form angle derivative dtheta/dchi shows native theta is not
    # geometric minor-radius angle even on a concentric finite-aspect circle.
    expected_derivative = np.sqrt(R0**2-radius**2)/(R0+radius*np.cos(chi))
    numerical_derivative = (theta[2:]-theta[:-2])/(chi[2:]-chi[:-2])
    np.testing.assert_allclose(numerical_derivative, expected_derivative[1:-1], rtol=2e-5)


def test_variable_F_true_area_flux_and_legacy_label_against_closed_form():
    R0, c, F0, alpha, radius = 165., 800., -2960000., .08, 27.
    R = np.linspace(R0-65, R0+65, 65)
    Z = np.linspace(-65., 65., 65)
    RR, ZZ = np.meshgrid(R, Z)
    psi_edge = -c*40**2/2
    data = dict(R=R, Z=Z, psi=-c*((RR-R0)**2+ZZ**2)/2,
                F=F0+alpha*np.linspace(0., psi_edge, 65),
                psi_axis=0., psi_edge=psi_edge, Raxis=R0, Zaxis=0.)
    eq = chart.Equilibrium(data)
    psi_surface = -c*radius**2/2
    _, _, _, _, chi, rho = eq.contour(psi_surface, [R0+radius, 0.], 512)
    s = np.sqrt(R0**2-radius**2)
    # Exact angular integral: integral dchi/(R0+rho*cos(chi))=2*pi/sqrt(R0^2-rho^2).
    # Its independent radial primitives are I1=R0-s and I3=(R0-s)^2*(2*R0+s)/3.
    I1, I3 = R0-s, (R0-s)**2*(2*R0+s)/3
    expected_true = abs(F0*I1-alpha*c*I3/2)
    expected_legacy = abs((F0+alpha*psi_surface)*I1)
    true_flux = eq.toroidal_flux_per_radian(chi, rho, nradial=64)
    legacy_flux = eq.toroidal_flux_per_radian(chi, rho, nradial=64,
                                            boundary_F=eq.F(psi_surface))
    np.testing.assert_allclose(true_flux, expected_true, rtol=2e-11)
    np.testing.assert_allclose(legacy_flux, expected_legacy, rtol=2e-11)
    assert legacy_flux > true_flux
    assert (legacy_flux-true_flux)/true_flux > 1e-3
