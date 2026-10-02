"""Independent geometry/Maxwell oracles for the periodic-cylinder adapter."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("periodic_cylinder", ROOT / "python/periodic_cylinder.py")
adapter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter
SPEC.loader.exec_module(adapter)


def circle(r, theta, R0=165.0, shear=0.0):
    r, theta = np.broadcast_arrays(r, theta)
    c, s = np.cos(theta), np.sin(theta)
    return dict(r=r, theta=theta, R=R0 + r * c, Z=-r * s,
                nu=shear * r * c / R0, R_r=c, R_theta=-r * s,
                Z_r=-s, Z_theta=-r * c,
                nu_r=shear * c / R0, nu_theta=-shear * r * s / R0)


def matrix(g, R0=165.0):
    return adapter.deformation(g["r"], g["R"], R0, *(g[key] for key in
                               ("R_r", "R_theta", "Z_r", "Z_theta", "nu_r", "nu_theta")))


def test_circle_components_against_closed_form():
    g = circle(np.array([10., 25., 40.]), np.array([0.2, 1.3, 3.4]))
    M, J = matrix(g)
    A = np.array([[1+2j, 3-4j, 5+6j]] * 3)
    B = A * (2-3j)
    c, s, metric = np.cos(g["theta"]), np.sin(g["theta"]), 165 / g["R"]
    expected_A = np.stack((A[:, 0]*c-A[:, 1]*s, metric*A[:, 2],
                           -A[:, 0]*s-A[:, 1]*c), axis=-1)
    expected_B = np.stack((metric*(B[:, 0]*c-B[:, 1]*s), B[:, 2],
                           metric*(-B[:, 0]*s-B[:, 1]*c)), axis=-1)
    np.testing.assert_allclose(J, g["R"] / 165, rtol=2e-15)
    np.testing.assert_allclose(adapter.push_one_form(M, A), expected_A, rtol=2e-15)
    np.testing.assert_allclose(adapter.push_magnetic(M, J, B), expected_B, rtol=2e-15)
    np.testing.assert_allclose(adapter.pull_one_form(M, expected_A), A, rtol=2e-15)
    np.testing.assert_allclose(adapter.pull_magnetic(M, J, expected_B), B, rtol=2e-15)


def test_general_chart_preserves_line_and_flux_pairings():
    g = circle(37., 0.73, shear=0.4)
    M, J = matrix(g)
    # Independently differentiate the Cartesian embedding, not the matrix.
    q = np.array([37., 0.73, 12.])
    def embedding(q):
        r, theta, z = q
        R = 165 + r*np.cos(theta)
        phi = z/165 + 0.4*r*np.cos(theta)/165
        return np.array([R*np.cos(phi), R*np.sin(phi), -r*np.sin(theta)])
    tangents = []
    for v in (np.array([0.2, 0.3/37., -0.4]), np.array([0.4, -0.1/37., 0.2])):
        h = 1e-4
        tangents.append((embedding(q+h*v)-embedding(q-h*v))/(2*h))
    phi = q[2]/165 + g["nu"]
    basis = np.array([[np.cos(phi), -np.sin(phi), 0],
                      [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
    A = np.array([0.3+0.2j, -0.7j, 1.2])
    B = np.array([0.9j, 0.4, -0.2+0.1j])
    At = basis @ adapter.push_one_form(M, A)
    Bt = basis @ adapter.push_magnetic(M, J, B)
    native_v = np.array([0.2, 0.3, -0.4])
    native_w = np.array([0.4, -0.1, 0.2])
    np.testing.assert_allclose(At @ tangents[0], A @ native_v, rtol=2e-9, atol=2e-10)
    np.testing.assert_allclose(Bt @ np.cross(*tangents), B @ np.cross(native_v, native_w),
                               rtol=2e-9, atol=2e-10)


def polynomial_mode(n=2):
    r = np.linspace(3., 60., 17)
    m, R0 = 3, 165.
    a, b = 0.002+0.003j, -0.004+0.001j
    kz = 2/R0
    B = np.stack((1j*m*a*r-1j*kz*b*r, -2*a*r, np.full_like(r, 2*b, dtype=complex)), axis=-1)
    return adapter.NativeMode(r, np.zeros_like(B), B, m, n, R0), a, b


def test_radial_gauge_against_polynomial_potential():
    mode, a, b = polynomial_mode()
    gauge = adapter.RadialGauge(mode)
    r = np.array([3., 7.3, 25.1, 59.2])
    A, residual = gauge.evaluate(r)
    r0, kz = mode.r[0], mode.n/mode.R0
    expected_Az = a*(r*r-r0*r0)
    expected_Atheta = b*r-mode.m*a*r0*r0/(kz*r)
    np.testing.assert_allclose(A[:, 2], expected_Az, atol=5e-14)
    np.testing.assert_allclose(A[:, 1], expected_Atheta, atol=5e-14)
    np.testing.assert_allclose(residual, 0, atol=2e-15)
    assert gauge.relative_radial_residual() < 1e-13


def test_radial_gauge_integrates_the_same_interpolated_Bz():
    r = np.linspace(3., 60., 17)
    c, kz = .0002+.0003j, 2/165.
    Bz = c*r**3
    B = np.stack((-1j*kz*c*r**4/5, np.zeros(r.size), Bz), axis=-1)
    mode = adapter.NativeMode(r, np.zeros_like(B), B, 3, 2, 165.)
    gauge = adapter.RadialGauge(mode)
    samples = np.array([7.3, 27.2, 59.3])
    A, _ = gauge.evaluate(samples)
    np.testing.assert_allclose(A[:, 1], c*samples**4/5, rtol=3e-15, atol=1e-12)
    # A separately fitted spline of nodal r*Bz does not satisfy this quartic
    # potential between nodes; Bz itself is an exactly interpolated cubic.
    h = 1e-4
    plus, _ = gauge.evaluate(samples+h)
    minus, _ = gauge.evaluate(samples-h)
    curl_z = ((samples+h)*plus[:, 1]-(samples-h)*minus[:, 1])/(2*h*samples)
    np.testing.assert_allclose(curl_z, c*samples**3, rtol=1e-9)


@pytest.mark.parametrize("shear", [0., 0.4])
def test_curl_and_divergence_in_physical_cylindrical_coordinates(shear):
    mode, _, _ = polynomial_mode()
    gauge = adapter.RadialGauge(mode)
    def fields(R, Z):
        r, theta = np.hypot(R-165, Z), np.arctan2(-Z, R-165)
        g = circle(r, theta, shear=shear)
        M, J = matrix(g)
        A, _ = gauge.evaluate(r)
        _, B = mode.interpolate(r)
        phase = np.exp(1j*(mode.m*theta-mode.n*g["nu"]))
        return adapter.push_one_form(M, A)*phase, adapter.push_magnetic(M, J, B)*phase
    R, Z, h, n = 180., -23., 1e-4, mode.n
    A, B = fields(R, Z)
    ApR, BpR = fields(R+h, Z)
    AmR, BmR = fields(R-h, Z)
    ApZ, BpZ = fields(R, Z+h)
    AmZ, BmZ = fields(R, Z-h)
    dR, dZ = (ApR-AmR)/(2*h), (ApZ-AmZ)/(2*h)
    curl = np.array([1j*n*A[2]/R-dZ[1], dZ[0]-dR[2],
                     dR[1]+A[1]/R-1j*n*A[0]/R])
    np.testing.assert_allclose(curl, B, rtol=3e-8, atol=2e-10)
    divergence = ((R+h)*BpR[0]-(R-h)*BmR[0])/(2*h*R) + 1j*n*B[1]/R + (BpZ[2]-BmZ[2])/(2*h)
    assert abs(divergence) < 2e-10


def test_wrong_mode_sign_and_non_solenoidal_input_are_rejected():
    mode, _, _ = polynomial_mode(n=-2)
    g = circle(20., 1.1)
    with pytest.raises(ValueError, match="radial-B reconstruction residual"):
        adapter.transform_mode(mode, g, max_relative_br_residual=1e-8)
    mode, _, _ = polynomial_mode()
    mode.magnetic[8, 0] += 0.1
    with pytest.raises(ValueError, match="radial-B reconstruction residual"):
        adapter.transform_mode(mode, g, max_relative_br_residual=1e-8)


def test_map_orientation_axis_and_extrapolation_are_rejected():
    g = circle(20., 1.1)
    g["Z_r"] *= -1
    g["Z_theta"] *= -1
    with pytest.raises(ValueError, match="positive oriented Jacobian"):
        matrix(g)
    with pytest.raises(ValueError, match="must be positive"):
        matrix(circle(0., 1.1))
    mode, _, _ = polynomial_mode()
    with pytest.raises(ValueError, match="outside EB.dat"):
        mode.interpolate(61.)
    mode.n = 0
    with pytest.raises(ValueError, match="nonzero toroidal"):
        adapter.RadialGauge(mode)


def test_eb_schema_and_conflicting_duplicates(tmp_path):
    rows = np.arange(4*13, dtype=float).reshape(4, 13)
    rows[:, 0] = [1., 2., 3., 4.]
    path = tmp_path/"EB.dat"
    np.savetxt(path, rows)
    mode = adapter.NativeMode.from_eb(path, m=7, n=2, R0=165)
    np.testing.assert_equal(mode.electric[0], [1+2j, 3+4j, 5+6j])
    np.testing.assert_equal(mode.magnetic[0], [7+8j, 9+10j, 11+12j])
    np.savetxt(path, np.vstack((rows, rows[1])))
    assert adapter.NativeMode.from_eb(path, m=7, n=2, R0=165).r.size == 4
    rows[1, 0] = 1.
    np.savetxt(path, rows)
    with pytest.raises(ValueError, match="conflicting fields"):
        adapter.NativeMode.from_eb(path, m=7, n=2, R0=165)
    # A domain beyond the conflicting interface retains its own intact rows.
    extra = np.array(rows[-1], copy=True)
    extra[0] = 5.
    np.savetxt(path, np.vstack((rows, extra)))
    # Need four continuous rows, so add a sixth independent radial node.
    extra[0] = 6.
    with path.open("a") as stream:
        np.savetxt(stream, extra[None, :])
    assert adapter.NativeMode.from_eb(path, m=7, n=2, R0=165, r_min=3.).r.size == 4


def electromagnetic_mode():
    r = np.linspace(3., 60., 19)
    m, n, R0, omega = 3, 2, 165., 2*np.pi
    a, b, d = .002+.003j, .0004+.0001j, 2e-13-5e-13j
    kz = n/R0
    Ar, Atheta, Phi = a*r, b*r*r, d*r*r
    B = np.stack((-1j*kz*Atheta, 1j*kz*Ar, 3*b*r-1j*m*a), axis=-1)
    E = np.stack((-2*d*r+1j*omega*Ar/adapter.C_CGS,
                  -1j*m*Phi/r+1j*omega*Atheta/adapter.C_CGS, -1j*kz*Phi), axis=-1)
    return adapter.NativeMode(r, E, B, m, n, R0), a, b, d, omega


def test_paired_axial_gauge_against_independent_polynomial_fields():
    mode, a, b, d, omega = electromagnetic_mode()
    gauge = adapter.AxialGauge(mode, omega)
    r = np.array([3., 7.7, 27.2, 60.])
    A, Phi = gauge.evaluate(r)
    np.testing.assert_allclose(A, np.stack((a*r, b*r*r, np.zeros(r.size)), axis=-1), atol=1e-14)
    np.testing.assert_allclose(Phi, d*r*r, atol=1e-23)
    diagnostics = gauge.diagnostics()
    assert max(value for key, value in diagnostics.items() if key.startswith("relative_")) < 2e-13
    assert gauge.integrated_electric_residual() < 2e-13
    g = circle(r, [.2, 1.1, 2.3, 5.7], shear=.4)
    fields = adapter.transform_electromagnetic_mode(mode, g, omega_lab=omega,
        max_relative_B_residual=1e-10, max_relative_E_residual=1e-10,
        max_relative_Faraday_residual=1e-10, max_relative_div_B=1e-10)
    expected = d*r*r*np.exp(1j*(mode.m*g["theta"]-mode.n*g["nu"]))
    np.testing.assert_allclose(fields["Phi_statV"], expected, atol=1e-23)
    assert fields["omega_lab_rad_per_s"] == 2*np.pi


@pytest.mark.parametrize("wrong_omega", [-2*np.pi, 1., 2*np.pi*192916.0825356307])
def test_lab_time_sign_hertz_conversion_and_frame_mistakes_fail(wrong_omega):
    mode, _, _, _, _ = electromagnetic_mode()
    with pytest.raises(ValueError, match="supplied laboratory E"):
        adapter.transform_electromagnetic_mode(mode, circle(20., .7), omega_lab=wrong_omega,
            max_relative_B_residual=1e-10, max_relative_E_residual=1e-10,
            max_relative_Faraday_residual=1e-10, max_relative_div_B=1e-10)


def test_mode_data_frequency_contract(tmp_path):
    path = tmp_path/"mode_data.dat"
    path.write_text("%m n\n%Re(flab) Im(flab)\n%Re(fmov) Im(fmov)\n%r_res 0\n"
                    "7 2\n1 0\n192916.0825356307 0\n58.49591397161032 0\n")
    assert adapter.laboratory_omega(path, m=7, n=2) == 2*np.pi
    with pytest.raises(ValueError, match="harmonics differ"):
        adapter.laboratory_omega(path, m=-7, n=2)


def test_surface_chart_against_analytic_circle_and_coordinate_derivatives():
    g = circle(np.linspace(3., 60., 9)[:, None],
               np.arange(128)[None, :]*2*np.pi/128, shear=.4)
    g["R0_cm"] = 165.
    surface = adapter.SurfaceChart(g)
    # Nodes retain the coordinate values; off-node derivatives are checked
    # against coordinate differences rather than a second derivative table.
    node = surface.geometry_at(g["r"][4, 17], g["theta"][4, 17])
    for key in ("R", "Z", "nu"):
        np.testing.assert_allclose(node[key], g[key][4, 17], atol=3e-14, rtol=0)
    for r, theta in ((7.3, .37), (27.1, 1.42), (55.2, 5.27)):
        exact = circle(r, theta, shear=.4)
        interpolated = surface.geometry_at(r, theta)
        for key in ("R", "Z", "nu"):
            np.testing.assert_allclose(interpolated[key], exact[key], atol=3e-6, rtol=0)
            h = 1e-4
            radial = (surface.geometry_at(r+h, theta)[key]
                      - surface.geometry_at(r-h, theta)[key])/(2*h)
            angular = (surface.geometry_at(r, theta+h)[key]
                       - surface.geometry_at(r, theta-h)[key])/(2*h)
            np.testing.assert_allclose(interpolated[key+"_r"], radial, atol=3e-9, rtol=0)
            np.testing.assert_allclose(interpolated[key+"_theta"], angular, atol=3e-7, rtol=0)
        recovered = surface.invert(exact["R"], exact["Z"], guess=(r+.01, theta+.001))
        np.testing.assert_allclose(recovered, [r, theta], atol=3e-6, rtol=0)
    with pytest.raises(ValueError, match="outside the supplied surface chart"):
        surface.geometry_at(61., .3)
