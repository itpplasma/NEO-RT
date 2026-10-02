"""Periodic-cylinder harmonics on an explicitly supplied toroidal chart.

Native orthonormal components use (er, etheta, ez), Gaussian units, and
exp(i*(m*theta+n*z/R0)). The chart is R(r,theta), Z(r,theta),
phi=z/R0+nu(r,theta). The radius is a supplied flux label; this module never
identifies it with geometric minor radius. Potentials and electric fields
are one-forms; magnetic fields are flux two-forms.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PPoly


PHASE_CONVENTION = "exp(i*(m*theta+n*z/R0))"
RADIAL_CONVENTION = "sqrt(2*abs(Phi_tor_per_radian/Btf))"
C_CGS = 29979245800.0


def deformation(r, R, R0, R_r, R_theta, Z_r, Z_theta, nu_r=0, nu_theta=0):
    """Return the orthonormal deformation M and its positive determinant.

    Lengths are cm, nu is radians, and radial derivatives use the native
    flux radius in cm. Counter-clockwise poloidal charts have negative J
    with phi=z/R0 and are rejected: their native convention must be resolved
    explicitly rather than changing the mode or field sign silently.
    """
    values = np.broadcast_arrays(r, R, R_r, R_theta, Z_r, Z_theta, nu_r, nu_theta)
    if not np.isfinite(R0) or R0 <= 0:
        raise ValueError("R0 must be positive and finite")
    if any(not np.all(np.isfinite(value)) for value in values):
        raise ValueError("geometry contains non-finite values")
    r, R, R_r, R_theta, Z_r, Z_theta, nu_r, nu_theta = values
    if np.any(r <= 0) or np.any(R <= 0):
        raise ValueError("native r and physical R must be positive")
    M = np.zeros(r.shape + (3, 3))
    M[..., 0, 0], M[..., 0, 1] = R_r, R_theta / r
    M[..., 1, 0], M[..., 1, 1] = R * nu_r, R * nu_theta / r
    M[..., 1, 2] = R / R0
    M[..., 2, 0], M[..., 2, 1] = Z_r, Z_theta / r
    J = np.linalg.det(M)
    if np.any(J <= 0) or np.any(~np.isfinite(J)):
        raise ValueError("chart must have a finite positive oriented Jacobian")
    return M, J


def push_one_form(M, value):
    """Transform physical A [G cm] or E [statV/cm] by M^{-T}."""
    value = np.asarray(value)
    if value.shape[-1:] != (3,) or not np.all(np.isfinite(value)):
        raise ValueError("one-form must have three finite components")
    return np.linalg.solve(np.swapaxes(M, -1, -2), value[..., None])[..., 0]


def push_magnetic(M, J, value):
    """Piola transform a magnetic field [G], preserving its flux/divergence."""
    value = np.asarray(value)
    if value.shape[-1:] != (3,) or not np.all(np.isfinite(value)):
        raise ValueError("magnetic field must have three finite components")
    return np.einsum("...ij,...j->...i", M, value) / np.asarray(J)[..., None]


def pull_one_form(M, value):
    """Return native components of a toroidal one-form."""
    return np.einsum("...ji,...j->...i", M, value)


def pull_magnetic(M, J, value):
    """Return native magnetic components from toroidal magnetic flux."""
    rhs = np.asarray(value) * np.asarray(J)[..., None]
    return np.linalg.solve(M, rhs[..., None])[..., 0]


@dataclass
class NativeMode:
    r: np.ndarray
    electric: np.ndarray
    magnetic: np.ndarray
    m: int
    n: int
    R0: float

    @classmethod
    def from_eb(cls, path, *, m, n, R0, r_min=None, r_max=None):
        """Read KiLCA EB.dat on one continuous, explicitly selected domain.

        Columns are r, Er, Etheta, Ez, Br, Btheta, Bz (complex). Zone
        interfaces can contain two different one-sided fields at one radius;
        select a continuous interval or reject them, never average them.
        """
        data = np.loadtxt(path, ndmin=2)
        if data.shape[1] != 13 or not np.all(np.isfinite(data)):
            raise ValueError("EB.dat must contain 13 finite columns")
        order = np.argsort(data[:, 0], kind="stable")
        data = data[order]
        for bound in (r_min, r_max):
            if bound is not None and (not np.isfinite(bound)
                                      or bound < data[0, 0] or bound > data[-1, 0]):
                raise ValueError("selected radial bound lies outside EB.dat")
        if r_min is not None:
            data = data[data[:, 0] >= r_min]
        if r_max is not None:
            data = data[data[:, 0] <= r_max]
        # Zone boundaries may repeat a radius. Never discard conflicting rows.
        r, indices, counts = np.unique(data[:, 0], return_index=True, return_counts=True)
        for start, count in zip(indices, counts):
            if not np.all(data[start:start + count] == data[start]):
                raise ValueError("EB.dat repeats a radius with conflicting fields")
        data = data[indices]
        if r.size < 4 or np.any(r <= 0):
            raise ValueError("EB.dat needs at least four distinct positive radii")
        if m != int(m) or n != int(n) or not np.isfinite(R0) or R0 <= 0:
            raise ValueError("integer harmonics and positive finite R0 required")
        fields = data[:, 1::2] + 1j * data[:, 2::2]
        return cls(r, fields[:, :3], fields[:, 3:], int(m), int(n), float(R0))

    def interpolate(self, r):
        """Interpolate without extrapolating beyond the measured radial range."""
        r = np.asarray(r)
        if not np.all(np.isfinite(r)) or np.any(r < self.r[0]) or np.any(r > self.r[-1]):
            raise ValueError("requested flux radius lies outside EB.dat")
        E = CubicSpline(self.r, self.electric, axis=0)(r)
        B = CubicSpline(self.r, self.magnetic, axis=0)(r)
        return E, B


@dataclass
class RadialGauge:
    """A_r=0 reconstruction, retaining the supplied radial-B residual.

    Anchors: Az(r_min)=0 and Atheta(r_min)=i*Br(r_min)/(n/R0).
    Btheta and r times the interpolating Bz spline are integrated exactly,
    not by finite differences. Only divergence-free input satisfies Br
    equation; the diagnostic must pass before exporting a physical potential.
    """

    mode: NativeMode

    def __post_init__(self):
        if self.mode.n == 0:
            raise ValueError("A_r=0 anchoring requires a nonzero toroidal harmonic")
        r, B = self.mode.r, self.mode.magnetic
        self._int_theta = CubicSpline(r, B[:, 1]).antiderivative()
        # Multiply the Bz polynomial by r before integrating: interpolating
        # nodal r*Bz independently would change the between-node curl Bz.
        Bz = CubicSpline(r, B[:, 2])
        coefficients = np.zeros((5, r.size-1), dtype=complex)
        coefficients[:-1] += Bz.c
        coefficients[1:] += r[:-1]*Bz.c
        self._int_z = PPoly(coefficients, r).antiderivative()

    def evaluate(self, r):
        _, B = self.mode.interpolate(r)
        r = np.asarray(r)
        anchor = self.mode.r[0]
        kz = self.mode.n / self.mode.R0
        Az = -(self._int_theta(r) - self._int_theta(anchor))
        Atheta = (self._int_z(r) - self._int_z(anchor)
                  + anchor * 1j * self.mode.magnetic[0, 0] / kz) / r
        A = np.stack((np.zeros_like(Az), Atheta, Az), axis=-1)
        reconstructed_Br = 1j * self.mode.m * Az / r - 1j * kz * Atheta
        return A, reconstructed_Br - B[..., 0]

    def native_curl(self, r):
        """Return curl A in the cylinder, including its supplied-Br mismatch."""
        _, B = self.mode.interpolate(r)
        _, residual = self.evaluate(r)
        B = np.array(B, copy=True)
        B[..., 0] += residual
        return B

    def relative_radial_residual(self):
        # Check nodes and between nodes independently of geometry conversion.
        r = np.sort(np.concatenate((self.mode.r, (self.mode.r[1:] + self.mode.r[:-1]) / 2)))
        _, residual = self.evaluate(r)
        _, B = self.mode.interpolate(r)
        scale = np.max(np.abs(B[:, 0]))
        error = np.max(np.abs(residual))
        return float(error / scale) if scale > 0 else (0.0 if error == 0 else float("inf"))


@dataclass
class AxialGauge:
    """Paired laboratory potentials for exp(i*m*theta+i*kz*z-i*omega*t).

    For n!=0, Az=0 fixes Phi=i*Ez/kz, Ar=-i*Btheta/kz and
    Atheta=i*Br/kz. The unforced Bz and both remaining electric equations
    must pass independent checks; no static potential is inferred from an
    arbitrary electric-field component.
    """

    mode: NativeMode
    omega_lab: complex

    def __post_init__(self):
        if self.mode.n == 0 or not np.isfinite(self.omega_lab):
            raise ValueError("axial gauge needs nonzero n and finite laboratory omega")
        self._electric = CubicSpline(self.mode.r, self.mode.electric, axis=0)
        self._magnetic = CubicSpline(self.mode.r, self.mode.magnetic, axis=0)

    def evaluate(self, r):
        E, B = self.mode.interpolate(r)
        kz = self.mode.n / self.mode.R0
        A = np.stack((-1j*B[..., 1]/kz, 1j*B[..., 0]/kz,
                      np.zeros_like(B[..., 0])), axis=-1)
        return A, 1j*E[..., 2]/kz

    def diagnostics(self):
        r = np.sort(np.concatenate((self.mode.r, (self.mode.r[1:] + self.mode.r[:-1])/2)))
        E, B = self.mode.interpolate(r)
        dE, dB = self._electric(r, 1), self._magnetic(r, 1)
        m, kz, omega = self.mode.m, self.mode.n/self.mode.R0, self.omega_lab
        A, Phi = self.evaluate(r)
        curl_A = np.stack((B[:, 0], B[:, 1],
                           1j*(dB[:, 0]+B[:, 0]/r)/kz-m*B[:, 1]/(kz*r)), axis=-1)
        gradient_Phi = np.stack((1j*dE[:, 2]/kz, 1j*m*Phi/r, 1j*kz*Phi), axis=-1)
        from_potentials_E = -gradient_Phi + 1j*omega*A/C_CGS
        curl_E = np.stack((1j*m*E[:, 2]/r-1j*kz*E[:, 1],
                           1j*kz*E[:, 0]-dE[:, 2],
                           dE[:, 1]+E[:, 1]/r-1j*m*E[:, 0]/r), axis=-1)
        faraday = curl_E-1j*omega*B/C_CGS
        curl_term_sum = np.stack((abs(m*E[:, 2]/r)+abs(kz*E[:, 1]),
                                 abs(kz*E[:, 0])+abs(dE[:, 2]),
                                 abs(dE[:, 1])+abs(E[:, 1]/r)+abs(m*E[:, 0]/r)), axis=-1)
        divergence_B = dB[:, 0]+B[:, 0]/r+1j*m*B[:, 1]/r+1j*kz*B[:, 2]
        def relative(error, reference):
            scale, value = np.max(np.abs(reference)), np.max(np.abs(error))
            return float(value/scale) if scale > 0 else (0. if value == 0 else float("inf"))
        return {"relative_curl_A_minus_B": relative(curl_A-B, B),
                "relative_E_from_potentials_minus_E": relative(from_potentials_E-E, E),
                # Scale the Faraday residual by both sides. This records large
                # cancellation at low lab frequency rather than hiding it.
                "relative_Faraday_residual": relative(faraday,
                                                       np.maximum(abs(curl_E), abs(omega*B/C_CGS))),
                "relative_div_B": relative(divergence_B, np.abs(dB[:, 0])+np.abs(B[:, 0]/r)
                                            +np.abs(m*B[:, 1]/r)+np.abs(kz*B[:, 2])),
                "Faraday_term_conditioning": relative(curl_term_sum, omega*B/C_CGS),
                "max_abs_Faraday_statV_per_cm2": float(np.max(abs(faraday))),
                "max_abs_inductive_curl_statV_per_cm2": float(np.max(abs(omega*B/C_CGS)))}

    def integrated_electric_residual(self):
        """Check radial E using a primitive instead of differentiating Ez data."""
        r = np.sort(np.concatenate((self.mode.r, (self.mode.r[1:]+self.mode.r[:-1])/2)))
        _, Phi = self.evaluate(r)
        kz = self.mode.n/self.mode.R0
        derivative = -self.mode.electric[:, 0]+self.omega_lab*self.mode.magnetic[:, 1]/(kz*C_CGS)
        primitive = CubicSpline(self.mode.r, derivative).antiderivative()
        Phi_from_integral = Phi[0]+primitive(r)-primitive(r[0])
        scale, error = np.max(abs(Phi)), np.max(abs(Phi_from_integral-Phi))
        return float(error/scale) if scale > 0 else (0. if error == 0 else float("inf"))


def transform_mode(mode, geometry, *, max_relative_br_residual):
    """Return a single-n toroidal amplitude at explicit chart sample points.

    Geometry requires arrays r, theta, R, Z, nu and their supplied radial and
    angular derivatives; all arrays broadcast to the same shape. The phase
    exp(i*(m*theta-n*nu)) converts native harmonics to exp(i*n*phi).
    E is a diagnostic vector only: finite-frequency input supplies no static
    electrostatic potential for the NTV Hamiltonian.
    """
    if not np.isfinite(max_relative_br_residual) or max_relative_br_residual < 0:
        raise ValueError("a finite nonnegative radial-B residual limit is required")
    if "R0_cm" in geometry and not np.isclose(mode.R0, geometry["R0_cm"], rtol=1e-12, atol=0):
        raise ValueError("native mode R0 differs from the supplied chart")
    gauge = RadialGauge(mode)
    relative_error = gauge.relative_radial_residual()
    if relative_error > max_relative_br_residual:
        raise ValueError(f"native radial-B reconstruction residual {relative_error:.6g} "
                         f"exceeds {max_relative_br_residual:.6g}")
    r, theta, R, Z, nu = np.broadcast_arrays(*(geometry[key] for key in
                                              ("r", "theta", "R", "Z", "nu")))
    if any(not np.all(np.isfinite(value)) for value in (theta, Z, nu)):
        raise ValueError("geometry contains non-finite coordinates or phase")
    M, J = deformation(r, R, mode.R0, *(geometry[key] for key in
                       ("R_r", "R_theta", "Z_r", "Z_theta", "nu_r", "nu_theta")))
    A, residual = gauge.evaluate(r)
    E, B = mode.interpolate(r)
    phase = np.exp(1j * (mode.m * theta - mode.n * nu))[..., None]
    return {"r_cm": r, "theta_rad": theta, "R_cm": R, "Z_cm": Z,
            "A_Gcm": push_one_form(M, A) * phase,
            "B_G": push_magnetic(M, J, B) * phase,
            "E_statVcm_diagnostic": push_one_form(M, E) * phase,
            "native_Br_residual_G": residual,
            "max_relative_Br_residual": relative_error,
            "n_tor": mode.n, "source_m": mode.m, "R0_cm": mode.R0}


def transform_electromagnetic_mode(mode, geometry, *, omega_lab,
                                  max_relative_B_residual, max_relative_E_residual,
                                  max_relative_Faraday_residual, max_relative_div_B):
    """Push paired axial-gauge A/Phi after validating supplied laboratory E/B."""
    if "R0_cm" in geometry and not np.isclose(mode.R0, geometry["R0_cm"], rtol=1e-12, atol=0):
        raise ValueError("native mode R0 differs from the supplied chart")
    for limit in (max_relative_B_residual, max_relative_E_residual,
                  max_relative_Faraday_residual, max_relative_div_B):
        if not np.isfinite(limit) or limit < 0:
            raise ValueError("finite nonnegative E/B residual limits are required")
    gauge = AxialGauge(mode, omega_lab)
    diagnostics = gauge.diagnostics()
    if diagnostics["relative_curl_A_minus_B"] > max_relative_B_residual:
        raise ValueError(f"curl A differs from supplied B: {diagnostics}")
    if diagnostics["relative_E_from_potentials_minus_E"] > max_relative_E_residual:
        raise ValueError(f"paired potentials differ from supplied laboratory E: {diagnostics}")
    if diagnostics["relative_Faraday_residual"] > max_relative_Faraday_residual:
        raise ValueError(f"supplied laboratory E/B violate Faraday: {diagnostics}")
    if diagnostics["relative_div_B"] > max_relative_div_B:
        raise ValueError(f"supplied magnetic field has nonzero divergence: {diagnostics}")
    r, theta, R, Z, nu = np.broadcast_arrays(*(geometry[key] for key in
                                              ("r", "theta", "R", "Z", "nu")))
    if any(not np.all(np.isfinite(value)) for value in (theta, Z, nu)):
        raise ValueError("geometry contains non-finite coordinates or phase")
    M, J = deformation(r, R, mode.R0, *(geometry[key] for key in
                       ("R_r", "R_theta", "Z_r", "Z_theta", "nu_r", "nu_theta")))
    A, Phi = gauge.evaluate(r)
    E, B = mode.interpolate(r)
    phase = np.exp(1j*(mode.m*theta-mode.n*nu))
    return {"r_cm": r, "theta_rad": theta, "R_cm": R, "Z_cm": Z,
            "A_Gcm": push_one_form(M, A)*phase[..., None],
            "Phi_statV": Phi*phase, "B_G": push_magnetic(M, J, B)*phase[..., None],
            "E_statVcm": push_one_form(M, E)*phase[..., None],
            "omega_lab_rad_per_s": omega_lab, "n_tor": mode.n,
            "source_m": mode.m, "R0_cm": mode.R0,
            "native_EB_diagnostics_json": json.dumps(diagnostics, sort_keys=True)}


def laboratory_omega(path, *, m, n):
    """Read KiLCA mode_data.dat; frequencies are Hz, fields are laboratory-frame."""
    data = np.loadtxt(path, comments="%", ndmin=2)
    if data.shape != (4, 2) or not np.all(np.isfinite(data)):
        raise ValueError("mode_data.dat must contain four finite two-column rows")
    if not np.array_equal(data[0], [m, n]):
        raise ValueError("mode_data.dat harmonics differ from the requested mode")
    omega = 2*np.pi*complex(*data[1])
    return omega.real if omega.imag == 0 else omega


def read_geometry(path):
    """Read an explicit native chart NPZ; require its coordinate provenance."""
    with np.load(Path(path), allow_pickle=False) as data:
        for key, value in (("radial_convention", RADIAL_CONVENTION),
                           ("phase_convention", PHASE_CONVENTION),
                           ("length_units", "cm")):
            if key not in data or str(data[key].item()) != value:
                raise ValueError(f"geometry must declare {key}={value}")
        if "provenance" not in data or not str(data["provenance"].item()).strip():
            raise ValueError("geometry must identify its native chart provenance")
        keys = ("r", "theta", "R", "Z", "nu", "R_r", "R_theta", "Z_r", "Z_theta",
                "nu_r", "nu_theta")
        if any(key not in data for key in keys):
            raise ValueError("geometry is missing chart coordinates or derivatives")
        if "R0_cm" not in data or data["R0_cm"].shape != ():
            raise ValueError("geometry must declare its scalar machine R0_cm")
        result = {key: np.asarray(data[key]) for key in keys}
        result["R0_cm"] = float(data["R0_cm"])
        return result


class SurfaceChart:
    """Continuous, integrable interpolation of a supplied r/theta chart.

    Coordinates are primary. Their derivatives come from the same splines,
    so curl/Piola identities can be checked between sample points without
    interpolating mutually inconsistent coordinate and derivative tables.
    """

    def __init__(self, geometry):
        from scipy.optimize import root
        self._root = root
        r, theta = np.asarray(geometry["r"]), np.asarray(geometry["theta"])
        if r.ndim != 2 or theta.shape != r.shape:
            raise ValueError("surface chart requires a separable two-dimensional r/theta grid")
        self.r, self.theta = r[:, 0], theta[0]
        if not np.all(r == self.r[:, None]) or not np.all(theta == self.theta[None, :]):
            raise ValueError("surface chart r/theta grid is not separable")
        if (self.r.size < 4 or self.theta.size < 4 or np.any(np.diff(self.r) <= 0)
                or np.any(np.diff(self.theta) <= 0) or self.theta[-1] >= self.theta[0]+2*np.pi):
            raise ValueError("surface chart needs increasing radial and periodic angular nodes")
        self.R0 = float(geometry["R0_cm"])
        th = np.append(self.theta, self.theta[0]+2*np.pi)
        self._angular = {}
        for key in ("R", "Z", "nu"):
            value = np.asarray(geometry[key])
            if value.shape != r.shape or not np.all(np.isfinite(value)):
                raise ValueError("surface coordinates must be finite and match the chart grid")
            self._angular[key] = CubicSpline(th, np.concatenate((value, value[:, :1]), axis=1),
                                              axis=1, bc_type="periodic")

    def geometry_at(self, r, theta):
        if not np.isfinite(r) or not self.r[0] <= r <= self.r[-1] or not np.isfinite(theta):
            raise ValueError("point lies outside the supplied surface chart")
        angle = (theta-self.theta[0]) % (2*np.pi)+self.theta[0]
        geometry = dict(r=r, theta=theta, R0_cm=self.R0)
        for key, angular in self._angular.items():
            radial = CubicSpline(self.r, angular(angle))
            geometry[key] = float(radial(r))
            geometry[key+"_r"] = float(radial(r, 1))
            geometry[key+"_theta"] = float(CubicSpline(self.r, angular(angle, 1))(r))
        return geometry

    def invert(self, R, Z, *, guess):
        def residual(point):
            geometry = self.geometry_at(*point)
            return [geometry["R"]-R, geometry["Z"]-Z]
        solution = self._root(residual, guess, tol=1e-10)
        if not solution.success or np.max(np.abs(residual(solution.x))) > 1e-8:
            raise ValueError("point does not invert within the supplied surface chart")
        return tuple(solution.x)
