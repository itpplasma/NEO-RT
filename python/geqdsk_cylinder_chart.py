#!/usr/bin/env python3
"""Regenerate a specified clockwise straight-field-line chart from GEQDSK.

The recovered equil_r_q_psi.dat supplies its declared approximate flux radius
and outboard phase anchor. This constructs an independent chart from the matched
equilibrium; it does not claim to recover a missing historical angle table.
"""

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline, PchipInterpolator, RectBivariateSpline
from scipy.optimize import brentq, root

from periodic_cylinder import PHASE_CONVENTION, RADIAL_CONVENTION, deformation


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_geqdsk(path):
    with Path(path).open() as stream:
        header, body = stream.readline(), stream.read()
    nw, nh = map(int, header.split()[-2:])
    pattern = r"[-+]?\d*\.\d+(?:[eEdD][-+]?\d+)?"
    values = np.array([float(x.replace("D", "e").replace("d", "e"))
                       for x in re.findall(pattern, body)])
    if values.size < 20+5*nw+nw*nh:
        raise ValueError("incomplete GEQDSK arrays")
    keys = ("rdim", "zdim", "R0", "rleft", "zmid", "Raxis", "Zaxis", "psi_axis",
            "psi_edge", "Btf", "current")
    data = dict(zip(keys, values[:11]))
    data["F"] = values[20:20+nw]*1e6  # T m -> G cm
    offset = 20+4*nw
    data["psi"] = values[offset:offset+nw*nh].reshape(nh, nw)*1e8
    data["q"] = values[offset+nw*nh:offset+nw*nh+nw]
    data["R"] = 100*(data["rleft"]+np.linspace(0, data["rdim"], nw))
    data["Z"] = 100*(data["zmid"]-data["zdim"]/2+np.linspace(0, data["zdim"], nh))
    data["psi_axis"] *= 1e8
    data["psi_edge"] *= 1e8
    data["Raxis"] *= 100
    data["Zaxis"] *= 100
    return data


class Equilibrium:
    def __init__(self, data):
        self.data = data
        self.psi = RectBivariateSpline(data["R"], data["Z"], data["psi"].T, kx=3, ky=3)
        grid = np.linspace(data["psi_axis"], data["psi_edge"], data["F"].size)
        order = np.argsort(grid)
        self.F = CubicSpline(grid[order], data["F"][order])
        solution = root(lambda x: [self.psi(x[0], x[1], dx=1)[0, 0],
                                    self.psi(x[0], x[1], dy=1)[0, 0]],
                        [data["Raxis"], data["Zaxis"]])
        if not solution.success:
            raise ValueError("GEQDSK magnetic-axis refinement failed")
        self.Raxis, self.Zaxis = solution.x

    def contour(self, psi_target, phase_anchor, ntheta):
        """Clockwise geometric rays; first crossing stays inside the grid."""
        chi0 = np.arctan2(-(phase_anchor[1]-self.Zaxis), phase_anchor[0]-self.Raxis)
        chi = chi0+np.linspace(0, 2*np.pi, ntheta+1)
        radius = np.empty(chi.size)
        axis_psi = self.psi(self.Raxis, self.Zaxis)[0, 0]
        sign = np.sign(psi_target-axis_psi)
        if sign == 0:
            raise ValueError("magnetic axis has no nonsingular poloidal chart")
        for j, angle in enumerate(chi[:-1]):
            c, s = np.cos(angle), -np.sin(angle)
            bounds = []
            for origin, direction, grid in ((self.Raxis, c, self.data["R"]),
                                             (self.Zaxis, s, self.data["Z"])):
                if abs(direction) > 1e-14:
                    edge = grid[-1] if direction > 0 else grid[0]
                    bounds.append((edge-origin)/direction)
            probe = np.linspace(0., .999*min(bounds), 160)
            values = sign*(self.psi(self.Raxis+probe*c, self.Zaxis+probe*s, grid=False)-psi_target)
            crossings = np.flatnonzero(values > 0)
            if not crossings.size or crossings[0] == 0:
                raise ValueError("flux surface has no first ray crossing inside GEQDSK")
            k = crossings[0]
            radius[j] = brentq(lambda rho: sign*(self.psi(self.Raxis+rho*c,
                              self.Zaxis+rho*s)[0, 0]-psi_target), probe[k-1], probe[k], xtol=1e-10)
        radius[-1] = radius[0]
        R = self.Raxis+radius*np.cos(chi)
        Z = self.Zaxis-radius*np.sin(chi)
        periodic_radius = CubicSpline(chi, radius, bc_type="periodic")
        radial_derivative = periodic_radius(chi, 1)
        dR = radial_derivative*np.cos(chi)-radius*np.sin(chi)
        dZ = -radial_derivative*np.sin(chi)-radius*np.cos(chi)
        dl = np.hypot(dR, dZ)
        psi_R = self.psi(R, Z, dx=1, grid=False)
        psi_Z = self.psi(R, Z, dy=1, grid=False)
        Bp_clockwise = (-psi_Z*dR+psi_R*dZ)/(R*dl)
        if np.any(Bp_clockwise == 0):
            raise ValueError("poloidal field vanishes on the requested contour")
        dphi_dchi = self.F(psi_target)*dl/(R*R*Bp_clockwise)
        integral = cumulative_trapezoid(dphi_dchi, chi, initial=0.)
        q = integral[-1]/(2*np.pi)
        theta = integral/q
        if not np.all(np.diff(theta) > 0):
            raise ValueError("straight poloidal angle is not monotone")
        return theta, R, Z, q, chi, radius

    def toroidal_flux_per_radian(self, chi, radius, *, nradial=32, boundary_F=None):
        """Independent area integral of Bphi=F(psi)/R, in G cm^2/radian."""
        nodes, weights = np.polynomial.legendre.leggauss(nradial)
        rho = radius[:, None]*(nodes[None, :]+1)/2
        R = self.Raxis+rho*np.cos(chi[:, None])
        Z = self.Zaxis-rho*np.sin(chi[:, None])
        psi = self.psi(R, Z, grid=False)
        Bphi = (self.F(psi) if boundary_F is None else boundary_F)/R
        radial_integral = np.sum(Bphi*rho*radius[:, None]*weights[None, :]/2, axis=1)
        return abs(np.trapz(radial_integral, chi)/(2*np.pi))


def build_chart(gfile, table_path, machine_path, *, r_min, r_max, nr=17, ntheta=128,
                native_q_sign, max_relative_q_error, max_relative_flux_error,
                max_relative_label_flux_error):
    """Preserve recovered r_eff and outboard phase; verify geometry independently."""
    for tolerance in (max_relative_q_error, max_relative_flux_error, max_relative_label_flux_error):
        if not np.isfinite(tolerance) or tolerance < 0:
            raise ValueError("finite nonnegative independent geometry tolerances required")
    table = np.loadtxt(table_path)
    machine = np.loadtxt(machine_path).ravel()
    if table.ndim != 2 or table.shape[1] != 11 or machine.size != 2:
        raise ValueError("expected 11-column equilibrium table and Btf/R0 machine pair")
    Btf, R0 = machine
    if not np.all(np.isfinite(table)) or not np.isfinite(Btf) or Btf == 0 or not np.isfinite(R0) or R0 <= 0:
        raise ValueError("invalid recovered equilibrium metadata")
    source_r = table[:, 0]
    if np.any(np.diff(source_r) <= 0) or np.any(source_r <= 0):
        raise ValueError("native toroidal-flux radius must be strictly increasing")
    reconstructed_r = np.sqrt(2*abs(table[:, 3]/Btf))
    if not np.allclose(source_r, reconstructed_r, rtol=1e-12, atol=1e-12):
        raise ValueError("equilibrium table radius violates its toroidal-flux definition")
    if not source_r[0] <= r_min < r_max <= source_r[-1] or nr < 4 or ntheta < 16:
        raise ValueError("chart needs >=4 radii and >=16 angles inside the recovered flux domain")
    if native_q_sign not in (-1, 1):
        raise ValueError("native q sign relative to the recovered table must be explicit")
    eq = Equilibrium(read_geqdsk(gfile))
    r, theta = np.linspace(r_min, r_max, nr), np.arange(ntheta)*2*np.pi/ntheta
    R, Z = np.empty((nr, ntheta)), np.empty((nr, ntheta))
    q_geometry, flux_geometry, label_flux_geometry = np.empty(nr), np.empty(nr), np.empty(nr)
    q_table = native_q_sign*PchipInterpolator(source_r, table[:, 1])(r)
    source_flux = PchipInterpolator(source_r, table[:, 3])(r)
    R_start = PchipInterpolator(source_r, table[:, 7])(r)
    Z_start = PchipInterpolator(source_r, table[:, 8])(r)
    for j in range(nr):
        target_psi = eq.psi(R_start[j], Z_start[j])[0, 0]
        th, Rc, Zc, q_geometry[j], chi, radius = eq.contour(
            target_psi, [R_start[j], Z_start[j]], max(256, 4*ntheta))
        Rc[-1], Zc[-1] = Rc[0], Zc[0]
        R[j] = CubicSpline(th, Rc, bc_type="periodic")(theta)
        Z[j] = CubicSpline(th, Zc, bc_type="periodic")(theta)
        flux_geometry[j] = eq.toroidal_flux_per_radian(chi, radius)
        # The legacy KAMEL accumulator uses F(surface) for the whole enclosed
        # weighted area. Preserve that label and report its finite-beta error.
        label_flux_geometry[j] = eq.toroidal_flux_per_radian(chi, radius,
                                                           boundary_F=eq.F(target_psi))
    relative_q_error = float(np.max(abs(q_geometry-q_table)/abs(q_table)))
    relative_flux_error = float(np.max(abs(flux_geometry-source_flux)/abs(source_flux)))
    relative_label_error = float(np.max(abs(label_flux_geometry-source_flux)/abs(source_flux)))
    diagnostics = dict(r_cm=r.tolist(), q_geometry=q_geometry.tolist(), q_native_table=q_table.tolist(),
                       relative_q_error=relative_q_error, relative_toroidal_flux_error=relative_flux_error,
                       relative_native_label_flux_error=relative_label_error,
                       native_label_flux_Gcm2=source_flux.tolist(),
                       reconstructed_native_label_flux_Gcm2=label_flux_geometry.tolist(),
                       true_toroidal_flux_Gcm2=flux_geometry.tolist())
    if (relative_q_error > max_relative_q_error or relative_flux_error > max_relative_flux_error
            or relative_label_error > max_relative_label_flux_error):
        raise ValueError(f"matched GEQDSK fails recovered native chart checks: {diagnostics}")
    def theta_derivative(values):
        # Use the same coordinate interpolant as SurfaceChart. Spectral
        # derivatives here and cubic derivatives in continuous evaluation
        # would define different one-form samples between the same nodes.
        angular = CubicSpline(np.append(theta, 2*np.pi),
                              np.concatenate((values, values[:, :1]), axis=1),
                              axis=1, bc_type="periodic")
        return angular(theta, 1)
    geometry = dict(r=np.broadcast_to(r[:, None], R.shape),
                    theta=np.broadcast_to(theta[None, :], R.shape), R=R, Z=Z,
                    nu=np.zeros_like(R), R_r=CubicSpline(r, R, axis=0)(r, 1),
                    Z_r=CubicSpline(r, Z, axis=0)(r, 1), R_theta=theta_derivative(R),
                    Z_theta=theta_derivative(Z), nu_r=np.zeros_like(R), nu_theta=np.zeros_like(R))
    _, J = deformation(geometry["r"], R, R0, *(geometry[key] for key in
                        ("R_r", "R_theta", "Z_r", "Z_theta", "nu_r", "nu_theta")))
    diagnostics["min_oriented_J"] = float(J.min())
    provenance = dict(method="independent clockwise straight-field-line GEQDSK reconstruction",
                      native_angle="clockwise from recovered outboard phase anchor",
                      toroidal_angle="geometric phi=z/R0; nu=0",
                      native_q_sign_relative_to_table=native_q_sign,
                      native_flux_label="legacy boundary-F contour integral; finite-beta approximation",
                      historical_angle_table_recovered=False,
                      GEQDSK_sha256=digest(gfile), equilibrium_table_sha256=digest(table_path),
                      machine_data_sha256=digest(machine_path), diagnostics=diagnostics)
    geometry.update(radial_convention=RADIAL_CONVENTION, phase_convention=PHASE_CONVENTION,
                    length_units="cm", R0_cm=R0, Btf_G=Btf,
                    provenance=json.dumps(provenance, sort_keys=True))
    return geometry, diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geqdsk", type=Path)
    parser.add_argument("equilibrium_table", type=Path)
    parser.add_argument("machine_data", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--r-min-cm", type=float, required=True)
    parser.add_argument("--r-max-cm", type=float, required=True)
    parser.add_argument("--nr", type=int, default=17)
    parser.add_argument("--ntheta", type=int, default=128)
    parser.add_argument("--native-q-sign", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--max-relative-q-error", type=float, required=True)
    parser.add_argument("--max-relative-flux-error", type=float, required=True)
    parser.add_argument("--max-relative-label-flux-error", type=float, required=True)
    args = parser.parse_args()
    geometry, diagnostics = build_chart(args.geqdsk, args.equilibrium_table, args.machine_data,
        r_min=args.r_min_cm, r_max=args.r_max_cm, nr=args.nr, ntheta=args.ntheta,
        native_q_sign=args.native_q_sign, max_relative_q_error=args.max_relative_q_error,
        max_relative_flux_error=args.max_relative_flux_error,
        max_relative_label_flux_error=args.max_relative_label_flux_error)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as stream:
        np.savez_compressed(stream, **geometry)
    print(json.dumps(diagnostics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
