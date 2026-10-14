"""Energy ladder: isolate phase, radial amplitude, and full-orbit deposition.

Decision thresholds are declared in registry_proposal.json before execution.
Incomplete or invalid data never produce a physics verdict.
"""
import json
import sys
from pathlib import Path

import numpy as np

HALVES = (.02, .04, .06)
MODES = ("thin", "dph", "radial", "both")
R0, A = 6.2, .62


def volume(rho):
    c = R0 - np.sqrt(R0**2 - A**2)
    return 2 * np.pi**2 * R0 * (R0**2 - (R0 - np.asarray(rho)**2 * c)**2)


def potato_window(path, rho, h):
    tb, lo, hi = np.loadtxt(path / "boxcounted_torque.dat", unpack=True)
    lo, hi = np.sqrt(lo), np.sqrt(hi)
    weights = np.maximum(0, np.minimum(hi, rho + h) - np.maximum(lo, rho - h)) / (hi - lo)
    return float(np.sum(tb * weights) * 1e-7 /
                 np.sum(weights * (volume(hi) - volume(lo))))


def neort_torque(path):
    row = np.loadtxt(path / "driftorbit_torque.out", ndmin=2)[-1]
    return float(.1 * sum(row[3:6]) / row[1])


def window(rho_grid, torque, rho, h):
    grid = np.linspace(rho - h, rho + h, 4001)
    values = np.exp(np.interp(grid, rho_grid, np.log(np.abs(torque))))
    dv = np.gradient(volume(grid))
    return float(np.sum(values * dv) / np.sum(dv) * np.sign(torque[0]))


def deposited(rho_grid, torque, rho, h, width, kernel):
    grid = np.linspace(max(rho_grid[0], rho - h - width),
                       min(rho_grid[-1], rho + h + width), 20001)
    if grid[0] > rho - h - width or grid[-1] < rho + h + width:
        raise ValueError("deposition kernel exceeds the declared NEO-RT grid")
    values = np.exp(np.interp(grid, rho_grid, np.log(np.abs(torque))))
    cdf = (lambda x: np.clip((x + 1) / 2, 0, 1)) if kernel == "tophat" else (
        lambda x: .5 + np.arcsin(np.clip(x, -1, 1)) / np.pi)
    fraction = cdf((rho + h - grid) / width) - cdf((rho - h - grid) / width)
    return float(np.sum(values * fraction * np.gradient(volume(grid))) /
                 (volume(rho + h) - volume(rho - h)) * np.sign(torque[0]))


def main(root):
    root = Path(root)
    scan = json.loads((root / "inputs/scan.json").read_text())
    rows, failures, profiles = {}, [], {}
    for c in scan["cases"]:
        cid, rho = c["id"], np.sqrt(c["s"])
        try:
            pd = root / "out" / ("p" + cid) / "run"
            boxes = np.loadtxt(pd / "boxcounted_torque.dat", ndmin=2)
            integral = float(np.loadtxt(pd / "integral_torque.dat"))
            closure = float(sum(boxes[:, 0]) / integral)
            row = {"energy_eV": c["T0_eV"], "s": c["s"],
                   "closure": closure, "windows": {}}
            if not np.isfinite(closure) or abs(closure - 1) > 1e-6:
                failures.append(f"{cid}: box/integral closure fails")
            if c["kind"] != "energy_resolution":
                for mode in MODES:
                    values = np.array([neort_torque(root / "out/neort" /
                        f"{cid}_{r:.2f}_{mode}") for r in c["rho_grid"]])
                    if not np.all(np.isfinite(values)) or np.any(values == 0):
                        raise ValueError("nonfinite or zero NEO-RT torque")
                    if not np.all(np.sign(values) == np.sign(values[0])):
                        raise ValueError("NEO-RT torque changes sign in analysis grid")
                    profiles[cid, mode] = values
            for h in HALVES:
                p = potato_window(pd, rho, h)
                w = {"P": p}
                if c["kind"] != "energy_resolution":
                    for mode in MODES:
                        w[mode] = window(c["rho_grid"], profiles[cid, mode], rho, h)
                        if np.sign(p) != np.sign(w[mode]):
                            failures.append(f"{cid}: signs disagree ({mode})")
                    w.update(phase_change=w["dph"] / w["thin"] - 1,
                             radial_change=w["both"] / w["dph"] - 1,
                             residual=p / w["both"] - 1)
                    if c["kind"] == "ladder":
                        factor = np.sqrt(c["T0_eV"] / 2000)
                        delta = .0348 if c["s"] == .09 else .0246
                        for kernel, coefficient in (("tophat", 2.64), ("arcsine", 2.19)):
                            k = deposited(c["rho_grid"], profiles[cid, "both"],
                                          rho, h, coefficient * delta * factor, kernel)
                            w[kernel] = k
                            w[kernel + "_residual"] = p / k - 1
                row["windows"][str(h)] = w
            rows[cid] = row
        except (OSError, ValueError, IndexError, ZeroDivisionError) as exc:
            failures.append(f"{cid}: {exc}")
    gates = {}
    for c in scan["cases"]:
        if c["id"] not in rows:
            continue
        row = rows[c["id"]]
        if c["kind"] in ("half_amplitude", "energy_resolution"):
            base = rows.get(c["base"])
            if base is None:
                failures.append(c["id"] + ": missing base control")
                continue
            multiplier = 4 if c["kind"] == "half_amplitude" else 1
            ratio = multiplier * row["windows"]["0.06"]["P"] / base["windows"]["0.06"]["P"]
            gates[c["id"] + "_P"] = {"ratio": ratio, "pass": abs(ratio - 1) <= .05}
            if c["kind"] == "half_amplitude":
                for mode in MODES:
                    ratio = 4 * row["windows"]["0.06"][mode] / base["windows"]["0.06"][mode]
                    gates[c["id"] + "_" + mode] = {"ratio": ratio, "pass": abs(ratio - 1) <= .05}
        elif c["T0_eV"] == 20:
            ratio = row["windows"]["0.06"]["P"] / row["windows"]["0.06"]["thin"]
            gates[c["id"] + "_cold"] = {"ratio": ratio, "pass": abs(ratio - 1) <= .05}
    for name, gate in gates.items():
        if not gate["pass"]:
            failures.append(name + ": validity gate fails")
    diagnostics = []
    for c in scan["cases"]:
        if c["kind"] == "ladder" and c["id"] in rows:
            ws = rows[c["id"]]["windows"]
            w = ws["0.04"]
            sigma = max(.03, abs(w["residual"] - ws["0.06"]["residual"]))
            diagnostics.append({"id": c["id"], "phase": w["phase_change"],
                "radial": w["radial_change"], "residual": w["residual"], "sigma": sigma,
                "phase_resolved": abs(w["phase_change"]) > 2 * sigma,
                "deposition_predicted": any(abs(w[k + "_residual"]) <= max(2 * sigma, .25)
                                             for k in ("tophat", "arcsine"))})
    result = {"status": "invalid/incomplete" if failures else "valid ladder",
              "failures": failures, "gates": gates, "diagnostics": diagnostics, "cases": rows}
    (root / "out").mkdir(exist_ok=True)
    (root / "out/analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"])
    for d in diagnostics:
        print(f"{d['id']}: phase {d['phase']:+.4f}; radial {d['radial']:+.4f}; "
              f"P/Nboth-1 {d['residual']:+.4f} +/- {d['sigma']:.4f}; "
              f"deposition prediction {d['deposition_predicted']}")
    for failure in failures:
        print(failure)
    return int(bool(failures))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
