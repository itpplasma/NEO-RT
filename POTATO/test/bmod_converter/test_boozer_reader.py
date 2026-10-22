"""An analytic synthetic chart and BC spectrum; no experiment data."""
import importlib.util
from pathlib import Path
import subprocess
import sys

import numpy as np


def main():
    binary, tool, work = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    precision = sys.argv[4:] == ["precision"]
    work.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("geometry", tool)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    s, theta = np.linspace(0, 1, 33), np.arange(64)*2*np.pi/64
    ss, tt = np.meshgrid(s, theta, indexing="ij")
    geometry = np.stack((300 + 100*ss*np.cos(tt), 100*ss*np.sin(tt),
                         -0.13 - 0.02*np.sin(tt)), axis=-1)
    modes = np.array([0, 1])
    inner = module.support_disk(modes, np.array([[300, 0], [1, -1j]]), True)
    outer = module.support_disk(modes, np.array([[300, 0], [99, -99j]]), False)
    assert 0 < inner[2] < 1 and outer[2] > 99**2
    module.write_chart(work/"geometry.nc", s, theta, geometry, [0.01, 0.99],
                       np.r_[inner, outer])
    if precision:
        s = np.linspace(0, 1, 4001)
        ss, tt = np.meshgrid(s, theta, indexing="ij")
        radius = 200*np.sqrt(ss)
        geometry = np.stack((650 + radius*np.cos(tt), radius*np.sin(tt),
                             -0.13 - 0.02*np.sin(tt)), axis=-1)
        module.write_chart(work/"geometry.nc", s, theta, geometry,
                           [0.00025, 0.99975], [650, 0, 1, 650, 0, 201**2])
    for n in (3, -3):
        lines = ["CC synthetic analytic spectrum"]*4
        lines += ["m n ns nfp flux a R", "2 0 9 1 0 0 0"]
        for sj in np.linspace(0.1, 0.9, 9):
            lines += ["s iota Jpol Itor pprime sqrtg", "units",
                      f"{sj:.17e} -1 0 0 0 0", "m n geometry bmnc bmns"]
            for m in (-1, 0, 1):
                lines += [f"{m} {n} 0 0 0 0 0 0 "
                          f"{(m+2)*(1+sj+sj*sj)*1e-4:.17e} "
                          f"{(m-2)*(1-sj)*1e-4:.17e}"]
        (work/"pert.bc").write_text("\n".join(lines) + "\n")
        (work/"potato.in").write_text(
            f"&potato_nml n_tor={n}, use_boozer_perturbation=.true., "
            "boozer_geometry_file='geometry.nc', boozer_perturbation_file='pert.bc' /\n")
        subprocess.run([str(binary)], cwd=work, check=True)


if __name__ == "__main__":
    main()
