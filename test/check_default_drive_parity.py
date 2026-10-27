"""Compare default-off outputs to a separately built main on this machine."""
import subprocess
import sys
import tempfile
from pathlib import Path

import f90nml

candidate, reference, inputs = map(Path, sys.argv[1:])
cases = [("s01", 0.1, False), ("s05", 0.5, False), ("s09", 0.9, False),
         ("krook", 0.5, True)]
with tempfile.TemporaryDirectory(prefix="default-drive-", dir=Path.cwd()) as tmp:
    for label, surface, krook in cases:
        outputs = []
        for name, exe in [("main", reference), ("pr", candidate)]:
            work = Path(tmp) / label / name
            work.mkdir(parents=True)
            for filename in ["in_file", "in_file_pert", "plasma.in", "profile.in"]:
                (work / filename).symlink_to(inputs / filename)
            config = f90nml.read(inputs / "template.in")
            config["params"]["s"] = surface
            if krook:
                config["params"].update(
                    krook_collisions=True, krook_nu=1e5, krook_model=0,
                    krook_eta_steps=16, krook_max_depth=12, krook_rtol=1e-3,
                    vsteps=8, mth_max_abs=1,
                )
            config.write(work / "case.in")
            result = subprocess.run([str(exe), "case"], cwd=work, capture_output=True,
                                    text=True, timeout=600)
            assert result.returncode == 0, result.stdout + result.stderr
            outputs.append({p.name: p.read_bytes() for p in work.glob("*.out")})
        assert outputs[0], f"no reference outputs for {label}"
        assert outputs[0].keys() == outputs[1].keys(), f"output file set changed: {label}"
        for filename in outputs[0]:
            assert outputs[0][filename] == outputs[1][filename], (label, filename)
        print(label, "byte-identical", len(outputs[0]), "output files", flush=True)
