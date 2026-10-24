#!/usr/bin/env python3
"""Regenerate a displacement golden record that names committed public inputs.

usage: gen_golden.py <potato.x> <golden.json> [threads]

The record's "inputs" and "files" entries are kept. The case is staged exactly
as test_potato_golden.py stages it, potato.x runs once, and the script rewrites
"integral_torque" and the rows file (every row of every listed output file).
Run it with a Release build and review the note fields by hand afterwards.
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from test_potato_golden import stage_inputs


def main():
    if len(sys.argv) not in (3, 4):
        print(__doc__, file=sys.stderr)
        return 2
    exe = Path(sys.argv[1]).resolve()
    golden_path = Path(sys.argv[2]).resolve()
    threads = sys.argv[3] if len(sys.argv) == 4 else "8"
    golden = json.loads(golden_path.read_text())
    with tempfile.TemporaryDirectory(prefix="potato-golden-gen-") as tmp:
        work = Path(tmp) / "case"
        stage_inputs(golden_path, golden["inputs"], work)
        subprocess.run([str(exe)], cwd=work, check=True,
                       env=dict(os.environ, OMP_NUM_THREADS=threads),
                       stdout=subprocess.DEVNULL)
        rows = {}
        for name, spec in golden["files"].items():
            data = np.loadtxt(work / name, ndmin=2)
            if not np.all(np.isfinite(data)) or data.size == 0:
                raise SystemExit(f"{name}: empty or nonfinite output")
            ncol = data.shape[1]
            columns = spec["key_columns"] + sum(spec.get("scale_groups", []), [])
            if max(columns) >= ncol:
                raise SystemExit(f"{name}: column spec exceeds {ncol} columns")
            rows[name] = data
        golden["integral_torque"] = float(np.loadtxt(work / "integral_torque.dat"))
    np.savez_compressed(golden_path.parent / golden["rows_file"], **rows)
    golden_path.write_text(json.dumps(golden, indent=2) + "\n")
    for name, data in rows.items():
        print(f"{name}: {data.shape}")
    print(f"integral_torque = {golden['integral_torque']!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
