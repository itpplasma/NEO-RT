"""Reconstruct native torque from root rows, including concurrent append calls."""
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

executable, fixture, mode = sys.argv[1:]
with tempfile.TemporaryDirectory(prefix="resonance-ledger-") as directory:
    work = Path(directory)
    for name in ("in_file", "plasma.in", "profile.in"):
        (work / name).symlink_to(Path(fixture) / name)
    ledger = work / "roots.txt"
    expected_path = work / "expected.txt"
    env = dict(os.environ, OMP_NUM_THREADS="2", OMP_THREAD_LIMIT="2",
               OPENBLAS_NUM_THREADS="1")
    env.pop("NEORT_RESONANCE_LEDGER", None)

    def run():
        return subprocess.run([executable, mode, str(expected_path)], cwd=work,
                              env=env, capture_output=True, text=True, timeout=90)

    if mode in ("krook", "no-torque"):
        env["NEORT_RESONANCE_LEDGER"] = str(ledger)
        result = run()
        message = ("requires collisionless transport" if mode == "krook"
                   else "requires comptorque=.true.")
        assert result.returncode != 0 and message in result.stdout + result.stderr, result
    else:
        baseline = run()
        assert baseline.returncode == 0, baseline.stdout + baseline.stderr
        expected = expected_path.read_text()
        assert not ledger.exists()
        env["NEORT_RESONANCE_LEDGER"] = str(ledger)
        # Repeat appends exercise both initial header creation and reopening.
        for repeat in range(2):
            result = run()
            assert result.returncode == 0, result.stdout + result.stderr
            assert expected_path.read_text() == expected, "ledger changes native result"
        lines = ledger.read_text().splitlines()
        headers = [line for line in lines if line.startswith("#")]
        assert len(headers) == 1, headers
        columns = headers[0].removeprefix("# ").split()
        rows = [dict(zip(columns, line.split(), strict=True))
                for line in lines if not line.startswith("#")]
        assert rows, "ledger must contain resonance roots"
        sums = {}
        counts = {}
        for row in rows:
            key = (round(float(row["s"]), 12), int(row["species"]), int(row["mth"]),
                   row["branch"])
            contribution = float(row["dTds"])
            assert math.isfinite(contribution), row
            assert int(row["istate"]) == 2, row
            if mode == "mono":
                assert float(row["sample_weight"]) == 1.0, row
                assert float(row["ux"]) == 1.5, row
            sums[key] = sums.get(key, 0.0) + contribution
            counts[key[:2]] = counts.get(key[:2], 0) + 1
        for line in expected.splitlines():
            s, species, mth, *torques = line.split()
            prefix = (round(float(s), 12), int(species))
            assert counts.get(prefix, 0) > 0, ("missing surface/species", prefix)
            for branch, value in zip(("cop", "ctr", "trap"), torques, strict=True):
                key = (*prefix, int(mth), branch)
                assert math.isclose(sums.pop(key, 0.0), 2*float(value),
                                    rel_tol=2e-13, abs_tol=1e-280), (key, value)
        assert any(float(line.split()[-1]) != 0.0 or
                   float(line.split()[-2]) != 0.0 or
                   float(line.split()[-3]) != 0.0 for line in expected.splitlines())
        assert not sums, ("unexpected attribution", sums)
