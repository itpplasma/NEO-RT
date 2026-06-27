#!/usr/bin/env python3
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SKIP = 77
DEFAULT_CASE = Path("potato_benchmarks/rung4_torque_30835/run/potato_zero_mid_e12_n16")
INPUT_FILES = (
    "efit.eqdsk",
    "eqmagprofs.dat",
    "potato.in",
    "profile_poly.in",
    "field_divB0.inp",
    "bmod_n.dat",
    "convexwall.dat",
)
PROBE_FIELDS = (
    "  probe_rho_pol = 0.9d0",
    "  probe_ux = 1.378773740d0",
    "  probe_eta = 4.119637570d-5",
    "  probe_m = 0",
    "  probe_n = 2",
)


def stage_case(case_dir, dst_dir, extra_lines):
    for name in INPUT_FILES:
        shutil.copy(case_dir / name, dst_dir / name)
    potato_in = dst_dir / "potato.in"
    lines = potato_in.read_text().splitlines()
    out = []
    inserted = False
    for line in lines:
        if not inserted and line.strip() == "/":
            out.extend(extra_lines)
            inserted = True
        out.append(line)
    if not inserted:
        raise SystemExit("no namelist terminator '/' in potato.in")
    potato_in.write_text("\n".join(out) + "\n")


def run_probe(exe, work):
    subprocess.run([str(exe)], cwd=work, check=True)
    return (work / "potato_resonance_probe.dat").read_text()


def bracket_classes(text):
    return sorted(
        {int(m.group(1)) for m in re.finditer(r"^# bracket\s+(\d+)", text, re.MULTILINE)}
    )


def main():
    if len(sys.argv) != 2:
        print("usage: test_probe_clip.py <potato_resonance_probe.x>", file=sys.stderr)
        return 2

    exe = Path(sys.argv[1])
    case_dir = Path(os.environ.get("POTATO_PROBE_CASE", DEFAULT_CASE))
    if not exe.is_file() or not case_dir.is_dir():
        return SKIP
    if not all((case_dir / name).is_file() for name in INPUT_FILES):
        return SKIP

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        clipped = tmp / "clipped"
        clipped.mkdir()
        stage_case(case_dir, clipped, PROBE_FIELDS)
        clipped_classes = bracket_classes(run_probe(exe, clipped))
        if clipped_classes:
            print(
                f"default clip must report no brackets, got classes {clipped_classes}",
                file=sys.stderr,
            )
            return 1

        full = tmp / "full"
        full.mkdir()
        stage_case(case_dir, full, ("  clip_resonance_classes = .false.", *PROBE_FIELDS))
        full_classes = bracket_classes(run_probe(exe, full))
        if not {3, 5}.issubset(full_classes):
            print(
                "full domain must report m=0 brackets in classes 3 and 5, "
                f"got classes {full_classes}",
                file=sys.stderr,
            )
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
