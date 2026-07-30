#!/usr/bin/env python3
"""Regenerate the Boozer Ampere kernel and compare it byte-for-byte."""

from pathlib import Path
import subprocess
import tempfile


HERE = Path(__file__).resolve().parent
EXPECTED = HERE.parents[1] / "src/generated/boozer_ampere_kernel.f90"


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="neo_rt_fortsym_") as directory:
        candidate = Path(directory) / EXPECTED.name
        subprocess.run(
            ["fpm", "run", "--profile", "release", "--", str(candidate)],
            cwd=HERE,
            check=True,
        )
        if candidate.read_bytes() != EXPECTED.read_bytes():
            raise SystemExit(
                "generated Boozer Ampere kernel is stale; regenerate it with "
                "'cd derivations/fortsym && fpm run --profile release'"
            )
    print("generated Boozer Ampere kernel: OK")


if __name__ == "__main__":
    main()
