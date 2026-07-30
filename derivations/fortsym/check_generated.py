#!/usr/bin/env python3
"""Regenerate the Boozer kernels and compare them byte-for-byte."""

from pathlib import Path
import subprocess
import tempfile


HERE = Path(__file__).resolve().parent
KERNELS = {
    "generate_boozer_ampere": (
        HERE.parents[1] / "src/generated/boozer_ampere_kernel.f90"
    ),
    "generate_boozer_jxb": (
        HERE.parents[1] / "src/generated/boozer_jxb_kernel.f90"
    ),
}


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="neo_rt_fortsym_") as directory:
        for executable, expected in KERNELS.items():
            candidate = Path(directory) / expected.name
            subprocess.run(
                [
                    "fpm",
                    "run",
                    "--profile",
                    "release",
                    executable,
                    "--",
                    str(candidate),
                ],
                cwd=HERE,
                check=True,
            )
            if candidate.read_bytes() != expected.read_bytes():
                raise SystemExit(
                    f"generated {expected.name} is stale; regenerate it with "
                    f"'cd derivations/fortsym && fpm run --profile release "
                    f"{executable}'"
                )
    print("generated Boozer kernels: OK")


if __name__ == "__main__":
    main()
