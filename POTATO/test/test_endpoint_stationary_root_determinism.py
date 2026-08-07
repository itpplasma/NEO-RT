#!/usr/bin/env python3
"""Run the focused topology oracle under serial and threaded OpenMP settings."""

from __future__ import annotations

import os
import subprocess
import sys


def run(executable: str, threads: str) -> bytes:
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": threads,
            "OMP_DYNAMIC": "FALSE",
            "OMP_PROC_BIND": "close",
            "OMP_PLACES": "cores",
        }
    )
    completed = subprocess.run(
        [executable],
        check=False,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode != 0:
        sys.stderr.buffer.write(completed.stdout)
        raise SystemExit(
            f"endpoint oracle failed under OMP_NUM_THREADS={threads}: "
            f"exit {completed.returncode}"
        )
    return completed.stdout


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} ENDPOINT_EXECUTABLE", file=sys.stderr)
        return 2
    serial = run(sys.argv[1], "1")
    threaded_first = run(sys.argv[1], "16")
    threaded_second = run(sys.argv[1], "16")
    if serial != threaded_first or threaded_first != threaded_second:
        print("endpoint oracle output is not deterministic across OMP settings", file=sys.stderr)
        return 1
    print("endpoint oracle OMP=1/16/16 output is deterministic")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
