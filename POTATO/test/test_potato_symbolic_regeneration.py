#!/usr/bin/env python3
"""Regenerate POTATO's Fortsym kernels and compare bytes with the tree."""

from __future__ import annotations

import subprocess
import sys
import tempfile
import os
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_potato_symbolic_regeneration.py REPOSITORY", file=sys.stderr)
        return 2

    repository = Path(sys.argv[1]).resolve()
    generator_dir = repository / "tools" / "gc_symbolics"
    expected_dir = repository / "POTATO" / "SRC" / "generated"
    expected = sorted(expected_dir.glob("potato_*.f90"))
    if not expected:
        print("no committed POTATO Fortsym outputs found", file=sys.stderr)
        return 1

    # Keep generated artifacts on the same storage policy as faep tests.  The
    # default is beside the checkout, never Python's implicit /tmp directory;
    # CI may select /var/tmp/ert explicitly with NEORT_TEST_TMPDIR.
    scratch_root = Path(
        os.environ.get("NEORT_TEST_TMPDIR") or repository.parent / ".neort-test-tmp"
    ).resolve()
    tmp_root = Path("/tmp").resolve()
    if scratch_root == tmp_root or tmp_root in scratch_root.parents:
        print("NEORT_TEST_TMPDIR must not point under /tmp", file=sys.stderr)
        return 2
    scratch_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="potato-fortsym-", dir=scratch_root) as temporary:
        generated_dir = Path(temporary)
        result = subprocess.run(
            ["fo", "exec", "gen_potato_kernels", str(generated_dir)],
            cwd=generator_dir,
            check=False,
        )
        if result.returncode != 0:
            return result.returncode

        generated = sorted(generated_dir.glob("potato_*.f90"))
        if [path.name for path in generated] != [path.name for path in expected]:
            print("Fortsym output file set differs", file=sys.stderr)
            print("expected:", [path.name for path in expected], file=sys.stderr)
            print("generated:", [path.name for path in generated], file=sys.stderr)
            return 1

        for expected_path, generated_path in zip(expected, generated):
            if expected_path.read_bytes() != generated_path.read_bytes():
                print(
                    f"Fortsym regeneration differs: {expected_path}",
                    file=sys.stderr,
                )
                return 1

    print("POTATO Fortsym regeneration matches committed bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
