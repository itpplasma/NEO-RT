"""Require a deliberate diagnostic, rather than a reader or bounds failure."""
import subprocess
import sys

messages = {
    "missing-argument": "pertfile_curvature = .true. needs a curvature file",
    "missing-file": "required if pertfile_curvature = .true.",
    "without-pertfile": "pertfile_curvature requires pertfile = .true.",
}
result = subprocess.run(sys.argv[1:], text=True, capture_output=True, timeout=20)
assert result.returncode != 0, "invalid curvature input accepted"
assert messages[sys.argv[2]] in result.stdout + result.stderr, result.stdout + result.stderr
