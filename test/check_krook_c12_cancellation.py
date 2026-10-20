"""Extract the exact C12 inputs and run the cancelling-channel regression."""
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

executable, fixture = sys.argv[1:]
with tempfile.TemporaryDirectory(prefix="krook-c12-") as directory:
    with tarfile.open(fixture) as archive:
        for name in ("driftorbit.in", "species.in", "plasma.in", "profile.in",
                     "in_file", "in_file_pert"):
            with archive.extractfile(name) as stream:
                (Path(directory) / name).write_bytes(stream.read())
    subprocess.run([executable], cwd=directory, check=True, timeout=60)
