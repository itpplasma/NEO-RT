"""Extract the circular-field inputs and run the trapped high-harmonic regression."""
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

executable, fixture = sys.argv[1:]
with tempfile.TemporaryDirectory(prefix="krook-trapped-") as directory:
    with tarfile.open(fixture) as archive:
        for name in ("driftorbit.in", "plasma.in", "profile.in", "in_file",
                     "in_file_pert"):
            with archive.extractfile(name) as stream:
                (Path(directory) / name).write_bytes(stream.read())
    subprocess.run([executable], cwd=directory, check=True, timeout=120)
