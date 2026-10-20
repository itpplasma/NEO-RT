"""Run the retained MARS node in an isolated directory with its physical inputs."""
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

executable, fixture = sys.argv[1:]
with tempfile.TemporaryDirectory(prefix="krook-reference-") as directory:
    with tarfile.open(fixture) as archive:
        for name in ("in_file", "in_file_pert", "r0p50000.in",
                     "r0p50000.species.in"):
            with archive.extractfile(name) as stream:
                (Path(directory) / name).write_bytes(stream.read())
    subprocess.run([executable], cwd=directory, check=True, timeout=60)
