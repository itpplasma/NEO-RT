"""The angle-map loader must reject the opposite stream-function convention.

Run through CTest with its normal angle-map executable and fixture environment.
The fixture geometry and direct field stay fixed; only the independent Boozer
stream-function column changes sign.
"""

import os
from pathlib import Path
import subprocess
import sys
import tempfile


def opposite_stream_function(source, destination):
    lines = Path(source).read_text().splitlines()
    for i, line in enumerate(lines):
        row = line.split()
        if len(row) != 10:
            continue
        try:
            values = [float(value.replace("D", "E")) for value in row]
        except ValueError:
            continue
        values[6] = -values[6]
        values[7] = -values[7]
        lines[i] = f"{int(values[0]):5d}{int(values[1]):5d}" + "".join(
            f"{value:17.8e}" for value in values[2:]
        )
    Path(destination).write_text("\n".join(lines) + "\n")


with tempfile.TemporaryDirectory(prefix="neort-stream-convention-") as work:
    fixture = Path(work) / "opposite_stream.bc"
    opposite_stream_function(os.environ["BOOZER_SOLOVEV_FILE"], fixture)
    env = dict(os.environ, BOOZER_SOLOVEV_FILE=str(fixture))
    result = subprocess.run(
        [sys.argv[1]], env=env, text=True, capture_output=True, timeout=120
    )
    diagnostic = "inconsistent Boozer stream-function convention"
    if result.returncode == 0 or diagnostic not in result.stdout + result.stderr:
        print((result.stdout + result.stderr)[-3000:])
        raise SystemExit("loader did not reject the opposite stream-function convention")
    print("PASS opposite stream-function convention rejected at load")
