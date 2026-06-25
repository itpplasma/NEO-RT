"""Generate a Boozer chartmap NetCDF from circ.bc for test_chartmap_input.

Writes the chartmap to the path given as the first command-line argument.

The script loads booz_xform_to_boozer_chartmap and bc_to_booz_xform from the
path given by the NEORT_LIBNEO_CONVERTERS env var (a directory containing the
libneo Python subpackage with those two modules), or from any directory on
PYTHONPATH that has a libneo/bc_to_booz_xform.py under it.  The main libneo
installed in the repo does not yet include bc_to_booz_xform; it lives in the
eqdsk-booz worktree.  We load it directly by file path to avoid name-shadowing.

Usage (CMake fixture injects env):
    NEORT_LIBNEO_CONVERTERS=/path/to/eqdsk-booz/python \\
        python3 test/gen_circ_chartmap.py <output.nc>
"""

import importlib.util
import os
import sys
import tempfile

if len(sys.argv) != 2:
    print("usage: gen_circ_chartmap.py <output.nc>", file=sys.stderr)
    sys.exit(1)

output_nc = sys.argv[1]

# circ.bc lives in ../../examples/circ.bc relative to this script.
script_dir = os.path.dirname(os.path.abspath(__file__))
bc_file = os.path.join(script_dir, "..", "examples", "circ.bc")
if not os.path.exists(bc_file):
    print(f"ERROR: circ.bc not found at {bc_file}", file=sys.stderr)
    sys.exit(1)


def _load_module_from_file(name, filepath):
    """Import a module directly from a .py path, bypassing package shadowing."""
    spec = importlib.util.spec_from_file_location(name, filepath)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# Locate the directory that contains libneo/bc_to_booz_xform.py.
# Search: NEORT_LIBNEO_CONVERTERS env, then PYTHONPATH entries.
def _find_converter_dir():
    candidates = []
    env_dir = os.environ.get("NEORT_LIBNEO_CONVERTERS", "")
    if env_dir:
        candidates.append(env_dir)
    pythonpath = os.environ.get("PYTHONPATH", "")
    for p in pythonpath.split(os.pathsep):
        if p:
            candidates.append(p)

    for d in candidates:
        candidate = os.path.join(d, "libneo", "bc_to_booz_xform.py")
        if os.path.isfile(candidate):
            return d
    return None


converter_dir = _find_converter_dir()
if converter_dir is None:
    print("ERROR: cannot locate libneo/bc_to_booz_xform.py.", file=sys.stderr)
    print("Set NEORT_LIBNEO_CONVERTERS to the directory containing the libneo/ "
          "subpackage (e.g. the eqdsk-booz worktree python/ dir).", file=sys.stderr)
    sys.exit(1)

# Load the two converter modules directly by file path.
_load_module_from_file(
    "libneo.bc_to_booz_xform",
    os.path.join(converter_dir, "libneo", "bc_to_booz_xform.py"),
)
_load_module_from_file(
    "libneo.booz_xform_to_boozer_chartmap",
    os.path.join(converter_dir, "libneo", "booz_xform_to_boozer_chartmap.py"),
)
_load_module_from_file(
    "libneo.boozer_chartmap_writer",
    os.path.join(converter_dir, "libneo", "boozer_chartmap_writer.py"),
)

from libneo.bc_to_booz_xform import convert_bc_to_boozmn  # noqa: E402
from libneo.booz_xform_to_boozer_chartmap import convert_boozmn_to_chartmap  # noqa: E402

with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
    boozmn_nc = tmp.name

try:
    convert_bc_to_boozmn(bc_file, boozmn_nc)
    # nzeta>=2 required by boozer_chartmap_io (endpoint-excluded grid needs at least
    # 2 points to be valid for a periodic spline).  For an axisymmetric tokamak 4
    # zeta points are sufficient; all phi slices are identical.
    convert_boozmn_to_chartmap(boozmn_nc, output_nc, nrho=50, ntheta=64, nzeta=4)
finally:
    if os.path.exists(boozmn_nc):
        os.unlink(boozmn_nc)

print(f"chartmap written to {output_nc}")
