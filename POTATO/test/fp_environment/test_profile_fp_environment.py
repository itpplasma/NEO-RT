"""The hot profile evaluators must not save and restore the FP environment.

gfortran wraps every procedure that can reach ieee_arithmetic, even through a
module imported with an unrelated ONLY list, in _gfortran_ieee_procedure_entry/
_exit.  phielec_of_psi and denstemp_of_psi run on every field evaluation, so a
module-scope IEEE import anywhere in their use chain costs ~45% of POTATO's CPU
time while changing no result.

usage: test_profile_fp_environment.py OBJECT...
"""
import re
import shutil
import subprocess
import sys

HOT = ("phielec_of_psi_", "denstemp_of_psi_")
SAVE = "_gfortran_ieee_procedure_entry"
SKIP = 77


def calls_by_function(objects):
    """Map each function symbol in the objects to its disassembly text."""
    bodies, current = {}, None
    for obj in objects:
        listing = subprocess.run(["objdump", "-dr", "--no-show-raw-insn", obj],
                                 capture_output=True, text=True, check=True).stdout
        for line in listing.splitlines():
            header = re.match(r"^[0-9a-f]+ <(.+)>:$", line)
            if header:
                current = header.group(1)
                bodies.setdefault(current, [])
            elif current:
                bodies[current].append(line)
    return {name: "\n".join(lines) for name, lines in bodies.items()}


def only_lto_bytecode(obj):
    """True for a -flto object, whose code is GIMPLE until the final link."""
    sections = subprocess.run(["objdump", "-h", obj], capture_output=True,
                              text=True, check=True).stdout
    return ".gnu.lto_." in sections


def main(objects):
    if shutil.which("objdump") is None:
        print("objdump not available")
        return SKIP
    bodies = calls_by_function(objects)
    missing = [name for name in HOT if name not in bodies]
    if missing and any(only_lto_bytecode(obj) for obj in objects):
        print("objects hold LTO bytecode only (POTATO_ENABLE_LTO); nothing to inspect")
        return SKIP
    if missing:
        print("profile evaluators not found in", objects, ":", missing)
        return 1
    saving = [name for name in HOT if SAVE in bodies[name]]
    for name in saving:
        print(f"{name} saves/restores the FP environment: a module in its use "
              "chain imports ieee_arithmetic at module scope")
    return 1 if saving else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
