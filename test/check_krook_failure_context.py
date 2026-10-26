"""Check the native transport rejection retains status and physical context."""
import math
import re
import subprocess
import sys


result = subprocess.run([sys.argv[1], "failure-context"], capture_output=True,
                        text=True, timeout=5)
assert result.returncode != 0, "invalid bounds were accepted"
output = result.stdout + result.stderr
assert "ERROR STOP finite-width Krook quadrature did not converge" in output, output
records = [line for line in output.splitlines()
           if line.startswith("KROOK_QUADRATURE_FAILURE ")]
assert len(records) == 1, output
fields = dict(re.findall(r"(\S+)=\s*([^\s]+)", records[0]))
expected = {
    "status": 1, "species_index": 3, "mth": -2, "mph": 3,
    "eta_steps": 16, "max_depth": 20, "s": .09,
    "mass_g": 9.1091e-28, "charge_statc": -4.803204e-10,
    "density_cm-3": 1e10, "omega_E_s-1": -708.125, "ux": .375,
    "vth_cm_s": 2.6524515369350795e9, "nu_s-1": 7,
    "eta_min": .25, "eta_max": .25, "eta_tp": .125, "sign_vpar": -1,
    "effective_harmonic": -2, "geometry": 2, "rtol": 1e-5, "atol": 0,
}
assert fields.keys() == expected.keys(), (fields, expected)
for key, value in expected.items():
    assert math.isclose(float(fields[key]), value, rel_tol=2e-16, abs_tol=0), (
        key, fields[key], value)
estimates = [line for line in output.splitlines()
             if line.startswith("KROOK_QUADRATURE_ESTIMATE ")]
assert len(estimates) == 1, output
estimate = dict(re.findall(r"(\S+)=\s*([^\s]+)", estimates[0]))
assert estimate == {
    "status": "1", "channels": "3", "accepted_leaves": "0",
    "exhausted_leaves": "0", "global_finite": "T", "leaf_finite": "T",
    "sum_abs_saturated": "F",
}, estimate
channels = [dict(re.findall(r"(\S+)=\s*([^\s]+)", line))
            for line in output.splitlines()
            if line.startswith("KROOK_QUADRATURE_CHANNEL ")]
assert len(channels) == 3, output
for index, channel in enumerate(channels, start=1):
    assert channel.pop("channel") == str(index), channel
    assert channel.keys() == {
        "value", "error", "threshold", "exhausted_error", "sum_abs_leaf"}, channel
    assert all(float(value) == 0 for value in channel.values()), channel
assert "KROOK_QUADRATURE_WORST recorded=F" in output, output
assert "KROOK_QUADRATURE_LEAF " not in output, output
assert "KROOK_QUADRATURE_INVERSE " not in output, output
print("native Krook rejection preserves status and physical context")
