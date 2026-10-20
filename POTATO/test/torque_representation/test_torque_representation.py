"""Behavioral oracles for signed harmonics and the thermal kinetic-energy floor."""
import array
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

EXE = Path(os.environ.get("POTATO_EXE", "POTATO/build/potato.x")).resolve()
CASE = Path(__file__).resolve().parent.parent / "golden_record_resonance"
FILES = ("potato.in", "field_divB0.inp", "convexwall.dat", "profile_poly.in",
         "circ.eqdsk", "bmod_n.dat")


def conjugate_perturbation(path):
    """Conjugate the complex harmonic for the equivalent negative-n convention."""
    data = bytearray(path.read_bytes())
    offset = 0
    for _ in range(2):
        offset += struct.unpack_from("=i", data, offset)[0] + 8
    size = struct.unpack_from("=i", data, offset)[0]
    start = offset + 4 + size // 2
    imag = array.array("d")
    imag.frombytes(data[start:offset + 4 + size])
    assert any(imag), "the sign oracle needs a complex perturbation"
    data[start:offset + 4 + size] = array.array("d", (-x for x in imag)).tobytes()
    path.write_bytes(data)


class TorqueRepresentation(unittest.TestCase):
    def run_case(self, conjugate=False, **options):
        with tempfile.TemporaryDirectory(prefix="potato-representation-") as tmp:
            work = Path(tmp)
            for name in FILES:
                shutil.copy(CASE / name, work / name)
            text = (work / "potato.in").read_text()
            for name, value in options.items():
                text, count = re.subn(rf"\b{name}\s*=\s*[^\n]+", f"{name} = {value}", text)
                if not count:
                    text = text.replace("\n/", f"\n  {name} = {value}\n/")
            (work / "potato.in").write_text(text)
            if conjugate:
                conjugate_perturbation(work / "bmod_n.dat")
            result = subprocess.run([EXE], cwd=work, capture_output=True, text=True,
                                    timeout=300, env=dict(os.environ, OMP_NUM_THREADS="2"))
            outputs = {name: (work / name).read_text() for name in
                       ("integral_torque.dat", "boxcounted_torque.dat", "fort.31415")
                       if (work / name).is_file()}
        return result, outputs

    def successful(self, **options):
        result, outputs = self.run_case(**options)
        self.assertEqual(result.returncode, 0, result.stdout[-1500:] + result.stderr[-1500:])
        return result.stdout, outputs

    def number(self, text, label):
        return float(re.search(rf"^{label} =\s*(\S+)", text, re.M).group(1))

    def samples(self, text):
        return sorted(float(x) for x in re.findall(r"^toten =\s*(\S+)", text, re.M))

    def test_signed_harmonic(self):
        options = dict(nenerg=1, enkin_min_over_temp="0.1d0")
        _, positive = self.successful(m_min=-1, m_max=3, n_tor=2, **options)
        _, negative = self.successful(m_min=-3, m_max=1, n_tor=-2,
                                      conjugate=True, **options)
        reference = float(positive["integral_torque.dat"])
        self.assertTrue(math.isfinite(reference) and abs(reference) > 1e-30)
        self.assertTrue(math.isclose(reference, float(negative["integral_torque.dat"]),
                                     rel_tol=1e-10))
        for p, n in zip(positive["boxcounted_torque.dat"].splitlines(),
                        negative["boxcounted_torque.dat"].splitlines(), strict=True):
            self.assertTrue(math.isclose(float(p.split()[0]), float(n.split()[0]),
                                         rel_tol=1e-10, abs_tol=abs(reference)*1e-12))
        def canonical_lines(output, sign):
            return sorted(tuple(map(float, row.split()[:3])) +
                          tuple(sign*int(x) for x in row.split()[3:])
                          for row in output["fort.31415"].splitlines())
        self.assertEqual(canonical_lines(positive, 1), canonical_lines(negative, -1))

    def test_kinetic_floor(self):
        text, outputs = self.successful(enkin_min_over_temp="0.1d0", nenerg=3)
        phi_max = self.number(text, "maximum potential energy")
        temperature = self.number(text, "maximum kinetic energy") / 3
        minimum = phi_max + 0.1*temperature
        self.assertAlmostEqual(self.number(text, "minimum total energy"), minimum, delta=2e-5)
        maximum = phi_max + 3*temperature
        expected = [minimum + (i + 0.5)*(maximum - minimum)/3 for i in range(3)]
        used = self.samples(text)
        self.assertEqual(len(used), 3)
        for actual, wanted in zip(used, expected):
            self.assertAlmostEqual(actual, wanted, delta=2e-5)
        slices = sorted({float(row.split()[0]) for row in outputs["fort.31415"].splitlines()})
        self.assertEqual(len(slices), len(used))
        for recorded, sampled in zip(slices, used):
            self.assertAlmostEqual(recorded, sampled, delta=1e-12)
        self.assertTrue(all(energy - phi_max >= 0.1*temperature for energy in used))

    def invalid_floor(self, floor):
        result, _ = self.run_case(enkin_min_over_temp=floor, nenerg=1)
        self.assertNotEqual(result.returncode, 0, "empty/backward thermal grid was accepted")
        self.assertIn("enkin_min_over_temp must be less than thermen_max", result.stderr)

    def test_equal_floor(self):
        self.invalid_floor("3d0")

    def test_excess_floor(self):
        self.invalid_floor("4d0")

    def test_legacy_grid(self):
        text, _ = self.successful(enkin_min_over_temp="0d0", nenerg=3)
        minimum = self.number(text, "minimum potential energy")
        maximum = self.number(text, "maximum total energy")
        self.assertEqual(self.number(text, "minimum total energy"), minimum)
        used = self.samples(text)
        self.assertEqual(len(used), 2)
        for actual, index in zip(used, (2, 3)):
            expected = minimum + (index - 0.5)*(maximum - minimum)/3
            self.assertAlmostEqual(actual, expected, delta=2e-5)

    def test_monoenergetic(self):
        # Thermal limits must not reject or shift the separate delta(x-x0) path.
        options = dict(monoenergetic_x="2d0", nenerg=7)
        baseline, first = self.successful(enkin_min_over_temp="0d0", **options)
        cutoff, second = self.successful(enkin_min_over_temp="4d0", **options)
        self.assertEqual(len(self.samples(baseline)), 1)
        self.assertEqual(self.samples(baseline), self.samples(cutoff))
        self.assertEqual(first, second)


if __name__ == "__main__":
    EXE = Path(sys.argv.pop(1)).resolve()
    unittest.main()
