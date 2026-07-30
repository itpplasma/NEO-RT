from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np


PATH = (
    Path(__file__).parents[1]
    / "python"
    / "plot_iter_mars_shielding_current.py"
)
SPEC = spec_from_file_location("plot_iter_mars_shielding_current", PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_spectral_norm_matches_explicit_complex_oracle() -> None:
    coefficients = np.array(
        [[3.0 + 4.0j, 0.0], [1.0j, 2.0 + 0.0j]], dtype=complex
    )
    assert np.allclose(
        MODULE.spectral_norm(coefficients), [5.0, np.sqrt(5.0)]
    )


def test_rational_surfaces_match_linear_q_oracle() -> None:
    stor = np.linspace(0.0, 1.0, 101)
    q = 1.0 + 3.0 * stor
    surfaces = MODULE.rational_surfaces(stor, q, -3)
    assert list(surfaces) == list(range(3, 13))
    for mode, surface in surfaces.items():
        assert np.isclose(surface, (mode / 3.0 - 1.0) / 3.0, atol=2.0e-16)
