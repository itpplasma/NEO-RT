from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pytest


PATH = (
    Path(__file__).parents[1]
    / "python"
    / "plot_mars_ampere_comparison.py"
)
SPEC = spec_from_file_location("plot_mars_ampere_comparison", PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_radial_profiles_use_complex_spectral_norm_oracle() -> None:
    actual = np.array(
        [[3.0 + 4.0j, 0.0], [1.0j, 2.0 + 0.0j]], dtype=complex
    )
    predicted = np.array(
        [[0.0, 0.0], [1.0j, -1.0 + 0.0j]], dtype=complex
    )
    actual_norm, predicted_norm, residual = MODULE.radial_profiles(
        actual, predicted
    )
    assert np.allclose(actual_norm, [5.0, np.sqrt(5.0)])
    assert np.allclose(predicted_norm, [0.0, np.sqrt(2.0)])
    assert np.allclose(residual, [1.0, 3.0 / np.sqrt(5.0)])


def test_zero_mars_reference_is_rejected() -> None:
    with pytest.raises(ValueError, match="identically zero"):
        MODULE.radial_profiles(
            np.zeros((3, 2), dtype=complex),
            np.ones((3, 2), dtype=complex),
        )
