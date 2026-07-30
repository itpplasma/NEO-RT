"""Behavioral tests for the full-radius jxb comparison plot."""

import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "python/plot_jxb_mars_comparison.py"
SPEC = importlib.util.spec_from_file_location("plot_jxb", SCRIPT)
PLOT_JXB = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(PLOT_JXB)


def test_linear_profile_integral_and_render(tmp_path: Path) -> None:
    edges, mars = PLOT_JXB.load_profiles(
        ROOT / "test/fixtures/jxb_linear_edges.out",
        ROOT / "test/fixtures/jxb_linear_profile.out",
    )
    cumulative = PLOT_JXB.cumulative_torque(edges, mars[:, 1])
    expected = 6.4 * np.pi**2
    assert np.isclose(cumulative[-1], expected, rtol=0.0, atol=1.0e-13)

    reconstruction = np.column_stack((mars[:, 0], mars[:, 1], mars[:, 1]))
    output = tmp_path / "comparison.png"
    metrics = PLOT_JXB.plot_comparison(
        edges, mars, reconstruction, output, "Manufactured profile"
    )
    assert output.is_file()
    assert output.with_suffix(".pdf").is_file()
    assert metrics["points"] == 3
    assert metrics["max_profile_residual"] == 0.0
