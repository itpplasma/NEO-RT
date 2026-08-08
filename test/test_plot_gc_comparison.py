"""Behavioral checks for comparison-plot model and signed-axis contracts."""

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest


MODULE_PATH = Path(__file__).parents[1] / "tools" / "plot_gc_comparison.py"
SPEC = importlib.util.spec_from_file_location("plot_gc_comparison", MODULE_PATH)
plot_gc_comparison = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = plot_gc_comparison
SPEC.loader.exec_module(plot_gc_comparison)


def test_plot_metadata_rejects_removed_model_and_preserves_signed_zero(tmp_path):
    assert plot_gc_comparison.canonical_model("legacy") == "boozer_thin"
    assert plot_gc_comparison.canonical_model("0") == "boozer_thin"
    assert plot_gc_comparison.canonical_model("neort_full") == "full_fow"
    assert plot_gc_comparison.canonical_model("2") == "full_fow"
    assert plot_gc_comparison.model_label("potato", torque=True) == (
        "POTATO (magnitude-only; conjugate-pair mode power)"
    )
    with pytest.raises(ValueError, match="real-space thin"):
        plot_gc_comparison.canonical_model("gc_thin")
    with pytest.raises(ValueError, match="real-space thin"):
        plot_gc_comparison.canonical_model("1")

    frequency = tmp_path / "frequency_inventory.out"
    frequency.write_text(
        "# frequency_inventory_v1\n"
        "0 passing 1 0.2 0.2 1.0 10.0 -2.0 -3.0 4.0 0 0\n"
    )
    assert plot_gc_comparison.read_frequency(frequency) == [
        ("boozer_thin", "passing", 1.0, 10.0, -2.0, -3.0, 4.0)
    ]

    run = tmp_path / "torque"
    run.mkdir()
    np.savetxt(run / "legacy_torque.out", [[0.25, 1.0, 0.1, -2.0, 3.0, -4.0]])
    np.savetxt(run / "gc_torque.out", [[0.25, 1.0, 0.1, -1.0, 2.0, -3.0]])
    (run / "integral_torque.dat").write_text("-5.0\n")
    figure = plot_gc_comparison.make_torque_figure(run)
    assert len(figure.axes) == 2
    assert figure.axes[0].get_yscale() == "symlog"
    assert figure.axes[1].get_ylim()[0] == 0.0
    signed_labels = [text.get_text() for text in figure.axes[0].get_legend().get_texts()]
    assert signed_labels == ["Boozer thin", "Full real-space FOW"]
    magnitude_labels = [
        text.get_text() for text in figure.axes[1].get_legend().get_texts()
    ]
    assert magnitude_labels == ["POTATO magnitude-only"]
    assert "not torque density" in figure.axes[1].get_title()
    assert any("sign gate" in text.get_text() for text in figure.axes[1].texts)
    assert any("conjugate-pair" in text.get_text() for text in figure.texts)
    plt.close(figure)

    (run / "integral_torque.dat").write_text("1.0 2.0\n")
    with pytest.raises(ValueError, match="one finite scalar"):
        plot_gc_comparison.make_torque_figure(run)

    reports = tmp_path / "reports"
    reports.mkdir()
    current_frequency = reports / "case_frequency_inventory.out"
    current_resonance = reports / "case_resonance_inventory.out"
    current_frequency.write_text("")
    current_resonance.write_text("")
    assert plot_gc_comparison.discover_report(
        reports, "*_frequency_inventory.out", "frequency_scan.out"
    ) == current_frequency
    assert plot_gc_comparison.discover_report(
        reports, "*_resonance_inventory.out", "frequency_resonances.out"
    ) == current_resonance
    (reports / "frequency_scan.out").write_text("")
    with pytest.raises(ValueError, match="ambiguous report"):
        plot_gc_comparison.discover_report(
            reports, "*_frequency_inventory.out", "frequency_scan.out"
        )

    (run / "gc_thin_torque.out").write_text("0 0 0 0 0 0\n")
    with pytest.raises(ValueError, match="deprecated ambiguous"):
        plot_gc_comparison.reject_deprecated_model_files(run)

    figure, axes = plt.subplots()
    plot_gc_comparison.configure_signed_axis(axes, np.array([-4.0, 2.0]))
    assert axes.get_yscale() == "symlog"
    lower, upper = axes.get_ylim()
    assert lower < 0.0 < upper
    assert any(np.allclose(line.get_ydata(), 0.0) for line in axes.lines)
    plt.close(figure)
