import importlib.util
from pathlib import Path

import pytest

script = Path(__file__).parents[1] / "MARS-JXB" / "analyze_gpec_torque.py"
spec = importlib.util.spec_from_file_location("analyze_gpec_torque", script)
assert spec is not None
assert spec.loader is not None
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def test_parses_independent_handwritten_gpec_oracle(tmp_path: Path) -> None:
    control = tmp_path / "gpec_control_n3.out"
    control.write_text(
        "GPEC_CONTROL\n"
        " plasma energy =  2.50000000E+001\n"
        " toroidal torque = -1.25000000E+000\n"
    )
    profile = tmp_path / "gpec_dw_profile_n3.out"
    profile.write_text(
        "GPEC_dw: Energy and torque profiles\n"
        " psi T_phi 2ndeltaW int(T_phi) int(2ndeltaW) dv/dpsi_n\n"
        " 0.10000000E+00 -2.00000000E+00 3.0 -0.20000000E+00 4.0 5.0\n"
        " 0.90000000E+00 -1.50000000E+00 6.0 -1.25000000E+00 7.0 8.0\n"
    )

    assert analysis.parse_control_torque(control) == pytest.approx(-1.25)
    rows = analysis.parse_dw_profile(profile)
    assert rows[0] == pytest.approx((0.1, -2.0, 3.0, -0.2, 4.0, 5.0))
    assert rows[1] == pytest.approx((0.9, -1.5, 6.0, -1.25, 7.0, 8.0))


def test_rejects_nonmonotonic_native_grid(tmp_path: Path) -> None:
    profile = tmp_path / "gpec_dw_profile_n3.out"
    profile.write_text(
        "0.2 1 2 3 4 5\n"
        "0.1 1 2 3 4 5\n"
    )

    with pytest.raises(ValueError, match="strictly increasing"):
        analysis.parse_dw_profile(profile)


def test_parses_handwritten_rational_surface_oracle(tmp_path: Path) -> None:
    dcon_log = tmp_path / "dcon.log"
    dcon_log.write_text(
        " psi = 4.894E-01, q = 1.333\n"
        " psi = 6.544E-01, q = 1.667\n"
        " psi = 7.626E-01, q = 2.000\n"
        " psi = 9.901E-01, q = 3.400\n"
    )

    assert analysis.parse_rational_surfaces(dcon_log, toroidal_mode=3) == [
        pytest.approx((0.4894, 4)),
        pytest.approx((0.6544, 5)),
        pytest.approx((0.7626, 6)),
    ]
