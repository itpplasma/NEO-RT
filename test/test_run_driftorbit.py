import importlib.util
import subprocess
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "python" / "run_driftorbit.py"
SPEC = importlib.util.spec_from_file_location("run_driftorbit", SCRIPT)
assert SPEC is not None
assert SPEC.loader is not None
run_driftorbit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_driftorbit)


def make_executable(path: Path, body: str) -> Path:
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(0o755)
    return path


def write_template(path: Path) -> Path:
    path.write_text("s = <S_TOKEN>\n")
    return path


def test_absolute_executable_receives_runname(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    marker = tmp_path / "argv.txt"
    executable = make_executable(
        tmp_path / "fake-success",
        f'printf "%s" "$1" > "{marker}"\n',
    )
    monkeypatch.chdir(tmp_path)

    run_driftorbit.run_single_flux_surface(
        str(executable), str(write_template(tmp_path / "template.in")),
        "surface0", 0.1, 0.2, 0.3, 1.0,
    )

    assert marker.read_text() == "surface0"


def test_single_failure_is_propagated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    executable = make_executable(
        tmp_path / "fake-failure",
        "echo 'solver failed' >&2\nexit 7\n",
    )
    monkeypatch.chdir(tmp_path)

    with pytest.raises(subprocess.CalledProcessError) as error:
        run_driftorbit.run_single_flux_surface(
            str(executable), str(write_template(tmp_path / "template.in")),
            "surface0", 0.1, 0.2, 0.3, 1.0,
        )

    assert error.value.returncode == 7


def test_multi_surface_failure_is_propagated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    executable = make_executable(
        tmp_path / "fake-failure",
        "echo 'solver failed' >&2\nexit 7\n",
    )
    profile = tmp_path / "profile.in"
    profile.write_text("0.1 0.2 0.3\n0.4 0.5 0.6\n")
    monkeypatch.chdir(tmp_path)

    with pytest.raises(subprocess.CalledProcessError) as error:
        run_driftorbit.run_multiple_flux_surfaces(
            str(executable), str(profile), str(write_template(tmp_path / "template.in")),
            "surface", 0, 2,
        )

    assert error.value.returncode == 7
