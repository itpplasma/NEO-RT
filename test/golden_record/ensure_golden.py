#!/usr/bin/env python3
"""
Ensure golden record data exists by building from main branch if needed.

This module handles:
1. Cloning/pulling the main branch to a reference directory
2. Building the main branch with the same CONFIG as the current build
3. Regenerating golden.h5 if the main branch or CONFIG has changed

The reference directory is cached locally to avoid unnecessary rebuilds.
The golden.h5 is regenerated when either the main branch commit or the
build CONFIG changes, ensuring tests compare against matching builds.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
MAIN_REF_DIR = SCRIPT_DIR / "main_ref"
GOLDEN_H5 = SCRIPT_DIR / "golden.h5"
LAST_BUILD_STATE_FILE = MAIN_REF_DIR / ".last_build_state"
REPO_URL = "https://github.com/itpplasma/NEO-RT"
DEFAULT_CONFIG = "Fast"


def run_cmd(cmd: list[str], cwd: Path | None = None, check: bool = True) -> str:
    """Run a command and return stdout."""
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=check)
    return result.stdout.strip()


def create_conftest() -> None:
    """Create conftest.py to exclude main_ref from pytest collection."""
    conftest = MAIN_REF_DIR / "conftest.py"
    conftest.write_text('collect_ignore_glob = ["**/*"]\n')


def clone_main_ref() -> None:
    """Clone the main branch to the reference directory."""
    print(f"Cloning main branch to {MAIN_REF_DIR}...")
    subprocess.run(
        [
            "git",
            "clone",
            "--branch",
            "main",
            "--single-branch",
            REPO_URL,
            str(MAIN_REF_DIR),
        ],
        check=True,
    )
    create_conftest()


def pull_main_ref() -> bool:
    """Pull latest changes. Returns True if there were updates."""
    print("Fetching latest changes from main...")
    run_cmd(["git", "fetch", "origin", "main"], cwd=MAIN_REF_DIR)

    local_head = run_cmd(["git", "rev-parse", "HEAD"], cwd=MAIN_REF_DIR)
    remote_head = run_cmd(["git", "rev-parse", "origin/main"], cwd=MAIN_REF_DIR)

    if local_head != remote_head:
        print(f"Updating from {local_head[:8]} to {remote_head[:8]}...")
        run_cmd(["git", "reset", "--hard", "origin/main"], cwd=MAIN_REF_DIR)
        return True

    print("Main branch is up to date.")
    return False


def detect_current_config() -> str:
    """Detect the CONFIG used in the current build from CMakeCache.txt."""
    # Look for build directory relative to the repo root
    repo_root = SCRIPT_DIR.parent.parent
    cmake_cache = repo_root / "build" / "CMakeCache.txt"

    if not cmake_cache.exists():
        print(f"No CMakeCache.txt found at {cmake_cache}, using {DEFAULT_CONFIG}")
        return DEFAULT_CONFIG

    content = cmake_cache.read_text()
    match = re.search(r"CMAKE_BUILD_TYPE:STRING=(\w+)", content)
    if match:
        config = match.group(1)
        print(f"Detected build CONFIG: {config}")
        return config

    print(f"Could not detect CONFIG from CMakeCache.txt, using {DEFAULT_CONFIG}")
    return DEFAULT_CONFIG


def get_last_config() -> str | None:
    """Get the CONFIG that was last used to build main_ref."""
    if LAST_BUILD_STATE_FILE.exists():
        return LAST_BUILD_STATE_FILE.read_text().strip()
    return None


def save_last_config(config: str) -> None:
    """Save the CONFIG that was just used to build main_ref."""
    LAST_BUILD_STATE_FILE.write_text(config)


def build_main_ref(config: str) -> Path:
    """Build the main branch with the given CONFIG and return path to executable."""
    print(f"Building main branch with CONFIG={config}...")
    # Clean first to ensure fresh build with correct config
    subprocess.run(["make", "clean"], cwd=MAIN_REF_DIR, check=False)
    subprocess.run(
        ["make", f"CONFIG={config}"],
        cwd=MAIN_REF_DIR,
        check=True,
    )
    return MAIN_REF_DIR / "build" / "neo_rt.x"


def regenerate_golden(executable: Path) -> None:
    """Regenerate golden.h5 using the given executable."""
    print("Regenerating golden.h5...")
    sys.path.insert(0, str(SCRIPT_DIR))
    from regenerate_golden import main as regenerate_main

    original_argv = sys.argv
    sys.argv = ["regenerate_golden.py", str(executable)]
    try:
        regenerate_main()
    finally:
        sys.argv = original_argv


def ensure_golden() -> Path:
    """
    Ensure golden.h5 exists and is up to date with main branch.

    The golden record is rebuilt when either:
    - The remote main branch has new commits
    - The build CONFIG has changed (Fast vs Debug)

    Returns the path to golden.h5.
    """
    current_config = detect_current_config()

    needs_clone = not MAIN_REF_DIR.exists()
    needs_build = False
    needs_regenerate = not GOLDEN_H5.exists()

    if needs_clone:
        clone_main_ref()
        needs_build = True
        needs_regenerate = True
    else:
        updated = pull_main_ref()
        if updated:
            needs_build = True
            needs_regenerate = True

    # Rebuild if config changed
    last_config = get_last_config()
    if last_config != current_config:
        print(f"CONFIG changed: {last_config} -> {current_config}")
        needs_build = True
        needs_regenerate = True

    executable = MAIN_REF_DIR / "build" / "neo_rt.x"

    if needs_build or not executable.exists():
        executable = build_main_ref(current_config)
        save_last_config(current_config)
        needs_regenerate = True

    if needs_regenerate:
        regenerate_golden(executable)

    if not GOLDEN_H5.exists():
        raise RuntimeError(f"Failed to create {GOLDEN_H5}")

    print(f"Golden record ready: {GOLDEN_H5} (CONFIG={current_config})")
    return GOLDEN_H5


def main() -> None:
    """CLI entry point."""
    ensure_golden()


if __name__ == "__main__":
    main()
