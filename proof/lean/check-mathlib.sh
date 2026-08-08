#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
if lake env lean MathlibGate.lean >/dev/null 2>&1; then
  echo "Mathlib import is available."
else
  echo "BLOCKED: Mathlib is unavailable in the configured Lean search path." >&2
  echo "Install and pin mathlib before enabling MathlibGate.lean." >&2
  exit 2
fi
