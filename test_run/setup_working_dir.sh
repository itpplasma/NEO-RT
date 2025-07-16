#!/bin/bash

# Script to set up NEO-RT working directory for real physics testing
# This script copies the necessary input files for running NEO-RT with real physics

echo "Setting up NEO-RT working directory for real physics testing..."

# Source files from examples
NEORT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXAMPLES_DIR="$NEORT_ROOT/examples/base"

# Required files
FILES=(
    "in_file"
    "driftorbit.in"
    "plasma.in"
)

# Copy files
for file in "${FILES[@]}"; do
    if [ -f "$EXAMPLES_DIR/$file" ]; then
        echo "Copying $file..."
        cp "$EXAMPLES_DIR/$file" .
    elif [ -f "$NEORT_ROOT/$file" ]; then
        echo "Copying $file from root..."
        cp "$NEORT_ROOT/$file" .
    else
        echo "ERROR: $file not found in $EXAMPLES_DIR or $NEORT_ROOT"
        exit 1
    fi
done

echo "Working directory setup complete!"
echo ""
echo "Files copied:"
for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file (failed to copy)"
    fi
done

echo ""
echo "To verify the setup, run:"
echo "  gfortran -o test_orbit_trajectory_comparison test_orbit_trajectory_comparison.f90"
echo "  ./test_orbit_trajectory_comparison"