#!/usr/bin/env bash
# Build (if needed) and run the three example tidal components from the
# paper (DW1, SW2, TW3), producing fortran/output/hough_s*.dat files that
# scripts/plot_paper_figures.py turns into Figures 1-3.
set -euo pipefail

cd "$(dirname "$0")/.."   # fortran/

if [ ! -x build/hough_main ]; then
  echo "== configuring build =="
  cmake -B build -S .
fi
echo "== building hough_main =="
cmake --build build

for preset in dw1 sw2 tw3; do
  echo ""
  echo "== running preset: $preset =="
  ./build/hough_main --preset="$preset"
done

echo ""
echo "Done. Now run: python3 scripts/plot_paper_figures.py"
