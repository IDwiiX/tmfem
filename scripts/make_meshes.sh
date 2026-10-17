#!/usr/bin/env bash
# Regenerates every committed .msh file from its .geo (Gmsh 4.12.1, ASCII msh 4.1).
# Run from anywhere:  bash scripts/make_meshes.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED="4.12.1"
VERSION="$(gmsh --version 2>&1 | tail -n 1)"
if [[ "$VERSION" != "$EXPECTED" ]]; then
  echo "make_meshes.sh: Gmsh $EXPECTED expected, found $VERSION (meshes would differ)" >&2
  exit 1
fi

mesh() {  # geo order clscale output
  gmsh "$1" -2 -order "$2" -clscale "$3" -format msh41 -v 2 -o "$4" > /dev/null
  echo "wrote ${4#"$ROOT"/}"
}

# Kirsch plate: 5 levels, each halving every element size; T3 (order 1) and T6 (order 2).
# Level 4 was added after the first P2-f run showed level 0 is pre-asymptotic (README, P2-f note).
SCALES=(1 0.5 0.25 0.125 0.0625)
for level in 0 1 2 3 4; do
  for order in 1 2; do
    t=$((order == 1 ? 3 : 6))
    mesh "$ROOT/cases/kirsch/kirsch.geo" "$order" "${SCALES[$level]}" \
      "$ROOT/cases/kirsch/msh/kirsch_t${t}_l${level}.msh"
  done
done

for order in 1 2; do
  t=$((order == 1 ? 3 : 6))
  mesh "$ROOT/cases/bar/bar.geo" "$order" 1 "$ROOT/cases/bar/bar_t${t}.msh"
done
