#!/bin/bash
# Larger-sample check of the deep-strip simulation: 20000 profiles (seed 5) at four depths for
# G = 0.015 and 0.03, to compare with the 1000-profile shares of engine/anaxia_deep.json.
# Run from anywhere: bash checks/c/deep_mc_check.sh > checks/c/deep_mc_check.txt
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
gcc -O2 -o "$TMP/deep_mc" "$DIR/deep_mc.c" -lm
for d in 0.04 0.06 0.08 0.09; do for g in 0.015 0.03; do
  echo "d=$d G=$g $("$TMP/deep_mc" "$d" 20000 "$g" 0 12 5 10 0)"
done; done
