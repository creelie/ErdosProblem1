#!/usr/bin/env bash
# run_all.sh -- builds the Lean project, checks its axioms and runs the independent checks.
#
# Usage: scripts/run_all.sh             (all steps)
#        scripts/run_all.sh --no-lean   (skip the Lean build)
# Environment:
#   CHECKS="python julia c cpp engine"   the checks to run (default: all five; "" runs none)
#   SKIP_LEAN=1                   same as --no-lean
# Requirements: elan (Lean 4) for the Lean build; python3, julia, a C compiler and a C++20
# compiler for the checks that are run; the engine check also needs network access to install the
# engine Anaxia (creelie/eBLLM, version 3.0.5) into a virtual environment.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(mktemp -d)"
trap 'rm -rf "$BUILD"' EXIT
CC_BIN="${CC:-cc}"
CXX_BIN="${CXX:-c++}"
CHECKS="${CHECKS-python julia c cpp engine}"
SKIP_LEAN="${SKIP_LEAN:-0}"
[ "${1:-}" = "--no-lean" ] && SKIP_LEAN=1

step() { printf '\n=== %s\n' "$1"; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "missing tool: $1" >&2; exit 2; }; }

for c in $CHECKS; do
  case "$c" in
    python) need python3 ;;
    julia) need julia ;;
    c) need "$CC_BIN" ;;
    cpp) need "$CXX_BIN" ;;
    engine) need python3 ;;
    *) echo "unknown check: $c" >&2; exit 2 ;;
  esac
done

if [ "$SKIP_LEAN" != 1 ]; then
  need lake
  step "Lean 4: lake build"
  ( cd "$ROOT/lean" && lake exe cache get >/dev/null 2>&1 || true )
  ( cd "$ROOT/lean" && lake build 2>&1 | tee "$BUILD/lean.log" | grep -E "error|Build completed|depends on axioms" ) || true
  if ! grep -q "Build completed successfully" "$BUILD/lean.log" || grep -q "error" "$BUILD/lean.log"; then
    echo "Lean build failed" >&2; exit 1
  fi
  step "Lean 4: axioms of the theorems"
  # every theorem listed in Axioms.lean may use only the three standard axioms
  count=0
  while IFS= read -r thm; do
    line="$(grep -F "'ErdosStraus.$thm' d" "$BUILD/lean.log" | head -n 1 | sed 's/^info: Axioms.lean:[0-9]*:[0-9]*: //' || true)"
    echo "$line"
    case "$line" in
      "'ErdosStraus.$thm' does not depend on any axioms") ;;
      "'ErdosStraus.$thm' depends on axioms: ["*"]")
        list="${line#*\[}"; list="${list%\]}"
        for ax in ${list//,/ }; do
          case "$ax" in propext|Classical.choice|Quot.sound) ;;
            *) echo "unexpected axiom $ax for $thm" >&2; exit 1 ;; esac
        done ;;
      *) echo "no axiom report for $thm" >&2; exit 1 ;;
    esac
    count=$((count + 1))
  done < <(sed -n 's/^#print axioms ErdosStraus\.//p' "$ROOT/lean/Axioms.lean")
  echo "$count theorems, each on the standard axioms at most"
  if grep -rnwE "sorry|admit|native_decide|bv_decide" "$ROOT/lean/ErdosStraus" "$ROOT/lean/ErdosStraus.lean"; then
    echo "sorry, admit, native_decide or bv_decide found" >&2; exit 1
  fi
  if grep -rnE "^[[:space:]]*axiom[[:space:]]" "$ROOT/lean/ErdosStraus" "$ROOT/lean/ErdosStraus.lean"; then
    echo "axiom declaration found" >&2; exit 1
  fi
  echo "no sorry, admit, native_decide, bv_decide or axiom declarations"
fi

for c in $CHECKS; do
  case "$c" in
    python)
      step "Python"
      python3 -I "$ROOT/checks/python/check.py"
      python3 -I "$ROOT/checks/python/exceptional.py" ;;
    julia)
      step "Julia"
      julia "$ROOT/checks/julia/check.jl"
      julia "$ROOT/checks/julia/tilt.jl" ;;
    c)
      step "C"
      "$CC_BIN" -O2 -std=c99 -o "$BUILD/check_c" "$ROOT/checks/c/check.c" -lm
      "$BUILD/check_c"
      "$CC_BIN" -O2 -std=c99 -o "$BUILD/constant_M" "$ROOT/checks/c/constant_M.c" -lm
      "$BUILD/constant_M" 100000000
      "$CC_BIN" -O2 -std=c99 -o "$BUILD/exceptional" "$ROOT/checks/c/exceptional.c" -lm
      "$BUILD/exceptional" ;;
    cpp)
      step "C++"
      "$CXX_BIN" -O2 -std=c++20 -o "$BUILD/check_cpp" "$ROOT/checks/cpp/check.cpp"
      "$BUILD/check_cpp"
      "$CXX_BIN" -O2 -std=c++20 -o "$BUILD/mean_values" "$ROOT/checks/cpp/mean_values.cpp"
      "$BUILD/mean_values" ;;
    engine)
      step "Engine Anaxia: the run behind the three corners"
      python3 -m venv "$BUILD/venv"
      "$BUILD/venv/bin/pip" install --quiet "git+https://github.com/creelie/eBLLM@v3.0.5"
      ( cd "$BUILD" && "$BUILD/venv/bin/python" "$ROOT/engine/anaxia_regions.py" 1000000 >/dev/null )
      python3 -I "$ROOT/engine/compare.py" "$BUILD/anaxia_regions.json" "$ROOT/engine/anaxia_regions.json" ;;
  esac
done

printf '\nAll steps passed.\n'
