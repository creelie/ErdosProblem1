#!/usr/bin/env bash
# build_paper.sh -- build the two papers and the release files.
#
# For each paper NAME, with its source in DIR:
#   ErdosStrausSieveConstant   paper/              "A power-saving asymptotic for the two-parameter
#                                                   Type II sieve sum of the Erdős–Straus equation"
#   ErdosStrausExceptionalSet  exceptional/paper/  "An explicit exponent for the exceptional set of
#                                                   the Erdős–Straus equation"
# it writes
#   dist/NAME_final2.pdf             the paper
#   dist/NAME_final2-tex.zip         LaTeX source (main.tex) with the figures as PNG
#   dist/NAME_final2-arxiv.tar.gz    the same source as an arXiv tarball
#   preprintNAME/NAME.pdf            the paper (the only file there)
#   MasterNAME_masterzipfinal2.zip   every tracked file with dist/ and both preprint folders
#
# Needs pdflatex (amsart for the first paper, revtex4-2 for the second), git, zip and tar. The figures
# are drawn by figures/make_figures.py and figures/make_exceptional.py.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
for t in pdflatex git zip tar zipinfo; do
  command -v "$t" >/dev/null || { echo "missing tool: $t" >&2; exit 2; }
done

PAPERS="ErdosStrausSieveConstant:paper ErdosStrausExceptionalSet:exceptional/paper"

latex() { pdflatex -interaction=nonstopmode -halt-on-error "$@" > /dev/null; }
build() { # build <jobname> <input>, run inside the paper directory
  rm -f "$1.aux" "$1.out" "$1.log" "$1.toc"
  for run in 1 2 3 4; do
    latex -jobname="$1" "$2"
    grep -q "Rerun to get" "$1.log" || break
  done
  if grep -E "Overfull|Underfull|undefined|Rerun to get" "$1.log"; then
    echo "build_paper.sh: fix the warnings above" >&2
    exit 1
  fi
  rm -f "$1.aux" "$1.log" "$1.out" "$1.toc"
}

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
rm -rf dist
mkdir -p dist

for entry in $PAPERS; do
  NAME="${entry%%:*}"
  DIR="${entry#*:}"
  PRE=preprint$NAME

  echo "== $NAME: paper"
  ( cd "$DIR" && build main main.tex )
  rm -rf "$PRE"
  mkdir -p "$PRE"
  cp "$DIR/main.pdf" "dist/${NAME}_final2.pdf"
  cp "$DIR/main.pdf" "$PRE/$NAME.pdf"

  echo "== $NAME: tex zip and arXiv tarball"
  SRC="$STAGE/${NAME}_final2-tex"
  mkdir -p "$SRC/figures"
  cp "$DIR/main.tex" "$SRC/"
  cp "$DIR"/figures/*.png "$SRC/figures/"
  ( cd "$SRC" && latex main.tex && latex main.tex && rm -f main.aux main.log main.out main.toc main.pdf mainNotes.bib )
  ( cd "$STAGE" && zip -qrX "${NAME}_final2-tex.zip" "${NAME}_final2-tex" )
  mv "$STAGE/${NAME}_final2-tex.zip" dist/
  entries="$(zipinfo -1 "dist/${NAME}_final2-tex.zip" | wc -l)"
  [ "$entries" -le 512 ] || { echo "build_paper.sh: the tex zip has $entries entries (at most 512)" >&2; exit 1; }
  tar -czf "dist/${NAME}_final2-arxiv.tar.gz" -C "$SRC" main.tex figures
done

for entry in $PAPERS; do
  NAME="${entry%%:*}"
  MASTER=Master${NAME}_masterzipfinal2.zip
  echo "== $MASTER"
  M="$STAGE/master-$NAME/Master$NAME"
  mkdir -p "$M"
  git ls-files -z | grep -zvE '^(dist/|preprint[^/]*/|Master[^/]*_masterzipfinal2\.zip$)' \
    | while IFS= read -r -d '' f; do
        [ -e "$f" ] || continue
        mkdir -p "$M/$(dirname "$f")"
        cp -p "$f" "$M/$f"
      done
  mkdir -p "$M/dist"
  cp dist/* "$M/dist/"
  for e in $PAPERS; do
    n="${e%%:*}"
    mkdir -p "$M/preprint$n"
    cp "preprint$n/$n.pdf" "$M/preprint$n/"
  done
  rm -f "$MASTER"
  ( cd "$STAGE/master-$NAME" && zip -qrX "$ROOT/$MASTER" "Master$NAME" )
done

echo "built: $(ls dist | tr '\n' ' ')"
