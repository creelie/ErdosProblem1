"""sum_K.py -- the certified lower bound for the measure K of the factorable pairs in the deep strips.

Reads checks/c/certify_K.out, one line per cell [k/W, (k+1)/W] of depths, as written by certify_K
(through the sweep that prints "cell k/W rc=... Vpass=..."), and adds up width * Vpass over the cells
that meet [1/42, 1/10].  Vpass is a lower bound for the density k(vartheta) on the whole cell; it is
printed with six decimals, so 10^-6 is subtracted from each value before it is used.  On the first
cell only the part above 1/42 is counted, and beyond 1/10 there is nothing to count (no triple of sizes
is admissible there).  Exact rational arithmetic throughout.

Then, with J < 0.000335 (Proposition of the paper) and pi < 3.141593,
    c1^3 = (27/8) (9 / (8 pi^2)) (40/49 - 6 J + 6 K) > (27/8) (9 / (8 * 3.141593^2)) (40/49 - 6 * 0.000335 + 6 Kb),
where Kb is K rounded down to five decimals; the script prints Kb, the bound for c1 rounded down to four
decimals, the sums by ranges of depth and the cells with their values of G.

Run: python3 sum_K.py [path to certify_K.out]
"""
import os
import re
import sys
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PATH = os.path.join(ROOT, "checks", "c", "certify_K.out")
LO, HI = F(1, 42), F(1, 10)
J_UP = F(335, 10 ** 6)
PI_UP = F(3141593, 10 ** 6)
RANGES = [(24, 36), (36, 46), (46, 62), (62, 82), (82, 103)]


def cells():
    out = {}
    pat = re.compile(r"cell (\d+)/(\d+) rc=(\d+) .*?G=([\d.]+) .*?KMAX=(\d+) Vpass=([\d.]+)")
    for line in open(PATH):
        m = pat.match(line)
        if not m:
            continue
        k, w, rc, g, kmax, v = m.groups()
        if rc != "0":
            continue
        k, w = int(k), int(w)
        val = max(F(v) - F(1, 10 ** 6), F(0))
        key = (F(k, w), F(k + 1, w))
        if key not in out or val > out[key][0]:
            out[key] = (val, g, int(kmax), k, w)
    return out


def bound():
    """The certified sum, its value rounded down to five decimals, the bound for c1 and the details."""
    cs = cells()
    keys = sorted(cs)
    for (a, b), (c, d) in zip(keys, keys[1:]):
        assert b <= c, "cells overlap"
    total = F(0)
    by_range = {r: F(0) for r in RANGES}
    covered = F(0)
    for (a, b) in keys:
        val, g, kmax, k, w = cs[(a, b)]
        lo, hi = max(a, LO), min(b, HI)
        if hi <= lo:
            continue
        part = (hi - lo) * val
        total += part
        covered += hi - lo
        for r in RANGES:
            if F(r[0], 1024) <= a and b <= F(r[1], 1024):
                by_range[r] += part
    Kb = F(int(total * 10 ** 5), 10 ** 5)
    share = F(40, 49) - 6 * J_UP + 6 * Kb
    c3 = F(27, 8) * F(9, 8) / PI_UP ** 2 * share
    cb = F(int((float(c3) ** (1 / 3)) * 10 ** 4), 10 ** 4)
    while cb ** 3 >= c3:
        cb -= F(1, 10 ** 4)
    assert cb ** 3 < c3
    return dict(cells=cs, keys=keys, total=total, by_range=by_range, covered=covered, Kb=Kb, share=share,
                c3=c3, cb=cb)


def main():
    r = bound()
    cs, keys, total, by_range, covered = r["cells"], r["keys"], r["total"], r["by_range"], r["covered"]
    Kb, share, c3, cb = r["Kb"], r["share"], r["c3"], r["cb"]
    print(f"cells: {len(keys)}, depths covered: {float(covered):.6f} of {float(HI - LO):.6f}")
    for r in RANGES:
        print(f"  cells {r[0]}-{r[1] - 1} (depth {r[0] / 1024:.5f} to {r[1] / 1024:.5f}): {float(by_range[r]):.7f}")
    print(f"K > {float(total):.7f}, so K > {Kb} = {float(Kb)}")
    print(f"share > 40/49 - 6*0.000335 + 6*{float(Kb)} = {float(share):.7f}")
    print(f"c1^3 > {float(c3):.7f} > {float(cb)}^3, so c1 > {float(cb)}")
    gs = {}
    for key in keys:
        gs.setdefault(cs[key][1], []).append(cs[key][3])
    for g, ks in sorted(gs.items()):
        print(f"  G = {g}: cells {min(ks)}-{max(ks)} ({len(ks)} cells)")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        PATH = sys.argv[1]
    main()
