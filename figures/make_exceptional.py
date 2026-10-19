"""make_exceptional.py -- the figures of "An explicit exponent for the exceptional set of the
Erdos-Straus equation", written to exceptional/paper/figures/ as opaque 256-colour PNG files about
1600 px wide.

fig_simplex.png:  the triangle log u + log a + log alpha = log P in 3D, with the corners A, B, C used by
                  the sieve, the normal form that separates the classes in each, and the centre.
fig_triples.png:  every Type II triple of the prime p = 4P - 1 = 342670648319, P = 2^8 3^2 5 7 11 13 17 19 23,
                  in the coordinates (log u, log a, log alpha)/log P, coloured by corner.
fig_ladder.png:   the constant c = (3/2) a^(1/3) for one, two and three corners and for the whole triangle.
fig_engine.png:   the 128 closed composites of the engine run in engine/anaxia_regions.json, by the
                  number of corners, the centre and the verdict.
Requires numpy, matplotlib and Pillow.
"""
import json
import math
import os
from math import gcd, isqrt

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "exceptional", "paper", "figures")
KAPPA = 9 / (2 * math.pi ** 2)
COL = {"A": "#d95f02", "B": "#1b9e77", "C": "#7570b3", "O": "#bbbbbb"}
plt.rcParams.update({"font.size": 10, "font.family": "serif"})


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=440, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    im = Image.open(path).convert("RGB")
    ink = np.asarray(im).min(axis=2) < 235
    rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
    m = 24
    im = im.crop((max(cols[0] - m, 0), max(rows[0] - m, 0),
                  min(cols[-1] + m + 1, im.width), min(rows[-1] + m + 1, im.height)))
    if im.width > 1600:
        im = im.resize((1600, round(im.height * 1600 / im.width)), Image.LANCZOS)
    im.quantize(colors=256, method=Image.Quantize.MEDIANCUT).save(path, optimize=True)


def tri(ax, pts, color, alpha):
    ax.add_collection3d(Poly3DCollection([pts], facecolors=color, edgecolors="k", linewidths=0.6, alpha=alpha))


def simplex_figure():
    fig = plt.figure(figsize=(5.6, 4.6))
    ax = fig.add_subplot(projection="3d")
    e = np.eye(3)
    m = lambda i, j: (e[i] + e[j]) / 2
    # corner of the vertex e_k is the triangle (e_k, m(k,i), m(k,j))
    tri(ax, [e[2], m(2, 0), m(2, 1)], COL["A"], 0.85)      # alpha large: u a <= X
    tri(ax, [e[1], m(1, 0), m(1, 2)], COL["B"], 0.85)      # a large: u alpha <= X
    tri(ax, [e[0], m(0, 1), m(0, 2)], COL["C"], 0.85)      # u large: a alpha <= X
    tri(ax, [m(0, 1), m(1, 2), m(0, 2)], COL["O"], 0.55)
    box = dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", alpha=0.95)
    ax.text2D(0.55, 0.90, "A: $ua\\leq X$, class $u/a$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.83, 0.40, "B: $u\\alpha\\leq X$\nclass $4u^2\\alpha$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.12, 0.30, "C: $a\\alpha\\leq X$\nclass $(4a^2\\alpha)^{-1}$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.40, 0.52, "centre: moduli\nabove $\\sqrt{p}$", transform=ax.transAxes, fontsize=9, bbox=box, zorder=20)
    ax.set_xlabel("$\\log u/\\log P$")
    ax.set_ylabel("$\\log a/\\log P$")
    ax.set_zlabel("$\\log\\alpha/\\log P$")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.view_init(elev=24, azim=40)
    ax.set_box_aspect((1, 1, 0.9))
    fig.tight_layout()
    save(fig, "fig_simplex.png")


def triples_figure():
    f = [(2, 8), (3, 2), (5, 1), (7, 1), (11, 1), (13, 1), (17, 1), (19, 1), (23, 1)]
    P = 1
    for q, k in f:
        P *= q ** k
    p = 4 * P - 1
    X = math.sqrt(p) / 16
    L = math.log(P)
    pts = {k: [] for k in COL}
    # coprime (u, a): at each prime the exponent goes to u, to a, or to neither; the rest to alpha
    choices = []
    for q, k in f:
        opts = [(1, 1)]
        for i in range(1, k + 1):
            opts += [(q ** i, 1), (1, q ** i)]
        choices.append(opts)
    for combo in np.ndindex(*[len(c) for c in choices]):
        u = a = 1
        for c, idx in zip(choices, combo):
            u *= c[idx][0]
            a *= c[idx][1]
        al = P // (u * a)
        if u >= 2 and a >= 2 and u * a <= X:
            k = "A"
        elif u * al <= X:
            k = "B"
        elif a * al <= X:
            k = "C"
        else:
            k = "O"
        pts[k].append((math.log(u) / L, math.log(a) / L, math.log(al) / L))
    fig = plt.figure(figsize=(5.6, 4.6))
    ax = fig.add_subplot(projection="3d")
    names = {"A": "A", "B": "B", "C": "C", "O": "centre"}
    for k in ("O", "A", "B", "C"):
        xs = np.array(pts[k])
        ax.scatter(xs[:, 0], xs[:, 1], xs[:, 2], s=0.6 if k == "O" else 2.0, c=COL[k], depthshade=False,
                   label=f"{names[k]}: {len(xs)} triples", alpha=0.35 if k == "O" else 0.9)
    ax.set_xlabel("$\\log u/\\log P$")
    ax.set_ylabel("$\\log a/\\log P$")
    ax.set_zlabel("$\\log\\alpha/\\log P$")
    ax.view_init(elev=24, azim=40)
    ax.legend(loc="upper left", markerscale=6, fontsize=9, frameon=False)
    ax.set_title(f"$p={p}$", fontsize=11)
    fig.tight_layout()
    save(fig, "fig_triples.png")
    return {k: len(v) for k, v in pts.items()}


def ladder_figure():
    labels = ["one\ncorner", "two\ncorners", "three\ncorners", "whole\ntriangle"]
    fr = [1, 2, 3, 4]
    cs = [1.5 * (j * KAPPA / 16) ** (1 / 3) for j in fr]
    fig = plt.figure(figsize=(5.6, 3.9))
    ax = fig.add_subplot(projection="3d")
    cols = [COL["A"], COL["B"], COL["C"], COL["O"]]
    for i, (c, col) in enumerate(zip(cs, cols)):
        ax.bar3d(i, 0, 0, 0.5, 0.15, c, color=col, alpha=0.95 if i < 3 else 0.45, shade=True, edgecolor="k")
        ax.text(i + 0.02, 0.07, c + 0.03, f"{c:.4f}", fontsize=10, zorder=30,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6"))
    ax.set_xticks([i + 0.25 for i in range(4)])
    ax.set_xticklabels(labels, fontsize=9)
    ax.tick_params(axis="z", labelsize=8)
    ax.set_yticks([])
    ax.set_ylim(0, 1.2)
    ax.set_zlabel("$c$ in $N\\exp(-c(\\log N)^{2/3})$")
    ax.set_zlim(0, 0.85)
    ax.view_init(elev=12, azim=-80)
    fig.tight_layout()
    save(fig, "fig_ladder.png")


def engine_figure():
    data = json.load(open(os.path.join(ROOT, "engine", "anaxia_regions.json")))
    grid = {}
    for r in data["composites"]:
        meas = r["measurement"]
        n = len(meas["corners"])
        centre = 1 if meas["centre"] == "eh" else 0
        verdict = "accepted" if r["outcome"] == "accepted" else (
            "distinctness" if "normal form" in r["reason"] else "other")
        grid[(n, centre, verdict)] = grid.get((n, centre, verdict), 0) + 1
    fig = plt.figure(figsize=(5.6, 4.3))
    ax = fig.add_subplot(projection="3d")
    vcol = {"accepted": "#1b9e77", "distinctness": "#d95f02", "other": "#7570b3"}
    elev, azim = 26, -38
    view = (math.cos(math.radians(azim)), math.sin(math.radians(azim)))
    ax.computed_zorder = False
    # bars far from the viewer first, so that nearer bars are drawn over them
    for (n, centre, v), cnt in sorted(grid.items(), key=lambda kv: (
            view[0] * kv[0][0] / 3.4 + view[1] * kv[0][1] / 1.3, list(vcol).index(kv[0][2]))):
        z0 = sum(c for (n2, c2, v2), c in grid.items() if n2 == n and c2 == centre and
                 list(vcol).index(v2) < list(vcol).index(v))
        ax.bar3d(n, centre, z0, 0.4, 0.3, cnt, color=vcol[v], edgecolor="k", alpha=1.0)
    for (n, centre) in sorted({(k[0], k[1]) for k in grid if k[1] == 0}):
        tot = sum(c for k, c in grid.items() if k[:2] == (n, centre))
        acc = grid.get((n, centre, "accepted"), 0)
        ax.text(n + 0.05, centre, tot + 1.5, f"{tot}" + (f" ({acc} acc.)" if acc else ""), fontsize=8.5,
                zorder=50, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.6"))
    for v, c in vcol.items():
        ax.bar3d(0, 0, 0, 0, 0, 0, color=c, label={"accepted": "accepted",
                 "distinctness": "rejected: form does not separate the classes",
                 "other": "rejected: coincidences in the centre, or nothing sieved"}[v])
    ax.set_xticks([0.2, 1.2, 2.2, 3.2])
    ax.set_xticklabels(["0", "1", "2", "3"])
    ax.set_xlabel("corners used")
    ax.set_yticks([0.15, 1.15])
    ax.set_yticklabels(["no centre", "centre"])
    ax.set_zlabel("composites")
    ax.legend(loc="upper left", fontsize=8.5, frameon=False)
    ax.view_init(elev=elev, azim=azim)
    fig.tight_layout()
    save(fig, "fig_engine.png")
    return grid


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    simplex_figure()
    print("triples by corner:", triples_figure())
    ladder_figure()
    print("engine grid:", engine_figure())
