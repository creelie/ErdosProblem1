"""make_exceptional.py -- the figures of "An explicit exponent for the exceptional set of the
Erdos-Straus equation", written to exceptional/paper/figures/ as opaque 256-colour PNG files about
1600 px wide.

fig_simplex.png:  the triangle log u + log a + log alpha = log P in 3D, with the corners A, B, C, the
                  normal form that separates the classes in each, the strips next to them, the deep
                  strips down to the depth 1/10 and the inner centre.
fig_triples.png:  every Type II triple of the prime p = 4P - 1 = 342670648319, P = 2^8 3^2 5 7 11 13 17 19 23,
                  in the coordinates (log u, log a, log alpha)/log P, coloured by corner.
fig_ladder.png:   the constant c = (3/2) a^(1/3) for one, two and three corners, for the corners with
                  the strips under the three rules, for the full strips, with the deep strips (Monte
                  Carlo value of K from engine/anaxia_deep.json), with every pair down to the depth 1/10,
                  and for the whole triangle.
fig_strips.png:   the uncovered density U on the strips and the loss of each rule, by the quadrature of
                  checks/python/centre.py.
fig_engine.png:   the closed composites of the three engine runs, engine/anaxia_regions.json (128),
                  engine/anaxia_centre.json (37) and engine/anaxia_deep.json (15), by construction and
                  verdict.
fig_deep.png:     a splitting of a modulus into three factors for Maynard's second theorem: the groups of
                  one assignment of an example profile, and the ranges of the size n of the long variable
                  served by the assignments that cover [2/3 - s, 1/3], computed exactly.
fig_cert2.png:    the density k of the factorable pairs in the deep strips: the strip width 1/2 - 3v, the
                  Monte Carlo value of the engine run and the certified bound of every cell of
                  checks/c/certify_K.out, and the Monte Carlo share by depth and G.
fig_cert.png:     the certified bound for the loss: the bound of every cell of checks/c/certify_J.out
                  (written by certify_J with VERBOSE=1), the quadrature value of the loss and the bound
                  1/2 + v - 2 theta(v), and the bound and the quadrature value by ranges of the depth.
Requires numpy, matplotlib and Pillow.
"""
import itertools
import json
import math
import os
import sys
from math import gcd, isqrt

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "checks", "python"))
import centre as quad  # noqa: E402
import sum_K as sumk  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "exceptional", "paper", "figures")
KAPPA = 9 / (2 * math.pi ** 2)
COL = {"A": "#d95f02", "B": "#1b9e77", "C": "#7570b3", "O": "#bbbbbb", "S": "#e6ab02", "R": "#e0e0e0",
       "D": "#a6761d"}
LIGHT = {"A": "#f4b183", "B": "#8fd3bd", "C": "#bcb8de"}
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
    D = 2.5 / 42  # the strips are drawn 2.5 times wider than their true width 1/42
    D2 = 1 / 10   # the deep strips end at the depth 1/10

    def inner(k, d=D):  # the point at depth d where the bands of the two other coordinates meet
        v = np.full(3, 0.5 - d)
        v[k] = 2 * d
        return v

    # corner of the vertex e_k is the triangle (e_k, m(k,i), m(k,j)); its strip lies along m(i,k)-m(j,k)
    for k, i, j, name in ((2, 0, 1, "A"), (1, 0, 2, "B"), (0, 1, 2, "C")):
        tri(ax, [e[k], m(k, i), m(k, j)], COL[name], 0.85)
        tri(ax, [m(k, j), m(k, i), inner(j), inner(i)], LIGHT[name], 0.9)
        tri(ax, [inner(j), inner(i), inner(i, D2), inner(j, D2)], COL["D"], 0.5)
    tri(ax, [inner(0, D2), inner(1, D2), inner(2, D2)], COL["O"], 0.55)
    box = dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", alpha=0.95)
    ax.text2D(0.55, 0.90, "A: $ua\\leq X$, class $u/a$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.83, 0.40, "B: $u\\alpha\\leq X$\nclass $4u^2\\alpha$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.10, 0.30, "C: $a\\alpha\\leq X$\nclass $(4a^2\\alpha)^{-1}$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.33, 0.66, "strips: Maynard I,\nclass $u/a$", transform=ax.transAxes, fontsize=9, bbox=box, zorder=20)
    ax.text2D(0.02, 0.56, "deep strips to depth $1/10$:\nMaynard II, split moduli", transform=ax.transAxes, fontsize=9,
              bbox=box, zorder=20)
    ax.text2D(0.42, 0.45, "inner centre", transform=ax.transAxes, fontsize=9, bbox=box, zorder=20)
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
    pts = {k: [] for k in ("A", "B", "C", "S", "D", "O", "R")}
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
        elif max(u, a, al) ** 2 < P:
            top = math.log(max(u, a, al)) / L
            k = "S" if top > 0.5 - 1 / 42 else ("D" if top >= 0.4 else "O")
        else:
            k = "R"
        pts[k].append((math.log(u) / L, math.log(a) / L, math.log(al) / L))
    fig = plt.figure(figsize=(5.6, 4.6))
    ax = fig.add_subplot(projection="3d")
    names = {"A": "A", "B": "B", "C": "C", "S": "strips", "D": "deep strips", "O": "inner centre", "R": "other"}
    for k in ("R", "O", "D", "S", "A", "B", "C"):
        xs = np.array(pts[k])
        faint = k in ("O", "R")
        ax.scatter(xs[:, 0], xs[:, 1], xs[:, 2], s=0.6 if faint else 2.0, c=COL[k], depthshade=False,
                   label=f"{names[k]}: {len(xs)} triples", alpha=0.35 if faint else 0.9)
    ax.set_xlabel("$\\log u/\\log P$")
    ax.set_ylabel("$\\log a/\\log P$")
    ax.set_zlabel("$\\log\\alpha/\\log P$")
    ax.view_init(elev=24, azim=40)
    ax.legend(loc="upper left", markerscale=6, fontsize=9, frameon=False)
    ax.set_title(f"$p={p}$", fontsize=11)
    fig.tight_layout()
    save(fig, "fig_triples.png")
    return {k: len(v) for k, v in pts.items()}


def strip_constants():
    """The loss J and the constant c for the rules, by the quadrature of checks/python/centre.py."""
    quad.configure(6, 16)
    return {r: quad.J(r) for r in ("var", "prime", "chain")}


def deep_constants():
    """The certified bound for K (checks/c/certify_K.out through checks/python/sum_K.py) and the Monte
    Carlo value of K of the best accepted composite of engine/anaxia_deep.json."""
    cert = sumk.bound()
    eng = json.load(open(os.path.join(ROOT, "engine", "anaxia_deep.json")))["summary"]["best"]
    return float(cert["Kb"]), eng["K_estimate"]


def ladder_figure(J, Kcert, Kmc):
    I = 1 / 84 - 3 / 3528
    labels = ["one\ncorner", "two\ncorners", "three\ncorners", "strips\nR1", "strips\nR1-R2",
              "strips\nR1-R3", "full\nstrips", "deep strips,\ncertified", "deep strips,\nMonte Carlo",
              "to depth\n1/10", "whole\ntriangle"]
    base = 0.75 + 6 * (I - J["chain"])
    shares = [1 / 4, 1 / 2, 3 / 4, 0.75 + 6 * (I - J["var"]), 0.75 + 6 * (I - J["prime"]),
              base, 40 / 49, base + 6 * Kcert, base + 6 * Kmc, 24 / 25, 1]
    cs = [1.5 * (s * KAPPA / 4) ** (1 / 3) for s in shares]
    cols = [COL["A"], COL["B"], COL["C"], "#f2d27a", "#ecbe3c", COL["S"], COL["S"], COL["D"], COL["D"],
            COL["D"], COL["O"]]
    alphas = [0.95] * 6 + [0.45, 0.95, 0.6, 0.45, 0.45]
    fig = plt.figure(figsize=(10.0, 4.3))
    ax = fig.add_subplot(projection="3d")
    for i, (c, col) in enumerate(zip(cs, cols)):
        ax.bar3d(1.5 * i, 0, 0, 0.6, 0.15, c, color=col, alpha=alphas[i], shade=True, edgecolor="k")
        ax.text(1.5 * i + 0.3, 0.07, 0.75 + 0.13 * (i % 2), f"{c:.4f}", fontsize=7, zorder=30, ha="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.6"))
    ax.set_xticks([1.5 * i + 0.3 for i in range(len(cs))])
    ax.set_xticklabels([])
    for i, lab in enumerate(labels):  # two staggered rows, so that neighbouring labels do not meet
        ax.text(1.5 * i + 0.3, 0, -0.06 - 0.24 * (i % 2), lab, fontsize=7.5, ha="center", va="top")
    ax.tick_params(axis="z", labelsize=8)
    ax.set_yticks([])
    ax.set_ylim(0, 1.2)
    ax.set_zlabel("$c$ in $N\\exp(-c(\\log N)^{2/3})$")
    ax.set_zlim(0, 0.95)
    ax.view_init(elev=12, azim=-82)
    ax.set_box_aspect((3.6, 0.45, 1))
    fig.tight_layout()
    save(fig, "fig_ladder.png")
    return dict(zip([l.replace("\n", " ") for l in labels], [round(c, 6) for c in cs]))


def strips_figure():
    """Left: the uncovered density U of the rule (R1)-(R3) on the strip of depth vartheta.  Right: the
    loss l(vartheta) of the three rules and of (R1) with the smooth part, on a logarithmic scale."""
    quad.configure(6, 16)
    fig = plt.figure(figsize=(7.4, 3.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ths = np.concatenate([np.linspace(0.0015, 0.021, 14), np.linspace(0.0215, 1 / 42 - 0.0003, 12)])
    vs = np.linspace(0, 1, 41)
    TH, VV = np.meshgrid(ths, vs, indexing="ij")
    SS = 2 * TH + VV * (0.5 - 3 * TH)
    U = np.zeros_like(TH)
    for i, d in enumerate(ths):
        setup = quad.chain_setup(d)
        be, S = quad.tau(d), 0.5 + d
        for k in range(len(vs)):
            s = SS[i, k]
            U[i, k] = quad.g_chain(s, S - s, setup) if be < s < S - be else 0.0
    ax.plot_surface(TH, SS, U, cmap="YlOrBr", edgecolor="k", linewidth=0.15, rstride=1, cstride=1, alpha=0.95)
    ax.set_xlabel("depth $\\vartheta$", labelpad=4)
    ax.set_ylabel("$\\log u/\\log x$", labelpad=4)
    ax.set_zlabel("uncovered density $U$", labelpad=2)
    ax.tick_params(labelsize=7)
    ax.view_init(elev=24, azim=-140)
    ax.set_title("(a)", fontsize=10)
    ax2 = fig.add_subplot(1, 2, 2)
    th2 = np.linspace(0.001, 1 / 42 - 0.0002, 30)
    rules = [("var", "(R1)", "#e6c25a", "-"), ("prime", "(R1), (R2)", "#c98f00", "--"),
             ("smooth", "(R1), (R2), smooth part", "#8c6d1f", "-."), ("chain", "(R1), (R2), (R3)", COL["A"], "-")]
    for r, lab, col, ls in rules:
        ax2.semilogy(th2, [max(quad.loss(d, r), 1e-9) for d in th2], ls, color=col, lw=2.0, label=lab)
    ax2.axvline(1 / 54, color="0.5", lw=0.6, ls=":")
    ax2.axvline(1 / 44, color="0.5", lw=0.6, ls=":")
    ax2.text(1 / 54, 1.3, "$1/54$", fontsize=7.5, color="0.35", ha="center")
    ax2.text(1 / 44, 1.3, "$1/44$", fontsize=7.5, color="0.35", ha="center")
    ax2.set_xlabel("depth $\\vartheta$")
    ax2.set_ylabel("loss $\\ell(\\vartheta)$")
    ax2.set_xlim(0, 1 / 42)
    ax2.set_ylim(1e-6, 3)
    ax2.tick_params(labelsize=8)
    ax2.legend(loc="lower right", fontsize=7.5, frameon=False)
    ax2.set_title("(b)", fontsize=10)
    fig.tight_layout()
    save(fig, "fig_strips.png")
    return float(U.max())


def engine_figure():
    fig = plt.figure(figsize=(11.6, 4.8))
    vcol = {"accepted": "#1b9e77", "distinctness": "#d95f02", "other": "#7570b3"}
    # (a) the run behind the three corners
    data = json.load(open(os.path.join(ROOT, "engine", "anaxia_regions.json")))
    grid = {}
    for r in data["composites"]:
        meas = r["measurement"]
        n = len(meas["corners"])
        centre = 1 if meas["centre"] == "eh" else 0
        verdict = "accepted" if r["outcome"] == "accepted" else (
            "distinctness" if "normal form" in r["reason"] else "other")
        grid[(n, centre, verdict)] = grid.get((n, centre, verdict), 0) + 1
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    elev, azim = 26, -38
    view = (math.cos(math.radians(azim)), math.sin(math.radians(azim)))
    ax.computed_zorder = False
    for (n, centre, v), cnt in sorted(grid.items(), key=lambda kv: (
            view[0] * kv[0][0] / 3.4 + view[1] * kv[0][1] / 1.3, list(vcol).index(kv[0][2]))):
        z0 = sum(c for (n2, c2, v2), c in grid.items() if n2 == n and c2 == centre and
                 list(vcol).index(v2) < list(vcol).index(v))
        ax.bar3d(n, centre, z0, 0.4, 0.3, cnt, color=vcol[v], edgecolor="k", alpha=1.0)
    for (n, centre) in sorted({(k[0], k[1]) for k in grid if k[1] == 0}):
        tot = sum(c for k, c in grid.items() if k[:2] == (n, centre))
        acc = grid.get((n, centre, "accepted"), 0)
        ax.text(n + 0.05, centre, tot + 1.5, f"{tot}" + (f" ({acc} acc.)" if acc else ""), fontsize=7.5,
                zorder=50, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.6"))
    for v, c in vcol.items():
        ax.bar3d(0, 0, 0, 0, 0, 0, color=c, label={"accepted": "accepted",
                 "distinctness": "form does not separate the classes",
                 "other": "a centre class equals a corner class, or nothing sieved"}[v])
    ax.set_xticks([0.2, 1.2, 2.2, 3.2])
    ax.set_xticklabels(["0", "1", "2", "3"], fontsize=7)
    ax.set_xlabel("corners used", fontsize=8)
    ax.set_yticks([0.15, 1.15])
    ax.set_yticklabels(["no centre", "centre"], fontsize=7)
    ax.set_zlabel("constructions", fontsize=8)
    ax.tick_params(axis="z", labelsize=7)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), fontsize=7, frameon=False, ncol=1)
    ax.view_init(elev=elev, azim=azim)
    ax.set_title("(a) corners: 128 constructions", fontsize=9)
    # (b) the run behind the strips
    data = json.load(open(os.path.join(ROOT, "engine", "anaxia_centre.json")))
    xs = ["none", "var", "prime", "smooth", "chain", "steps", "eh"]
    ys = ["farey", "M", "I"]
    wcol = {"accepted": "#1b9e77", "conditional": "#66a6d9", "form": "#d95f02", "coincidence": "#7570b3",
            "rule": "#b2182b"}
    cells = {}
    for r in data["composites"]:
        meas = r["measurement"]
        if meas["centre"] == "none":
            x, y = 0, 0
        else:
            x, y = xs.index(meas["rule"]), ys.index(meas["form"])
        if r["outcome"] in ("accepted", "conditional"):
            v = r["outcome"]
        elif "form" in r["reason"]:
            v = "form"
        elif "coincidences" in r["reason"]:
            v = "coincidence"
        else:
            v = "rule"
        cells.setdefault((x, y), []).append(v)
    ax = fig.add_subplot(1, 3, 2, projection="3d")
    ax.computed_zorder = False
    order = list(wcol)
    for (x, y), vv in sorted(cells.items(), key=lambda kv: (-kv[0][1], kv[0][0])):
        z0 = 0
        for v in sorted(vv, key=order.index):
            ax.bar3d(x, y, z0, 0.5, 0.5, 1, color=wcol[v], edgecolor="k", alpha=1.0)
            z0 += 1
    for v, c in wcol.items():
        ax.bar3d(0, 0, 0, 0, 0, 0, color=c, label={"accepted": "accepted", "conditional": "conditional (EH)",
                 "form": "form does not separate the centre", "coincidence": "coincidences not allowed",
                 "rule": "rule claims a divisor that is absent"}[v])
    ax.set_xticks([i + 0.25 for i in range(len(xs))])
    ax.set_xticklabels(["left\nout", "R1", "R1-2", "smooth", "R1-3", "steps", "EH"], fontsize=6.5)
    ax.set_yticks([i + 0.25 for i in range(len(ys))])
    ax.set_yticklabels(["$u/a$", "$4u^2\\alpha$", "$(4a^2\\alpha)^{-1}$"], fontsize=6.5)
    ax.set_zlabel("constructions", fontsize=8)
    ax.set_zticks([0, 1, 2])
    ax.tick_params(axis="z", labelsize=7)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), fontsize=7, frameon=False, ncol=2)
    ax.view_init(elev=30, azim=-58)
    ax.set_title("(b) strips: 37 constructions", fontsize=9)
    # (c) the run behind the deep strips
    data = json.load(open(os.path.join(ROOT, "engine", "anaxia_deep.json")))
    xs3 = ["none", "window", "eh", "two", "one", "noreserve"]
    gl = sorted({tuple(r["measurement"]["G"]) for r in data["composites"] if r["measurement"]["deep"] == "split"},
                key=lambda g: (len(g), -g[0]))
    dcol = {"accepted": "#1b9e77", "conditional": "#66a6d9", "window": "#7570b3", "fill": "#d95f02"}
    deep = {}
    for r in data["composites"]:
        meas = r["measurement"]
        if meas["deep"] == "split":
            x, y = xs3.index(meas["rule"]), gl.index(tuple(meas["G"]))
        else:
            x, y = xs3.index(meas["deep"]), 0
        if r["outcome"] in ("accepted", "conditional"):
            v = r["outcome"]
        else:
            v = "window" if meas["deep"] == "window" else "fill"
        deep.setdefault((x, y), []).append(v)
    ax = fig.add_subplot(1, 3, 3, projection="3d")
    ax.computed_zorder = False
    order = list(dcol)
    for (x, y), vv in sorted(deep.items(), key=lambda kv: (-kv[0][1], kv[0][0])):
        z0 = 0
        for v in sorted(vv, key=order.index):
            ax.bar3d(x, y, z0, 0.5, 0.5, 1, color=dcol[v], edgecolor="k", alpha=1.0)
            z0 += 1
    for v, c in dcol.items():
        ax.bar3d(0, 0, 0, 0, 0, 0, color=c, label={"accepted": "accepted", "conditional": "conditional (EH)",
                 "window": "window empty past $1/42$", "fill": "fill misses the targets"}[v])
    ax.set_xticks([i + 0.25 for i in range(len(xs3))])
    ax.set_xticklabels(["left\nout", "win-\ndow", "EH", "two\ngroups", "one\ngroup", "no\nreserve"], fontsize=6)
    ax.set_box_aspect((1.5, 1, 0.9))
    ax.set_yticks([i + 0.25 for i in range(len(gl))])
    ax.set_yticklabels([f"$G={g[0]:g}$" if len(g) == 1 else "best $G$" for g in gl], fontsize=6.5)
    ax.set_zlabel("constructions", fontsize=8)
    ax.set_zticks([0, 1])
    ax.tick_params(axis="z", labelsize=7)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), fontsize=7, frameon=False, ncol=2)
    ax.view_init(elev=30, azim=-58)
    ax.set_title(f"(c) deep strips: {len(data['composites'])} constructions", fontsize=9)
    fig.tight_layout()
    save(fig, "fig_engine.png")
    return grid, {k: len(v) for k, v in cells.items()}, {k: len(v) for k, v in deep.items()}


def served_intervals(Z, t, s, G, step=0.0025):
    """Every assignment of the example profile (Z, t) of size s, with targets on a grid of the given step,
    and the open interval of n that it serves: the four conditions of Maynard's second theorem at every
    vertex of the box of sizes, with d_L = s minus the others.  Exact for the given targets."""
    out = []
    for lab in itertools.product(range(3), repeat=len(Z)):
        b = [0.0, 0.0, 0.0]
        for z, g in zip(Z, lab):
            b[g] += z
        for L in range(3):
            others = [i for i in range(3) if i != L]
            for F in ([], [others[0]], [others[1]], others):
                res = max(len(F) - 1, 0) * G
                grids = [[b[i] + step * k for k in range(int((t - res) / step) + 1)] for i in F]
                for xs in (itertools.product(*grids) if F else [()]):
                    x = dict(zip(F, xs))
                    if sum(x[i] - b[i] for i in F) + res >= t:
                        continue
                    lo, hi = -1.0, 2.0
                    for corner in (itertools.product(*[(x[i], x[i] + G) for i in F]) if F else [()]):
                        d = list(b)
                        for i, v in zip(F, corner):
                            d[i] = v
                        d[L] = s - sum(d[i] for i in range(3) if i != L)
                        d1, d2, d3 = d
                        lo = max(lo, d1)
                        hi = min(hi, (1 - d2 - 2 * d3) / 2, (2 - d1 - 4 * d2 - 3 * d3) / 2, 2 - d1 - 5 * d2 - 2 * d3)
                    if hi > lo:
                        out.append(dict(lo=lo, hi=hi, lab=lab, L=L, F=tuple(F), x=xs, b=tuple(b)))
    return out


def greedy_cover(iv, a, b):
    reach, chosen = a, []
    while reach < b:
        cand = [r for r in iv if r["lo"] < reach and r["hi"] > reach]
        if not cand:
            return None
        best = max(cand, key=lambda r: r["hi"])
        chosen.append(best)
        reach = best["hi"]
    return chosen


def deep_figure():
    """An example profile of size s = 0.55 (depth 0.05) with G = 0.03: four prime factors above x^G and
    small primes of total size 0.115 (dust 0.002 included).  Left: the groups of one assignment that fills
    two groups; right: assignments that together serve every n in [2/3 - s, 1/3]."""
    s, G = 0.55, 0.03
    Z = [0.20, 0.12, 0.07, 0.045]
    small = [0.002, 0.005, 0.007, 0.009, 0.012, 0.015, 0.018, 0.022, 0.025]
    t = s - sum(Z)
    assert abs(t - sum(small)) < 1e-12
    iv = served_intervals(Z, t, s, G)
    pick = next(r for r in iv if r["lab"] == (0, 1, 1, 2) and r["L"] == 2 and r["F"] == (0, 1)
                and abs(r["x"][0] - 0.2525) < 1e-9 and abs(r["x"][1] - 0.21) < 1e-9)
    lo_n, hi_n = 2 / 3 - s, 1 / 3
    left = greedy_cover(iv, lo_n, pick["lo"] + 1e-9)
    right = greedy_cover(iv, pick["hi"] - 1e-9, hi_n)
    chosen = left + [pick] + right
    # the fill of the picked assignment: items in increasing order into groups 1 and 2
    groups = [[("big", z)] for z in (0.20, 0.19)] + [[("big", 0.045)]]
    groups[1] = [("big", 0.12), ("big", 0.07)]
    sizes = [0.20, 0.19, 0.045]
    k = 0
    for gi, target in ((0, pick["x"][0]), (1, pick["x"][1])):
        while sizes[gi] < target:
            groups[gi].append(("small", small[k]))
            sizes[gi] += small[k]
            k += 1
    for z in small[k:]:
        groups[2].append(("small", z))
        sizes[2] += z
    fig = plt.figure(figsize=(9.6, 4.4))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    gcol = [COL["A"], COL["B"], COL["C"]]
    for gi, items in enumerate(groups):
        z0 = 0.0
        for kind, z in items:
            col = gcol[gi] if kind == "big" else "#f3e3b0"
            ax.bar3d(1.6 * gi, 0, z0, 0.8, 0.8, z, color=col, edgecolor="k", linewidth=0.4, alpha=1.0, shade=True)
            z0 += z
        ax.text(1.6 * gi + 0.1, 0.4, z0 + 0.025, f"$d_{gi + 1}={z0:.3f}$", fontsize=8, zorder=40,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.6"))
    for gi in (0, 1):
        x = pick["x"][gi]
        for h, ls in ((x, "-"), (x + G, "--")):
            xs = [1.6 * gi - 0.15, 1.6 * gi + 0.95, 1.6 * gi + 0.95, 1.6 * gi - 0.15, 1.6 * gi - 0.15]
            ys = [-0.15, -0.15, 0.95, 0.95, -0.15]
            ax.plot(xs, ys, [h] * 5, ls, color="k", lw=0.8, zorder=30)
        ax.text(1.6 * gi + 1.0, 0.95, x, f"$x_{gi + 1}$", fontsize=8, zorder=40)
    ax.set_xticks([0.4, 2.0, 3.6])
    ax.set_xticklabels(["$q_1$ (filled)", "$q_2$ (filled)", "$q_3$ (rest)"], fontsize=8)
    ax.set_yticks([])
    ax.set_zlabel("size (exponent of $x$)", fontsize=8)
    ax.tick_params(axis="z", labelsize=7)
    ax.set_zlim(0, 0.32)
    ax.view_init(elev=18, azim=-62)
    ax.set_box_aspect((1.6, 0.6, 1.2))
    ax.set_title("(a) one assignment: prime factors above $x^G$ (dark), below $x^G$ (light)", fontsize=8.5)
    bx = fig.add_subplot(1, 2, 2, projection="3d")
    for j, r in enumerate(chosen):
        col = COL["S"] if r is pick else "#9ecae1"
        bx.bar3d(r["lo"], j, 0, r["hi"] - r["lo"], 0.55, 0.25, color=col, edgecolor="k", linewidth=0.4, alpha=1.0)
    for v in (lo_n, hi_n):
        bx.plot([v, v], [-0.3, len(chosen)], [0, 0], color="k", lw=1.2)
        bx.plot([v, v], [len(chosen), len(chosen)], [0, 0.6], color="k", lw=1.2)
    bx.set_xlabel("size $n$ of the long variable", fontsize=8)
    bx.set_yticks([j + 0.3 for j in range(len(chosen))])
    bx.set_yticklabels([f"{j + 1}" for j in range(len(chosen))], fontsize=7)
    bx.set_ylabel("assignment", fontsize=8)
    bx.set_zticks([])
    bx.set_xlim(0.1, 0.36)
    bx.set_xticks([lo_n, 0.2, 0.25, 0.3, hi_n])
    bx.set_xticklabels(["$2/3-s$", "0.2", "0.25", "0.3", "$1/3$"], fontsize=7)
    bx.set_zlim(0, 0.6)
    bx.view_init(elev=48, azim=-72)
    bx.set_box_aspect((2.0, 1.4, 0.35))
    bx.set_title("(b) the ranges of $n$ served by six assignments", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig_deep.png")
    return dict(chosen=[(round(r["lo"], 4), round(r["hi"], 4), r["lab"], r["L"], r["F"]) for r in chosen],
                sizes=[round(z, 4) for z in sizes])


def cert2_figure():
    """Left: on the depths [1/42, 1/10], the width 1/2 - 3v of the strip, the Monte Carlo density k of the
    best accepted composite of the engine run and the certified bound of every cell.  Right: the Monte
    Carlo share of the factorable pairs (rule (two)) by depth and G."""
    eng = json.load(open(os.path.join(ROOT, "engine", "anaxia_deep.json")))["summary"]
    gbest = eng["best"]["G"]
    mc = {}
    for key, v in eng["mc"].items():
        d, g, m = key.split(":")
        if m == "0":
            mc[(float(d), float(g))] = v
    depths = [d for d in eng["depths"] if d <= 0.1 + 1e-12]
    kmc = [max(mc[(d, g)] for g in gbest) * (0.5 - 3 * d) for d in depths]
    cert = sumk.bound()
    fig = plt.figure(figsize=(10.4, 4.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")

    def ribbon(axis, xs, zs, y, col, alpha):
        for x0, x1, z0, z1 in zip(xs[:-1], xs[1:], zs[:-1], zs[1:]):
            axis.add_collection3d(Poly3DCollection([[(x0, y, 0), (x1, y, 0), (x1, y, z1), (x0, y, z0)]],
                                                   facecolors=col, edgecolors="none", alpha=alpha))
        axis.plot(xs, [y] * len(xs), zs, color="k", lw=0.6)

    xs = np.linspace(1 / 42, 0.1, 40)
    sx, sz = [], []
    for (a, b) in cert["keys"]:
        lo, hi = max(float(a), 1 / 42), min(float(b), 0.1)
        if hi > lo:
            sx += [lo, hi]
            sz += [float(cert["cells"][(a, b)][0])] * 2
    if sx:
        ribbon(ax, 1000 * np.array(sx), np.array(sz), 0, COL["D"], 0.9)
    ribbon(ax, 1000 * np.array(depths), np.array(kmc), 1, COL["B"], 0.85)
    ribbon(ax, 1000 * xs, 0.5 - 3 * xs, 2, COL["O"], 0.6)
    ax.set_xlabel("depth $10^3\\vartheta$", labelpad=8)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["cell bound", "Monte Carlo", "$1/2-3\\vartheta$"], fontsize=7)
    ax.set_zlabel("density $k(\\vartheta)$")
    ax.set_xlim(1000 / 42, 100)
    ax.set_ylim(-0.3, 2.3)
    ax.set_zlim(0, 0.45)
    ax.view_init(elev=18, azim=-62)
    ax.set_box_aspect((2.0, 1.0, 1))
    ax.set_title(f"(a) the density $k$: $\\mathcal{{K}}>{float(cert['Kb']):.5f}$ certified", fontsize=10)
    bx = fig.add_subplot(1, 2, 2, projection="3d")
    gs = sorted({g for (_, g) in mc}, reverse=True)
    gcol = plt.cm.viridis(np.linspace(0.85, 0.1, len(gs)))
    for j, g in enumerate(gs):
        dd = [d for d in depths if (d, g) in mc]
        ribbon(bx, 1000 * np.array(dd), np.array([mc[(d, g)] for d in dd]), j, gcol[j], 0.85)
    bx.set_xlabel("depth $10^3\\vartheta$", labelpad=8)
    bx.set_yticks(list(range(len(gs))))
    bx.set_yticklabels([f"$G={g:g}$" for g in gs], fontsize=7)
    bx.set_zlabel("share of the strip")
    bx.set_xlim(1000 / 42, 100)
    bx.set_ylim(-0.3, len(gs) - 0.7)
    bx.set_zlim(0, 0.8)
    bx.view_init(elev=18, azim=-62)
    bx.set_box_aspect((2.0, 1.2, 1))
    bx.set_title("(b) Monte Carlo share by depth and $G$", fontsize=10)
    fig.tight_layout()
    save(fig, "fig_cert2.png")
    return dict(Kb=float(cert["Kb"]), cb=float(cert["cb"]), K_mc=eng["best"]["K_estimate"], G=gbest)


def read_cert():
    """The cells and the parts of checks/c/certify_J.out."""
    cells, parts = [], []
    with open(os.path.join(ROOT, "checks", "c", "certify_J.out")) as f:
        for line in f:
            if line.startswith("cell ["):
                a, b = (float(x) for x in line[6:line.index("]")].split(","))
                lbar = float(line.split("lbar=")[1].split()[0])
                cells.append((a, b, lbar))
            elif " cells " in line:
                parts.append(float(line.split()[-1]))
    return cells, parts


def cert_figure():
    cells, parts = read_cert()
    theta = lambda v: min(0.1 - 1.4 * v, 0.5 - 19 * v)
    lo, hi = 0.005, 1 / 42 - 2e-4
    quad.configure(6, 16)
    vs = np.concatenate([np.linspace(lo, 1 / 44, 26), 1 / 42 - (1 / 42 - 1 / 44) * 0.82 ** np.arange(1, 17)])
    vs = vs[vs <= hi]
    ls = np.array([quad.loss(v, "chain") for v in vs])
    r1 = np.array([0.5 + v - 2 * theta(v) for v in vs])
    fig = plt.figure(figsize=(10.4, 4.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")

    def ribbon(xs, zs, y, col, alpha):
        for x0, x1, z0, z1 in zip(xs[:-1], xs[1:], zs[:-1], zs[1:]):
            ax.add_collection3d(Poly3DCollection([[(x0, y, -4), (x1, y, -4), (x1, y, z1), (x0, y, z0)]],
                                                 facecolors=col, edgecolors="none", alpha=alpha))
        ax.plot(xs, [y] * len(xs), zs, color="k", lw=0.6)

    sx, sz = [], []
    for a, b, lbar in cells:
        if a >= lo and b <= hi:
            sx += [a, b]
            sz += [math.log10(lbar)] * 2
    sx, sz = 1000 * np.array(sx), np.array(sz)
    ribbon(1000 * vs, np.log10(ls), 0, COL["B"], 0.8)
    ribbon(sx, sz, 1, COL["A"], 0.75)
    ribbon(1000 * vs, np.log10(r1), 2, COL["O"], 0.6)
    ax.set_xlabel("depth $10^3\\vartheta$", labelpad=9)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["quadrature", "cell bound", "$1/2+\\vartheta-2\\theta$"], fontsize=7)
    ax.set_zticks([-3, -2, -1, 0])
    ax.set_zticklabels(["$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "1"])
    ax.set_zlim(-4, 0)
    ax.set_xlim(1000 * lo, 1000 / 42)
    ax.set_ylim(-0.3, 2.3)
    ax.set_zlabel("loss $\\ell(\\vartheta)$")
    ax.view_init(elev=13, azim=-72)
    ax.set_box_aspect((2.0, 1.0, 1))
    ax.set_title("(a) bounds for the loss", fontsize=10)

    bx = fig.add_subplot(1, 2, 2, projection="3d")
    names = ["$[0,0.005]$", "$[0.005,1/54]$", "$[1/54,1/44]$", "$[1/44,1/42]$"]
    bound = [parts[0] + parts[1], parts[2], parts[3], parts[4] + parts[5]]
    cuts = [(0, 0.005), (0.005, 1 / 54), (1 / 54, 1 / 44), (1 / 44, 1 / 42)]
    qv = []
    for a, b in cuts:
        dd, dw = quad.nodes(a, b, [])
        qv.append(sum(wi * quad.loss(di, "chain") for di, wi in zip(dd, dw)))
    for i in range(4):
        for j, (val, col) in enumerate(((qv[i], COL["B"]), (bound[i], COL["A"]))):
            bx.bar3d(i, 0.6 * j, 0, 0.5, 0.4, 1e4 * val, color=col, alpha=0.9, shade=True, edgecolor="k")
        bx.text(i - 0.1, 1.2, 1e4 * bound[i] + 0.15, f"{1e4 * bound[i]:.2f}", fontsize=7, zorder=30,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.6"))
    bx.set_xticks([i + 0.25 for i in range(4)])
    bx.set_xticklabels(names, fontsize=7)
    bx.set_yticks([0.2, 0.8])
    bx.set_yticklabels(["quadrature", "bound"], fontsize=7)
    bx.set_zlabel("$10^4\\times\\int\\ell$")
    bx.view_init(elev=20, azim=-58)
    bx.set_box_aspect((1.8, 0.9, 1))
    bx.set_title("(b) the integral by ranges", fontsize=10)
    fig.tight_layout()
    save(fig, "fig_cert.png")
    return {"bound": [round(x, 9) for x in bound], "quadrature": [round(x, 9) for x in qv]}


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    simplex_figure()
    print("triples by corner:", triples_figure())
    J = strip_constants()
    Kcert, Kmc = deep_constants()
    print("constants:", ladder_figure(J, Kcert, Kmc))
    print("largest U on the grid:", strips_figure())
    print("engine grid:", engine_figure())
    print("deep strips example:", deep_figure())
    print("deep strips bound:", cert2_figure())
    print("certified bound by ranges:", cert_figure())
