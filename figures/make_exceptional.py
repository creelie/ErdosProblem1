"""make_exceptional.py -- the figures of "An explicit exponent for the exceptional set of the
Erdos-Straus equation", written to exceptional/paper/figures/ as opaque 256-colour PNG files about
1600 px wide.

fig_simplex.png:  the triangle log u + log a + log alpha = log P in 3D, with the corners A, B, C, the
                  normal form that separates the classes in each, the strips next to them and the inner
                  centre.
fig_triples.png:  every Type II triple of the prime p = 4P - 1 = 342670648319, P = 2^8 3^2 5 7 11 13 17 19 23,
                  in the coordinates (log u, log a, log alpha)/log P, coloured by corner.
fig_ladder.png:   the constant c = (3/2) a^(1/3) for one, two and three corners, for the corners with
                  the strips under the three rules, for the full strips and for the whole triangle.
fig_strips.png:   the uncovered density U on the strips and the loss of each rule, by the quadrature of
                  checks/python/centre.py.
fig_engine.png:   the closed composites of the two engine runs, engine/anaxia_regions.json (128) and
                  engine/anaxia_centre.json (37), by construction and verdict.
Requires numpy, matplotlib and Pillow.
"""
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "exceptional", "paper", "figures")
KAPPA = 9 / (2 * math.pi ** 2)
COL = {"A": "#d95f02", "B": "#1b9e77", "C": "#7570b3", "O": "#bbbbbb", "S": "#e6ab02", "R": "#e0e0e0"}
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

    def inner(k):  # the point where the strips of the two other coordinates meet the inner triangle
        v = np.full(3, 0.5 - D)
        v[k] = 2 * D
        return v

    # corner of the vertex e_k is the triangle (e_k, m(k,i), m(k,j)); its strip lies along m(i,k)-m(j,k)
    for k, i, j, name in ((2, 0, 1, "A"), (1, 0, 2, "B"), (0, 1, 2, "C")):
        tri(ax, [e[k], m(k, i), m(k, j)], COL[name], 0.85)
        tri(ax, [m(k, j), m(k, i), inner(j), inner(i)], LIGHT[name], 0.9)
    tri(ax, [inner(0), inner(1), inner(2)], COL["O"], 0.55)
    box = dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", alpha=0.95)
    ax.text2D(0.55, 0.90, "A: $ua\\leq X$, class $u/a$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.83, 0.40, "B: $u\\alpha\\leq X$\nclass $4u^2\\alpha$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.10, 0.30, "C: $a\\alpha\\leq X$\nclass $(4a^2\\alpha)^{-1}$", transform=ax.transAxes, fontsize=10, bbox=box, zorder=20)
    ax.text2D(0.36, 0.60, "strips: Maynard,\nclass $u/a$", transform=ax.transAxes, fontsize=9, bbox=box, zorder=20)
    ax.text2D(0.40, 0.44, "inner centre", transform=ax.transAxes, fontsize=9, bbox=box, zorder=20)
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
    pts = {k: [] for k in ("A", "B", "C", "S", "O", "R")}
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
            k = "S" if math.log(max(u, a, al)) / L > 0.5 - 1 / 42 else "O"
        else:
            k = "R"
        pts[k].append((math.log(u) / L, math.log(a) / L, math.log(al) / L))
    fig = plt.figure(figsize=(5.6, 4.6))
    ax = fig.add_subplot(projection="3d")
    names = {"A": "A", "B": "B", "C": "C", "S": "strips", "O": "inner centre", "R": "other"}
    for k in ("R", "O", "S", "A", "B", "C"):
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


def ladder_figure(J):
    I = 1 / 84 - 3 / 3528
    labels = ["one\ncorner", "two\ncorners", "three\ncorners", "strips\nR1", "strips\nR1-R2",
              "strips\nR1-R3", "full\nstrips", "whole\ntriangle"]
    shares = [1 / 4, 1 / 2, 3 / 4, 0.75 + 6 * (I - J["var"]), 0.75 + 6 * (I - J["prime"]),
              0.75 + 6 * (I - J["chain"]), 40 / 49, 1]
    cs = [1.5 * (s * KAPPA / 4) ** (1 / 3) for s in shares]
    cols = [COL["A"], COL["B"], COL["C"], "#f2d27a", "#ecbe3c", COL["S"], COL["S"], COL["O"]]
    fig = plt.figure(figsize=(7.4, 4.0))
    ax = fig.add_subplot(projection="3d")
    for i, (c, col) in enumerate(zip(cs, cols)):
        ax.bar3d(1.5 * i, 0, 0, 0.6, 0.15, c, color=col, alpha=0.45 if i >= 6 else 0.95, shade=True, edgecolor="k")
        ax.text(1.5 * i - 0.25, 0.07, c + 0.03 + 0.08 * (i % 2), f"{c:.4f}", fontsize=7, zorder=30,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.6"))
    ax.set_xticks([1.5 * i + 0.3 for i in range(len(cs))])
    ax.set_xticklabels(labels, fontsize=7)
    ax.tick_params(axis="z", labelsize=8)
    ax.set_yticks([])
    ax.set_ylim(0, 1.2)
    ax.set_zlabel("$c$ in $N\\exp(-c(\\log N)^{2/3})$")
    ax.set_zlim(0, 0.9)
    ax.view_init(elev=12, azim=-82)
    ax.set_box_aspect((2.2, 0.5, 1))
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
    fig = plt.figure(figsize=(7.4, 4.6))
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
    ax = fig.add_subplot(1, 2, 1, projection="3d")
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
    ax = fig.add_subplot(1, 2, 2, projection="3d")
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
    fig.tight_layout()
    save(fig, "fig_engine.png")
    return grid, {k: len(v) for k, v in cells.items()}


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    simplex_figure()
    print("triples by corner:", triples_figure())
    J = strip_constants()
    print("constants:", ladder_figure(J))
    print("largest U on the grid:", strips_figure())
    print("engine grid:", engine_figure())
