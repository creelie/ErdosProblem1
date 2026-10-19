"""make_figures.py -- draws the two figures of the paper as PNG files (about 1600 px wide).

fig_ratio.png:  W(J)/(log J)^2 for 10 <= J <= 20000, the three-term formula, the limit
                9/(2 pi^2) and the constant 3/(2 pi^2) of the earlier lower bound.
fig_terms.png:  the summand 1/phi(4ua) over coprime pairs u, a <= 30, as a 3D bar chart.
fig_region.png: the surface log Z(sigma, tau) of the Dirichlet series of the inverted weight,
                with the point (0, 0) and the two slopes T, as a 3D schematic.
Requires numpy and matplotlib.
"""
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))
K = 9 / (2 * math.pi ** 2)
K1 = 3 / (2 * math.pi ** 2)
G = 0.57721566490153286
T = 0.3131549610
M = -1.0429405079


def phis(n):
    ph = np.arange(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if ph[p] == p:
            ph[p::p] -= ph[p::p] // p
    return ph


def ratio_figure():
    J = 20000
    ph = phis(J).astype(float)
    xs, ys = [], []
    marks = set(np.unique(np.round(np.logspace(1, math.log10(J), 120)).astype(int)))
    W = 0.0
    for m in range(1, J + 1):
        x = np.arange(1, m + 1)
        cop = np.gcd(m, x) == 1
        x = x[cop]
        c = np.where((m * x) % 2 == 0, 4.0, 2.0)
        t = 1.0 / (c * ph[m] * ph[x])
        W += 2 * t.sum() - (t[-1] if x[-1] == m else 0.0)
        if m in marks:
            xs.append(m)
            ys.append(W / math.log(m) ** 2)
    xs = np.array(xs, float)
    L = np.log(xs)
    f = K * ((L + G + T) ** 2 + M) / L ** 2
    fig, ax = plt.subplots(figsize=(8, 4.6), dpi=200)
    ax.plot(xs, ys, "o", ms=3.2, color="#1f4e79", label=r"$W(J)/(\log J)^2$, computed")
    ax.plot(xs, f, "-", lw=1.6, color="#c55a11", label=r"$\kappa\,((\log J+\gamma+T)^2+M)/(\log J)^2$")
    ax.axhline(K, ls="--", lw=1.2, color="#2e7d32", label=r"limit $\kappa=9/(2\pi^2)$")
    ax.axhline(K1, ls=":", lw=1.4, color="#7f7f7f", label=r"earlier bound $3/(2\pi^2)$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$J$")
    ax.set_ylim(0.1, 1.15)
    ax.set_xlim(10, J * 1.3)
    ax.grid(alpha=0.25, which="both")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_ratio.png"), facecolor="white")


def terms_figure():
    n = 30
    ph = phis(4 * n * n)
    fig = plt.figure(figsize=(8, 5.6), dpi=200)
    ax = fig.add_subplot(111, projection="3d")
    us, as_, hs, cs = [], [], [], []
    for u in range(1, n + 1):
        for a in range(1, n + 1):
            if math.gcd(u, a) != 1:
                continue
            us.append(u)
            as_.append(a)
            hs.append(1.0 / ph[4 * u * a])
            cs.append("#c55a11" if (u * a) % 2 else "#1f4e79")
    ax.bar3d(np.array(us) - 0.4, np.array(as_) - 0.4, 0, 0.8, 0.8, hs, color=cs, shade=True, linewidth=0)
    ax.set_xlabel(r"$u$")
    ax.set_ylabel(r"$a$")
    ax.set_zlabel(r"$1/\varphi(4ua)$")
    ax.view_init(elev=28, azim=-58)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_terms.png"), facecolor="white")


def primes(n):
    sieve = np.ones(n + 1, bool)
    sieve[:2] = False
    for q in range(2, int(n ** 0.5) + 1):
        if sieve[q]:
            sieve[q * q::q] = False
    return np.nonzero(sieve)[0].astype(float)


def logZ(s, t, ps):
    v = math.log(1 - 2.0 ** (-2 - s - t))
    q = ps[1:]
    e = 1 + 1 / ((q - 1) * q ** (1 + s)) + 1 / ((q - 1) * q ** (1 + t)) - (q + 1) / ((q - 1) * q ** (2 + s + t))
    return v + np.log(e).sum()


def region_figure():
    ps = primes(200000)
    g = np.linspace(-0.4, 0.8, 41)
    S, Tt = np.meshgrid(g, g)
    Zs = np.vectorize(lambda a, b: logZ(a, b, ps))(S, Tt)
    fig = plt.figure(figsize=(8, 5.8), dpi=200)
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    ax.plot_surface(S, Tt, Zs, cmap="viridis", alpha=0.75, zorder=1, linewidth=0, rstride=1, cstride=1)
    z0 = math.log(9 / math.pi ** 2)
    ax.scatter([0], [0], [z0], color="#c00000", s=40, depthshade=False, zorder=5)
    d = np.linspace(-0.35, 0.35, 2)
    ax.plot(d, 0 * d, z0 + T * d, color="#c00000", lw=2, zorder=5)
    ax.plot(0 * d, d, z0 + T * d, color="#c00000", lw=2, zorder=5)
    ax.text(0.05, 0.05, z0 + 0.12, r"$\log(9/\pi^2)$, slopes $T$", color="#c00000", zorder=6)
    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel(r"$\tau$")
    ax.set_zlabel(r"$\log Z(\sigma,\tau)$")
    ax.view_init(elev=24, azim=-128)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_region.png"), facecolor="white")


if __name__ == "__main__":
    region_figure()
    ratio_figure()
    terms_figure()
    for name in ("fig_ratio.png", "fig_terms.png", "fig_region.png"):
        print("wrote", os.path.join(OUT, name))
