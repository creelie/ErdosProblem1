"""centre.py -- the constant of the centre strips in "An explicit exponent for the exceptional set of
the Erdos-Straus equation", by deterministic quadrature, with a Monte Carlo check.

For a depth delta in (0, 1/42) put al = 2 delta, be = tau(delta) = min(1/10 - 7 delta/5, 1/2 - 19 delta)
(Maynard's window) and S = 1/2 + delta.  A pair (u, a) with u ~ x^s1, a ~ x^s2, s1 + s2 = S, is
counted when the modulus 4ua has a divisor in (x^al, x^be).  The rules below are sufficient
conditions for this; the loss J is the measure of the pairs that a rule does not count, and the
constant is c = (3/2) ((kappa/4) (3/4 + 6 (1/84 - 3/3528 - J)))^(1/3), kappa = 9/(2 pi^2).

    var    u or a lies in the window
    prime  ua has a prime factor in the window, or u or a lies in it
    smooth ... or the x^al-smooth part of ua lies in the window
    chain  ... or q T > x^al, where ga = min(be - al, al), T is the x^ga-smooth part of ua and q the
           largest prime factor of ua in (x^ga, x^al] (q = 1 if there is none)   [the paper's rule]

Densities: the x^g-smooth part of a number has density rho(t/g) (Dickman), a number with all prime
factors at least x^b has density omega(r/b)/b (Buchstab), and a number all of whose prime factors
lie in A = (ga, M] u [be, oo) has the density phi(r; M) with r phi(r) = 1_A(r) + int_{z in A, z < r}
phi(r - z) dz.  For the chain rule the uncovered density at (s1, s2) is
    g(s1, s2) = int int_{t1 + t2 <= al} rho(t1/ga) rho(t2/ga) phi(s1 - t1; M) phi(s2 - t2; M) dt1 dt2,
    M = max(al - t1 - t2, ga).
Usage: python3 centre.py [NG NH MC]   (Gauss nodes per piece, grid steps per ga, Monte Carlo samples)
"""
import math
import sys

import numpy as np

NG = 6
NH = 16
KAPPA = 9 / (2 * math.pi ** 2)
ICEIL = 1 / 84 - 3 / 3528
KMAX = 16

# ---- Dickman rho, Buchstab omega and Psi(u) = int_1^u omega on a fine grid -----------------------
H0 = 2e-5
U = 40.0
N0 = int(round(U / H0)) + 1
X0 = np.arange(N0) * H0
K1 = int(round(1 / H0))


def cumtrap(y, h):
    c = np.zeros_like(y)
    c[1:] = np.cumsum((y[1:] + y[:-1]) * h / 2)
    return c


RHO = np.ones(N0)
for _k in range(1, int(U)):
    _a, _b = _k * K1, min((_k + 1) * K1, N0 - 1)
    RHO[_a:_b + 1] = RHO[_a] - cumtrap(RHO[_a - K1:_b - K1 + 1] / X0[_a:_b + 1], H0)
OM = np.zeros(N0)
OM[K1:2 * K1 + 1] = 1 / X0[K1:2 * K1 + 1]
for _k in range(2, int(U)):
    _a, _b = _k * K1, min((_k + 1) * K1, N0 - 1)
    OM[_a:_b + 1] = (X0[_a] * OM[_a] + cumtrap(OM[_a - K1:_b - K1 + 1], H0)) / X0[_a:_b + 1]
PSI = np.zeros(N0)
PSI[K1:] = cumtrap(OM[K1:], H0)


def rho(t):
    t = np.asarray(t, float)
    return np.where(t <= 1, 1.0, np.where(t <= 2, 1 - np.log(np.maximum(t, 1)), np.interp(t, X0, RHO, right=0.0)))


def omega(t):
    t = np.asarray(t, float)
    return np.where(t < 1, 0.0, np.where(t <= 2, 1 / np.maximum(t, 1), np.interp(t, X0, OM)))


def psi(t):
    t = np.asarray(t, float)
    return np.where(t < 1, 0.0, np.where(t <= 2, np.log(np.maximum(t, 1)), np.interp(t, X0, PSI)))


GX, GW = np.polynomial.legendre.leggauss(NG)


def configure(ng, nh):
    """Set the number of Gauss nodes per piece and of grid steps per ga."""
    global NG, NH, GX, GW
    NG, NH = ng, nh
    GX, GW = np.polynomial.legendre.leggauss(NG)


def nodes(a, b, cuts):
    """Gauss-Legendre nodes on [a, b], split at the cuts (the kinks of the integrand)."""
    pts = sorted(set([a, b] + [c for c in cuts if a < c < b]))
    xs, ws = [], []
    for lo, hi in zip(pts[:-1], pts[1:]):
        if hi - lo < 1e-14:
            continue
        xs.append((hi - lo) / 2 * GX + (hi + lo) / 2)
        ws.append((hi - lo) / 2 * GW)
    if not xs:
        return np.zeros(0), np.zeros(0)
    return np.concatenate(xs), np.concatenate(ws)


def tau(d):
    return min(0.1 - 1.4 * d, 0.5 - 19 * d)


def s_nodes(d, extra=()):
    al, be, S = 2 * d, tau(d), 0.5 + d
    cuts = []
    for k in range(1, KMAX):
        cuts += [k * be, k * be + al, S - k * be, S - k * be - al]
        cuts += [k * be + e for e in extra] + [S - k * be - e for e in extra]
    return nodes(be, S - be, cuts)


# ---- the rules -----------------------------------------------------------------------------------
def f_noprime(s, al, be):
    """Density of n ~ x^s with no prime factor in (x^al, x^be): x^al-smooth part times be-rough part."""
    if s < be:
        return float(rho(s / al))
    cuts = [s - k * be for k in range(1, KMAX)] + [k * al for k in range(1, 30)]
    v, w = nodes(0.0, s - be, cuts)
    return float(rho(s / al) + (w * rho(v / al) * omega((s - v) / be)).sum() / be)


def d_smooth_in(s1, s2, al, be):
    """Mass with no prime factor in the window and x^al-smooth part of ua of size in (al, be)."""
    cuts1 = [s1 - k * be for k in range(1, KMAX)] + [al, be - al]
    v1, w1 = nodes(0.0, be, cuts1)
    tot = 0.0
    for a1, b1 in zip(v1, w1):
        lo, hi = max(0.0, al - a1), be - a1
        if hi <= lo:
            continue
        v2, w2 = nodes(lo, hi, [s2 - k * be for k in range(1, KMAX)] + [al])
        tot += b1 * rho(a1 / al) * omega((s1 - a1) / be) * (w2 * rho(v2 / al) * omega((s2 - v2) / be)).sum()
    return float(tot) / be ** 2


class Tables:
    """C[j, k] = int_0^{r_j} phi(r; M_k): phi is the density of numbers > 1 all of whose prime factors lie
    in A_k = (ga, M_k] u [be, oo), from r phi(r) = 1_A(r) + int_{z in A, z < r} phi(r - z) dz."""

    def __init__(self, ga, be, Ms, rmax):
        self.ga, self.be, self.Ms = ga, be, np.asarray(Ms)
        h = ga / NH
        n = int(math.ceil(rmax / h)) + 2
        K = len(Ms)
        C = np.zeros((n + 1, K))
        r = np.arange(n + 1) * h
        self.r = r
        bo = be / h
        bi, bf = int(math.floor(bo)), bo - math.floor(bo)
        mo = self.Ms / h
        mi = np.floor(mo).astype(int)
        mf = mo - mi
        kk = np.arange(K)

        def at(j, off_i, off_f):  # C at r_j - off * h by linear interpolation, 0 below 0
            jj = j - off_i
            a = np.where(jj >= 0, C[np.clip(jj, 0, n), kk], 0.0)
            b = np.where(jj - 1 >= 0, C[np.clip(jj - 1, 0, n), kk], 0.0)
            return (1 - off_f) * a + off_f * b

        def q(j):  # int over z in A, z < r_j, of phi(r_j - z)
            return at(j, NH, 0.0) - at(j, mi, mf) + at(j, bi, bf)

        qprev = np.zeros(K)
        for j in range(NH, n):
            r0, r1 = r[j], r[j + 1]
            a1, b1 = np.clip(r0, ga, self.Ms), np.clip(r1, ga, self.Ms)
            ind = np.log(b1 / a1) + (math.log(r1 / max(r0, be)) if r1 > be else 0.0)
            qn = q(j + 1)
            C[j + 1] = C[j] + ind + h / 2 * (qprev / r0 + qn / r1)
            qprev = qn
        self.C = C

    def phi(self, rr, k):
        rr = np.asarray(rr, float)
        C = self.C[:, k]
        M = self.Ms[k]
        cx = lambda xv: np.interp(xv, self.r, C, left=0.0)
        q = cx(rr - self.ga) - cx(rr - M) + cx(rr - self.be)
        ind = ((rr > self.ga) & (rr <= M)) | (rr >= self.be)
        return np.where(rr > self.ga, (ind + q) / np.maximum(rr, 1e-300), 0.0)


def chain_setup(d):
    al, be = 2 * d, tau(d)
    ga = min(be - al, al)
    tk, tw = nodes(0.0, min(al, 24 * ga), [k * ga for k in range(1, 25)] + [al - ga])
    T = Tables(ga, be, np.maximum(al - tk, ga), 0.5 + d)
    t1n = [nodes(0.0, t, [k * ga for k in range(1, 25)] + [t - k * ga for k in range(1, 25)]) for t in tk]
    return ga, tk, tw, T, t1n


def g_chain(s1, s2, setup):
    ga, tk, tw, T, t1n = setup
    g = 0.0
    for k, (t, wt) in enumerate(zip(tk, tw)):
        t1, w1 = t1n[k]
        g += wt * (w1 * rho(t1 / ga) * rho((t - t1) / ga) * T.phi(s1 - t1, k) * T.phi(s2 - t + t1, k)).sum()
    return g


def loss(d, rule):
    """The uncovered measure at depth d, integrated over s1 in (be, S - be)."""
    al, be, S = 2 * d, tau(d), 0.5 + d
    if rule == "var":
        return max(0.0, S - 2 * be)
    if rule == "chain":
        setup = chain_setup(d)
        ga = setup[0]
        s, w = s_nodes(d, extra=(ga,))
        return sum(wi * g_chain(si, S - si, setup) for si, wi in zip(s, w))
    s, w = s_nodes(d)
    tot = 0.0
    for si, wi in zip(s, w):
        g = f_noprime(si, al, be) * f_noprime(S - si, al, be)
        if rule == "smooth":
            g -= d_smooth_in(si, S - si, al, be)
        tot += wi * g
    return tot


def J(rule):
    tot = 0.0
    for a, b in ((0, 1 / 54), (1 / 54, 1 / 44), (1 / 44, 1 / 42)):
        dd, dw = nodes(a, b, [])
        tot += sum(wi * loss(di, rule) for di, wi in zip(dd, dw))
    return tot


def constant(Jval):
    return 1.5 * ((KAPPA / 4) * (0.75 + 6 * (ICEIL - Jval))) ** (1 / 3)


def j_var_exact():
    """J for the rule var in closed form: int (S - 2 be)_+ d delta, with S - 2 be = 3 delta + 3 delta/... ."""
    # S - 2 be = 1/2 + delta - 2/10 + 14 delta/5 = 3/10 + 19 delta/5 on (0, 1/44],
    # S - 2 be = 1/2 + delta - 1 + 38 delta = -1/2 + 39 delta on (1/44, 1/42), positive for delta > 1/78.
    part1 = 3 / 10 * (1 / 44) + 19 / 10 * (1 / 44) ** 2
    part2 = -1 / 2 * (1 / 42 - 1 / 44) + 39 / 2 * ((1 / 42) ** 2 - (1 / 44) ** 2)
    return part1 + part2


# ---- Monte Carlo check of the chain density -------------------------------------------------------
RNG = np.random.default_rng(5)


def pieces(s, eps):
    """Prime factors of a random number of size x^s: stick breaking (Poisson-Dirichlet), dust < eps."""
    out, rem = [], s
    while rem > eps:
        v = rem * RNG.random()
        out.append(v)
        rem -= v
    return out, rem


def mc_point(d, s1, n):
    al, be = 2 * d, tau(d)
    s2 = 0.5 + d - s1
    ga = min(be - al, al)
    hits = 0
    for _ in range(n):
        ps, dust = [], 0.0
        for s in (s1, s2):
            p, r = pieces(s, 1e-6)
            ps += p
            dust += r
        if any(al < p < be for p in ps):
            continue
        t = dust + sum(p for p in ps if p <= ga)
        q = max([p for p in ps if ga < p <= al], default=0.0)
        if q + t <= al:
            hits += 1
    return hits / n


def main(nmc):
    print(f"quadrature: {NG} Gauss nodes per piece, {NH} grid steps per ga")
    jv = J("var")
    assert abs(jv - j_var_exact()) < 1e-12, (jv, j_var_exact())
    results = {"var": jv}
    for rule in ("prime", "smooth", "chain"):
        results[rule] = J(rule)
    print(f"ceiling (every pair of the strips): share 40/49, c = {constant(0.0):.10f}")
    for rule, jval in results.items():
        print(f"rule {rule:6s}: J = {jval:.10f}, share = {0.75 + 6 * (ICEIL - jval):.7f}, c = {constant(jval):.9f}")
    assert abs(constant(0.0) - 1.5 * (45 / (49 * math.pi ** 2)) ** (1 / 3)) < 1e-12
    assert results["var"] > results["prime"] > results["smooth"] > results["chain"] > 0
    c1 = constant(results["chain"])
    assert 0.67920 < c1 < 0.67922, c1
    # Monte Carlo against the quadrature at four points
    for d, s1 in ((0.02, 0.2), (0.022, 0.3), (0.0235, 0.15), (0.019, 0.26)):
        setup = chain_setup(d)
        gq = g_chain(s1, 0.5 + d - s1, setup)
        gm = mc_point(d, s1, nmc)
        sd = math.sqrt(max(gm * (1 - gm), 1e-12) / nmc)
        print(f"  g at delta = {d}, s1 = {s1}: quadrature {gq:.5f}, Monte Carlo {gm:.5f} +- {sd:.5f}")
        assert abs(gq - gm) < 5 * sd + 2e-3
    print(f"OK: c1 = {c1:.6f} (chain rule), above c0 = {1.5 * (27 / (32 * math.pi ** 2)) ** (1 / 3):.6f}")


if __name__ == "__main__":
    configure(int(sys.argv[1]) if len(sys.argv) > 1 else 6, int(sys.argv[2]) if len(sys.argv) > 2 else 16)
    main(int(sys.argv[3]) if len(sys.argv) > 3 else 100000)
