"""The deep strips of the centre for the Erdos-Straus large sieve, assembled with the typed crossing,
the memory and the learn-from-failure log of the engine Anaxia (creelie/eBLLM, version 3.1.0).

The run behind the strips (anaxia_centre.py) counts the centre up to the depth 1/42, where the window
of Maynard's first theorem closes.  This run asks what can count the rest of the centre, the depths
1/42 <= delta < 1/6.  A pair (u, a) of the centre gives the modulus q = 4ua of size x^(1/2 + delta).
Maynard's second theorem counts primes p = -1 (mod q) for moduli whose weight is a convolution of
three sequences of sizes Q1, Q2, Q3 that meet four inequalities for the size N of the long variable
of Heath-Brown's identity.  No triple of sizes meets them beyond the depth 1/10, so the grid of depths
stops there.  The library holds:

    the large sieve, the corners of the first run and the strips of the second;
    the deep strips left out, counted by Maynard's window continued past 1/42, split for Maynard's
        second theorem, or counted under the Elliott-Halberstam conjecture;
    a splitting rule: the prime factors of q above x^G go to three groups, and the primes below x^G
        fill groups to prescribed sizes, which they overshoot by less than G:
        one      one group is filled;
        two      two groups are filled, and the budget keeps a reserve G for the overshoot of the first;
        noreserve  two groups are filled, without the reserve;
    the set of values of G (0.04, 0.03, 0.02, or the best of 0.015, 0.02, 0.025, 0.03 and 0.04 at
        each depth);
    the inputs (Maynard's first theorem, the second theorem in the form proved in the paper, the
        Elliott-Halberstam conjecture).

Every closed composite is checked.  A splitting rule is tested on integers: for every n in
[N0, N0 + NINT], every depth on a grid and random assignments that the rule allows, the primes of n
below x^G are added in increasing order to the filled groups, and every filled group must end with
a size in [target, target + G); a rule that allows an assignment the fill cannot reach is refuted.
The window rule is tested by the length of the window at each depth.  The constant of an accepted
composite comes from a Monte Carlo estimate of the share of factorable pairs (checks/c/deep_mc.c);
the proved bound is that of the certificate (checks/c/certify_K.c).  Verdicts go to the memory of
the engine, the failures to its learn-from-failure log, and both are exported to anaxia_deep.json.

Run:  pip install "git+https://github.com/creelie/eBLLM@v3.1.0"; python3 anaxia_deep.py [NINT NMC]
"""
import json
import math
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import ebllm
from ebllm.bllm.memory import Memory
from ebllm.bllm.selfmod import Change, Record
from ebllm.bllm.targets import Con, cross

HERE = os.path.dirname(os.path.abspath(__file__))
NINT = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
NMC = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
N0 = 10 ** 9
PROB = "erdos-straus"
KAPPA = 9 / (2 * math.pi ** 2)
J_QUAD = 0.000306          # the loss of the strips, by quadrature (Section 9 of the paper)
DEPTHS = [1 / 42] + [round(0.025 + 0.0025 * k, 4) for k in range(31)]   # 1/42 .. 0.1; nothing splits beyond 1/10
RULES = ("one", "two", "noreserve")
GSETS = {"g04": (0.04,), "g03": (0.03,), "g02": (0.02,), "gbest": (0.015, 0.02, 0.025, 0.03, 0.04)}


def node(i, gives, needs=(), source=""):
    return Con(i, PROB, gives, frozenset(needs), source, (i,))


LIB = [node("es.sieve", "bound", ["region:corners", "region:centre", "region:deep"],
            "Montgomery's large sieve with Q = sqrt(N) and the tilted lower bound for L(Q)"),
       node("es.corners", "region:corners", [], "the three corners with the Farey, M and I forms (first run)"),
       node("es.strips", "region:centre", [], "the strips of depth < 1/42 with the chain rule (second run)"),
       node("es.deep.none", "region:deep", [], "the deep strips left out"),
       node("es.deep.window", "region:deep", ["level:maynard"],
            "a divisor of q in Maynard's window (x^(2 delta), x^theta(delta)), continued past 1/42"),
       node("es.deep.split", "region:deep", ["split:rule", "split:G", "level:maynard2"],
            "the moduli split into three factors for Maynard's second theorem"),
       node("es.deep.eh", "region:deep", ["level:eh"], "every pair of the deep strips")]
for r in RULES:
    LIB.append(node(f"es.split.{r}", "split:rule", [], f"the {r} splitting rule"))
for g, vals in GSETS.items():
    LIB.append(node(f"es.G.{g}", "split:G", [], "G in " + ", ".join(str(v) for v in vals)))
LIB += [node("es.maynard", "level:maynard", [], "Maynard's first theorem: a divisor in the window, delta < 1/42"),
        node("es.maynard2", "level:maynard2", [],
             "Maynard's second theorem for weights that are convolutions of three factors of sizes "
             "admissible for every size of the long variable (Theorem 7.1 of the paper)"),
        node("es.eh", "level:eh", [], "Elliott-Halberstam conjecture (moduli up to p^(2/3))")]


def breed(lib):
    """Typed crossing to a fixed point (as in anaxia_regions.py)."""
    key = lambda c: (c.gives, tuple(sorted(c.parts)))
    pool = {key(c): c for c in lib}
    front = list(lib)
    while front:
        new = []
        for a in front:
            for b in lib:
                c = cross(a, b)
                if c is not None and key(c) not in pool:
                    pool[key(c)] = c
                    new.append(c)
        front = new
    return [c for c in pool.values() if c.closed and c.gives == "bound"]


# -- the integers n in [N0, N0 + NINT] with their prime factors ---------------------------------------
FACT = None


def factored():
    global FACT
    if FACT is None:
        rest = list(range(N0, N0 + NINT + 1))
        fac = [[] for _ in rest]
        lim = int(math.isqrt(N0 + NINT)) + 1
        sieve = bytearray([1]) * (lim + 1)
        for p in range(2, lim + 1):
            if not sieve[p]:
                continue
            sieve[p * p::p] = bytearray(len(sieve[p * p::p]))
            for i in range((-N0) % p, NINT + 1, p):
                while rest[i] % p == 0:
                    rest[i] //= p
                    fac[i].append(p)
        for i, r in enumerate(rest):
            if r > 1:
                fac[i].append(r)
        FACT = fac
    return FACT


def fill_failures(rule, gvals, seed=1):
    """Assignments allowed by the rule whose fill by the small primes misses its box."""
    rng = random.Random(seed)
    tried = fails = 0
    first = None
    for idx, ps in enumerate(factored()):
        n = N0 + idx
        for _ in range(2):
            delta = rng.choice(DEPTHS)
            G = rng.choice(gvals)
            lx = math.log(n) / (0.5 + delta)
            sizes = sorted(math.log(p) / lx for p in ps)
            small = [z for z in sizes if z < G]
            big = [z for z in sizes if z >= G]
            t = sum(small)
            groups = [0.0, 0.0, 0.0]
            for z in big:
                groups[rng.randrange(3)] += z
            L = rng.randrange(3)
            others = [i for i in range(3) if i != L]
            nf = rng.choice((1, 2)) if rule != "one" else 1
            F = others[:nf] if rng.random() < 0.5 else others[::-1][:nf]
            reserve = G if (rule == "two" and nf == 2) else 0.0
            room = t - reserve
            if room <= 0:
                continue
            w = [rng.random() for _ in F]
            scale = rng.random() * room / max(sum(w), 1e-12)
            targets = {i: groups[i] + wi * scale for i, wi in zip(F, w)}
            if sum(targets[i] - groups[i] for i in F) + reserve >= t:
                continue
            tried += 1
            items = list(small)
            got = dict(groups=list(groups))
            k = 0
            ok = True
            for i in F:
                size = groups[i]
                while size < targets[i] and k < len(items):
                    size += items[k]
                    k += 1
                if not (targets[i] <= size < targets[i] + G):
                    ok = False
                got[i] = size
            if not ok:
                fails += 1
                first = first or dict(n=n, depth=delta, G=G, F=F, targets=targets, small=small)
    return tried, fails, first


# -- Monte Carlo share of the factorable pairs --------------------------------------------------------
MC = {}
BIN = os.path.join(HERE, "deep_mc")


def mc_fraction(delta, G, mode):
    key = (delta, G, mode)
    if key not in MC:
        out = subprocess.run([BIN, str(delta), str(NMC), str(G), "0", "12", "11", "10", str(mode)],
                             capture_output=True, text=True, check=True).stdout
        MC[key] = float(out.split("rule=")[1].split()[0])
    return MC[key]


def mc_prefetch(keys):
    todo = [k for k in keys if k not in MC]
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 1) as pool:
        list(pool.map(lambda k: mc_fraction(*k), todo))


def k_estimate(rule, gvals):
    mode = {"two": 0, "one": 1, "noreserve": 2}[rule]
    mc_prefetch([(d, G, mode) for d in DEPTHS for G in gvals])
    vals = [max(mc_fraction(d, G, mode) for G in gvals) * (0.5 - 3 * d) for d in DEPTHS]
    return sum((DEPTHS[i + 1] - DEPTHS[i]) * (vals[i] + vals[i + 1]) / 2 for i in range(len(DEPTHS) - 1))


def constant(K):
    return 1.5 * (KAPPA / 4 * (40 / 49 - 6 * J_QUAD + 6 * K)) ** (1 / 3)


def window_length(delta):
    return min(0.1 - 1.4 * delta, 0.5 - 19 * delta) - 2 * delta


def check(c):
    parts = set(c.parts)
    if "es.deep.none" in parts:
        return "accepted", f"c = {constant(0):.6f}", dict(deep="none", K=0, c=round(constant(0), 6))
    if "es.deep.eh" in parts:
        cc = 1.5 * (KAPPA / 4) ** (1 / 3)
        return "conditional", f"c = {cc:.6f} if the Elliott-Halberstam conjecture holds", dict(deep="eh", c=round(cc, 6))
    if "es.deep.window" in parts:
        longest = max(window_length(d) for d in DEPTHS)
        meas = dict(deep="window", longest_window=round(longest, 6))
        return "rejected", f"the window is empty at every depth from 1/42 on (longest {longest:.4f})", meas
    rule = next(r for r in RULES if f"es.split.{r}" in parts)
    gset = next(g for g in GSETS if f"es.G.{g}" in parts)
    tried, fails, first = fill_failures(rule, GSETS[gset])
    meas = dict(deep="split", rule=rule, G=list(GSETS[gset]), assignments_tested=tried, fill_failures=fails,
                first_failure=first, nint=NINT, n0=N0)
    if fails:
        return "rejected", f"the {rule} rule allows {fails} of {tried} assignments that the fill cannot reach", meas
    K = k_estimate(rule, GSETS[gset])
    meas.update(K_estimate=round(K, 6), c=round(constant(K), 6), nmc=NMC)
    return "accepted", f"K = {K:.5f}, c = {constant(K):.6f} (Monte Carlo)", meas


def main():
    src = os.path.join(HERE, "..", "checks", "c", "deep_mc.c")
    subprocess.run(["gcc", "-O2", "-o", BIN, src, "-lm"], check=True)
    mem = Memory()
    mem.add_nodes(LIB)
    closed = breed(LIB)
    out = []
    for c in sorted(closed, key=lambda c: c.id):
        outcome, reason, meas = check(c)
        mem.record(c.parts, PROB, outcome, node="ES.deep", statement=reason, reason=json.dumps(meas),
                   survived=(outcome == "accepted"))
        out.append(dict(composite=c.id, parts=list(c.parts), outcome=outcome, reason=reason, measurement=meas))
    acc = sorted((r for r in out if r["outcome"] == "accepted"), key=lambda r: -r["measurement"]["c"])
    best = acc[0]["measurement"]
    window = next(r for r in out if r["measurement"].get("deep") == "window")
    norez = [r for r in out if r["measurement"].get("rule") == "noreserve"]
    records = [
        Record("Explicit exponent for E(N), the deep strips", 1, 0, "failed", "Maynard's window past 1/42",
               window["reason"], "missing-mode",
               Change("input for the deep strips", "library", "primes to moduli with a divisor in a window",
                      "Maynard's second theorem for moduli that split into three factors of admissible sizes",
                      "share of the centre counted"), "partial"),
        Record("Explicit exponent for E(N), the deep strips", 2, 1, "partial", "the noreserve splitting rule",
               norez[0]["reason"] if norez else "", "overclaim",
               Change("fill of two groups", "check", "fill two groups up to their targets",
                      "keep a reserve G for the overshoot of the first group", "assignments the fill reaches"),
               "partial", note=f"best accepted: rule {best['rule']}, G {best['G']}, K = {best['K_estimate']}"),
    ]
    for rec in records:
        mem.record_attempt(rec.as_dict())
    mem.commit()
    summary = dict(engine=f"Anaxia (ebllm) {ebllm.__version__}", library=len(LIB), closed=len(out),
                   accepted=len(acc), rejected=sum(r["outcome"] == "rejected" for r in out),
                   conditional=sum(r["outcome"] == "conditional" for r in out), best=best,
                   accepted_constants=sorted({r["measurement"]["c"] for r in acc}), memory=mem.stats(),
                   modlog=[r.as_dict() for r in records], depths=DEPTHS, nint=NINT, nmc=NMC,
                   mc={f"{d}:{G}:{m}": v for (d, G, m), v in sorted(MC.items())})
    json.dump(dict(summary=summary, composites=out), open(os.path.join(HERE, "anaxia_deep.json"), "w"), indent=1)
    print(json.dumps({k: summary[k] for k in ("engine", "library", "closed", "accepted", "rejected",
                                                "conditional", "best", "accepted_constants")}, indent=1))
    assert best["rule"] in ("two", "one") and best["fill_failures"] == 0
    print("OK: best accepted composite splits the moduli with rule", best["rule"], "and G in", best["G"],
          ", K about", best["K_estimate"], ", c about", best["c"])


if __name__ == "__main__":
    main()
