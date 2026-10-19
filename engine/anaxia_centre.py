"""The centre of the triangle for the Erdos-Straus large sieve, assembled with the typed crossing and
the memory of the engine Anaxia (creelie/eBLLM, version 3.0.5).

The run behind the three corners (anaxia_regions.py) left the centre out: there the moduli 4ua,
4u alpha, 4a alpha exceed sqrt(p), and two triples of the centre share a class at p = 531359.  This
run adds the centre as strips next to the three corners.  A centre triple (u, a, alpha) has
u, a, alpha < sqrt(P); it is assigned to the pair of its two smaller coordinates, whose product is
P^(1/2 + delta) with depth delta, and the strip of depth delta < 1/42 is counted with Maynard's
theorem on primes in progressions to moduli up to x^(1/2 + delta) that have a divisor in the window
(x^(2 delta), x^tau(delta)), tau(delta) = min(1/10 - 7 delta/5, 1/2 - 19 delta).

The library holds:
    the large sieve, the three corners of the earlier run, the centre left out or taken as strips;
    a normal form of the class in the centre (farey u/a, M = 4u^2 alpha, I = (4a^2 alpha)^-1);
    a rule that decides which pairs of a strip have a divisor in the window
        var     u or a lies in the window
        prime   ... or ua has a prime factor in the window
        smooth  ... or the x^(2 delta)-smooth part of ua lies in the window
        chain   ... or q T > x^(2 delta)  (T the smooth part below the step bound, q the largest
                prime factor of ua between the step bound and x^(2 delta))
        steps   ... or the x^(2 delta)-smooth part of ua exceeds x^(2 delta)  (no step bound)
        eh      every pair (Elliott-Halberstam conjecture)
    the input for the rules (Maynard's theorem, or the Elliott-Halberstam conjecture);
    the treatment of coincidences between the centre and the corners: none allowed, or those of
    the five-integer family (a'', alpha'', k, D, j) of the paper, whose number has a convergent
    weighted sum.

Every closed composite is checked on data: the normal form must separate the classes of the centre
for the primes p = 3 (mod 4) up to PMAX; every coincidence of a centre class with another class
must be allowed by the coincidence node (and, for the family, satisfy its identities); the rule must
never claim a divisor in a window when there is none, for all n up to NRULE and 64 windows.  The
constant is computed by the quadrature of checks/python/centre.py.  Verdicts go to the memory of
the engine, which is exported to anaxia_centre.json.

Run:  pip install "git+https://github.com/creelie/eBLLM@v3.0.5"; python3 anaxia_centre.py [PMAX NRULE]
"""
import json
import math
import os
import sys

from ebllm.bllm.targets import Con, cross
from ebllm.bllm.memory import Memory
import ebllm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "checks", "python"))
import centre as quad  # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
NRULE = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
KAPPA = 9 / (2 * math.pi ** 2)
C0 = 1.5 * (27 / (32 * math.pi ** 2)) ** (1 / 3)
FORMS = ("farey", "M", "I")
RULES = ("var", "prime", "smooth", "chain", "steps")
PROB = "erdos-straus"


def node(i, gives, needs=(), source=""):
    return Con(i, PROB, gives, frozenset(needs), source, (i,))


LIB = [node("es.sieve", "bound", ["region:corners", "region:centre"],
            "Montgomery's large sieve with Q = sqrt(N) and the tilted lower bound for L(Q)"),
       node("es.corners", "region:corners", [], "the three corners with the Farey, M and I forms (earlier run)"),
       node("es.centre.none", "region:centre", [], "the centre left out")]
for f in FORMS:
    LIB.append(node(f"es.centre.{f}", "region:centre", [f"cform:{f}", "count:strips", "sparse:centre"],
                    f"strips of depth < 1/42 next to the corners, classes compared in the {f} form"))
    LIB.append(node(f"es.cform.{f}", f"cform:{f}", [], f"the {f} normal form on the centre"))
for r in RULES:
    LIB.append(node(f"es.rule.{r}", "count:strips", ["level:maynard"], f"the pairs of a strip that the {r} rule counts"))
LIB.append(node("es.rule.eh", "count:strips", ["level:eh"], "every pair of the centre"))
LIB.append(node("es.maynard", "level:maynard", [],
                "Maynard: moduli up to x^(1/2+delta) with a divisor in (x^(2delta), x^tau(delta)), delta < 1/42"))
LIB.append(node("es.eh", "level:eh", [], "Elliott-Halberstam conjecture (moduli up to p^(2/3))"))
LIB.append(node("es.sparse.none", "sparse:centre", [], "no class of the centre may be a class of another triple"))
LIB.append(node("es.sparse.family", "sparse:centre", [],
                "coincidences of the centre (u<a) with C and (u>a) with B in the five-integer family"))


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


# -- data: the triples of the centre and their coincidences ----------------------------------------
def primes_3mod4(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [p for p in range(3, n + 1, 4) if s[p]]


def divisors(n):
    small, large = [], []
    d = 1
    while d * d <= n:
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
        d += 1
    return small + large[::-1]


def region(u, a, al, P, p):
    if u * u < P and a * a < P and al * al < P:
        return "centre"
    if u >= 2 and a >= 2 and 16 * (u * a) ** 2 <= p:
        return "A"
    if 16 * (u * al) ** 2 <= p:
        return "B"
    if 16 * (a * al) ** 2 <= p:
        return "C"
    return "rest"


def in_family(p, P, t, t2):
    """Is the coincidence of the centre triple t (u < a) with the triple t2 of C in the family?"""
    u, a, al = t
    u2, a2, al2 = t2
    s = a2 * al2
    num = u2 * a - u * a2
    if num <= 0 or num % p:
        return False
    k = num // p
    j = a - 4 * k * s
    D = s * a2 - a * al * j
    return (k >= 1 and j >= 1 and D >= 1 and P * j == s * (u * a2 - k) and u * D == s * k
            and al * j * (j + 4 * k * s) == s * a2 - D)


DATA = None


def data():
    """For every prime: the largest u and a in the centre, the M and I heights, and the coincidences."""
    global DATA
    if DATA is None:
        DATA = []
        for p in primes_3mod4(PMAX):
            P = (p + 1) // 4
            cls = {}
            hmax = {"farey": (0, 0), "M": 0, "I": 0}
            for u in divisors(P):
                for a in divisors(P // u):
                    if math.gcd(u, a) != 1:
                        continue
                    al = P // (u * a)
                    reg = region(u, a, al, P, p)
                    if reg == "rest":
                        continue
                    if reg == "centre":
                        mu, ma = hmax["farey"]
                        hmax["farey"] = (max(mu, u), max(ma, a))
                        hmax["M"] = max(hmax["M"], 4 * u * u * al)
                        hmax["I"] = max(hmax["I"], 4 * a * a * al)
                    cls.setdefault(u * pow(a, -1, p) % p, []).append((reg, (u, a, al)))
            coin = []
            for v in cls.values():
                cen = [t for r, t in v if r == "centre"]
                if not cen or len(v) < 2:
                    continue
                for r, t in v:
                    for t0 in cen:
                        if t0 == t or (r == "centre" and t0 > t):
                            continue
                        coin.append((t0, r, t))
            DATA.append((p, P, hmax, coin))
    return DATA


def form_ok(form):
    bad = 0
    for p, P, hmax, coin in data():
        if form == "farey":
            mu, ma = hmax["farey"]
            bad += mu * ma >= p
        else:
            bad += hmax[form] >= p
    return bad


def coincidences(sparse):
    """Returns (number of coincidences, number not allowed by the node, sum over them of 1/p)."""
    total = bad = 0
    weight = 0.0
    for p, P, hmax, coin in data():
        for t0, r, t in coin:
            total += 1
            weight += 1 / p
            ok = False
            if sparse == "family":
                if r == "C" and t0[0] < t0[1]:
                    ok = in_family(p, P, t0, t)
                elif r == "B" and t0[0] > t0[1]:
                    ok = in_family(p, P, (t0[1], t0[0], t0[2]), (t[1], t[0], t[2]))
            bad += not ok
    return total, bad, weight


# -- the rules on integers ---------------------------------------------------------------------------
WINDOWS = [(A, A * r) for A in (3.0, 7.5, 20.0, 55.0, 150.0, 400.0, 1100.0, 3000.0)
           for r in (1.1, 1.3, 1.7, 2.5, 4.0, 7.0, 15.0, 40.0)]
SPF = None


def factor(n):
    out = []
    while n > 1:
        p = SPF[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        out.append((p, e))
    return out


def rule_claims(rule, fac, A, B):
    """Does the rule claim a divisor of n in (A, B)?  (The pair is u a = n; the rules for u and a
    alone are among these, since u and a are divisors of n.)"""
    if any(A < p < B for p, _ in fac):
        return True
    if rule == "prime":
        return False
    smooth = math.prod(p ** e for p, e in fac if p <= A)
    if rule == "smooth":
        return A < smooth < B
    if rule == "steps":
        return smooth > A
    G = min(B / A, A)
    T = math.prod(p ** e for p, e in fac if p < G)
    q = max([p for p, _ in fac if G <= p <= A], default=1)
    return q * T > A  # chain


RULE_CACHE = {}


def rule_failures(rule):
    """Counterexamples: n <= NRULE and a window in which the rule claims a divisor that n lacks."""
    global SPF
    if rule in RULE_CACHE:
        return RULE_CACHE[rule]
    if SPF is None:
        SPF = list(range(NRULE + 1))
        for i in range(2, int(NRULE ** 0.5) + 1):
            if SPF[i] == i:
                for j in range(i * i, NRULE + 1, i):
                    if SPF[j] == j:
                        SPF[j] = i
    fails, first = 0, None
    for n in range(2, NRULE + 1):
        fac = factor(n)
        divs = None
        for A, B in WINDOWS:
            if rule == "var":
                continue  # u or a in the window: they divide ua, nothing to test
            if rule_claims(rule, fac, A, B):
                if divs is None:
                    divs = divisors(n)
                if not any(A < d < B for d in divs):
                    fails += 1
                    first = first or (n, A, B)
    RULE_CACHE[rule] = (fails, first)
    return RULE_CACHE[rule]


CONST = {}


def rule_constant(rule):
    if rule not in CONST:
        CONST[rule] = quad.constant(quad.J(rule)) if rule != "eh" else 1.5 * (KAPPA / 4) ** (1 / 3)
    return CONST[rule]


def check(c):
    parts = set(c.parts)
    if "es.centre.none" in parts:
        meas = dict(centre="none", c=C0, pmax=PMAX)
        return "accepted", f"c = {C0:.6f}", meas
    form = next(f for f in FORMS if f"es.centre.{f}" in parts)
    rule = next(r for r in RULES + ("eh",) if f"es.rule.{r}" in parts)
    sparse = "family" if "es.sparse.family" in parts else "none"
    fbad = form_ok(form)
    ncoin, nbad, wcoin = coincidences(sparse)
    rfail, rfirst = rule_failures(rule) if rule != "eh" else (0, None)
    const = rule_constant(rule) if not rfail else None
    meas = dict(centre="strips", form=form, rule=rule, sparse=sparse, form_failures=fbad,
                coincidences=ncoin, coincidences_not_allowed=nbad, coincidence_weight=round(wcoin, 9),
                rule_failures=rfail, rule_first_failure=rfirst, c=round(const, 9) if const else None,
                pmax=PMAX, nrule=NRULE)
    if fbad:
        return "rejected", f"the {form} form does not separate the classes of the centre on {fbad} primes", meas
    if nbad:
        return "rejected", f"{nbad} coincidences of the centre not allowed by the {sparse} node", meas
    if rfail:
        n, A, B = rfirst
        return "rejected", f"the {rule} rule claims a divisor in a window that has none ({rfail} cases, first n={n})", meas
    if rule == "eh":
        return "conditional", f"c = {const:.6f} if the Elliott-Halberstam conjecture holds", meas
    return "accepted", f"c = {const:.6f}", meas


def main():
    quad.configure(6, 16)
    mem = Memory()
    mem.add_nodes(LIB)
    closed = breed(LIB)
    out = []
    for c in sorted(closed, key=lambda c: c.id):
        outcome, reason, meas = check(c)
        mem.record(c.parts, PROB, outcome, node="ES.centre", statement=reason, reason=json.dumps(meas),
                   survived=(outcome == "accepted"))
        out.append(dict(composite=c.id, parts=list(c.parts), outcome=outcome, reason=reason, measurement=meas))
    mem.commit()
    acc = sorted((r for r in out if r["outcome"] == "accepted"), key=lambda r: -r["measurement"]["c"])
    summary = dict(engine=f"Anaxia (ebllm) {ebllm.__version__}", library=len(LIB), closed=len(out),
                   accepted=len(acc), rejected=sum(r["outcome"] == "rejected" for r in out),
                   conditional=sum(r["outcome"] == "conditional" for r in out),
                   best=acc[0]["measurement"] if acc else None, memory=mem.stats(), pmax=PMAX, nrule=NRULE)
    summary["accepted_constants"] = sorted({round(r["measurement"]["c"], 6) for r in acc})
    json.dump(dict(summary=summary, composites=out), open("anaxia_centre.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
    best = summary["best"]
    assert best["form"] == "farey" and best["rule"] == "chain" and best["sparse"] == "family"
    assert 0.67920 < best["c"] < 0.67922
    print("OK: best accepted composite is the corners with the strips, the Farey form, the chain rule and "
          "the coincidence family, c =", best["c"])


if __name__ == "__main__":
    main()
