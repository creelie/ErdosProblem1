"""Witness regions for the Erdos-Straus large sieve, assembled with the typed crossing and the memory
of the engine Anaxia (creelie/eBLLM, version 3.0.5).

A Type II triple for a prime p = 3 (mod 4) is (u, a, alpha) with u a alpha = P = (p + 1)/4 and
gcd(u, a) = 1; every n = -u/a (mod p) solves 4/n = 1/x + 1/y + 1/z.  In the coordinates
(log u, log a, log alpha) the triples fill a triangle.  A witness region is a corner of it in which one
product of two variables is at most sqrt(p)/16 (so that the primes can be counted with moduli of
level 1/2), together with a normal form of the class that is meant to prove the classes distinct:

    farey  u/a                 distinct when the products u a' stay below p,
    M      4 u^2 alpha         distinct when u^2 alpha < p,
    I      (4 a^2 alpha)^-1    distinct when a^2 alpha < p.

The library holds the three corners with each normal form, an empty corner, the centre of the
triangle (counted only under the Elliott-Halberstam conjecture), the three normal forms, the
counting inputs and the large sieve with the tilted lower bound.  The engine crosses typed nodes
until the sieve is closed; the check below decides every closed composite on the primes up to
PMAX and returns the constant c of the bound N exp(-(c - eps)(log N)^(2/3)).  Every verdict goes
to the engine's memory, which is exported to anaxia_regions.json.

Run:  pip install "git+https://github.com/creelie/eBLLM@v3.0.5"; python3 anaxia_regions.py
"""
import json
import math
import sys

from ebllm.bllm.targets import Con, cross
from ebllm.bllm.memory import Memory
import ebllm

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
KAPPA = 9 / (2 * math.pi ** 2)
PAIRS = ("ua", "ualpha", "aalpha")
FORMS = ("farey", "M", "I")
PROB = "erdos-straus"


def node(i, gives, needs=(), source=""):
    return Con(i, PROB, gives, frozenset(needs), source, (i,))


LIB = [node("es.sieve", "bound", [f"region:{q}" for q in PAIRS] + ["region:centre"],
            "Montgomery's large sieve with Q = sqrt(N) and the tilted lower bound for L(Q)")]
for q in PAIRS:
    for f in FORMS:
        LIB.append(node(f"es.{q}.{f}", f"region:{q}", [f"count:{q}", f"form:{f}"],
                        f"triples with 16*{q} <= sqrt(p), classes compared in the {f} form"))
    LIB.append(node(f"es.{q}.none", f"region:{q}", [], "this corner left out"))
    LIB.append(node(f"es.bv.{q}", f"count:{q}", [], f"Bombieri-Vinogradov for the moduli 4*{q} <= p^(1/2-delta)"))
LIB.append(node("es.centre.none", "region:centre", [], "the centre left out"))
LIB.append(node("es.centre.eh", "region:centre", ["count:centre"], "the centre of the triangle"))
LIB.append(node("es.eh", "count:centre", [], "Elliott-Halberstam conjecture (moduli up to p^(2/3))"))
for f in FORMS:
    LIB.append(node(f"es.form.{f}", f"form:{f}", [], f"the {f} normal form of the class -u/a"))


def breed(lib):
    """Typed crossing to a fixed point: feed library nodes into every composite whose input they fill.
    Composites with the same parts are the same construction, whatever the order of the crossings."""
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


# -- the check ------------------------------------------------------------------------------------
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


def corner(u, a, al, p):
    r = math.isqrt(p)
    if 16 * u * a <= r and u >= 2 and a >= 2:
        return "ua"
    if 16 * u * al <= r:
        return "ualpha"
    if 16 * a * al <= r:
        return "aalpha"
    return "centre"


HEIGHT = {"farey": lambda u, a, al: (u, a), "M": lambda u, a, al: u * u * al, "I": lambda u, a, al: a * a * al}


def certified(form, trip, p):
    """Does the normal form prove the classes of these triples pairwise distinct modulo p?"""
    if form == "farey":
        mu = max((t[0] for t in trip), default=0)
        ma = max((t[1] for t in trip), default=0)
        return mu * ma < p
    return all(HEIGHT[form](*t) < p for t in trip)


DATA = None


def data():
    global DATA
    if DATA is None:
        DATA = []
        for p in primes_3mod4(PMAX):
            P = (p + 1) // 4
            trip = {}
            for u in divisors(P):
                for a in divisors(P // u):
                    if math.gcd(u, a) == 1:
                        al = P // (u * a)
                        trip.setdefault(corner(u, a, al, p), []).append((u, a, al))
            DATA.append((p, trip))
    return DATA


def check(c):
    chosen = {}
    for part in c.parts:
        bits = part.split(".")
        if len(bits) == 3 and bits[1] in PAIRS + ("centre",) and bits[2] in FORMS + ("none", "eh"):
            chosen[bits[1]] = bits[2]
    used = [q for q in PAIRS if chosen.get(q, "none") != "none"]
    conditional = chosen.get("centre") == "eh"
    fails = {q: 0 for q in used}
    col_inner = col_edge = 0
    s_trip = 0.0
    for p, trip in data():
        for q in used:
            if not certified(chosen[q], trip.get(q, []), p):
                fails[q] += 1
        cls = {}
        for q in used + (["centre"] if conditional else []):
            for t in trip.get(q, []):
                u, a, al = t
                cls.setdefault(u * pow(a, -1, p) % p, []).append(t)
                s_trip += 1 / p
        for v in cls.values():
            if len(v) > 1:
                if any(t[0] == 1 or t[1] == 1 for t in v):
                    col_edge += len(v) - 1
                else:
                    col_inner += len(v) - 1
    area = len(used) / 4 + (0.25 if conditional else 0.0)
    a_const = KAPPA * area / 4 if area else 0.0
    cconst = 1.5 * a_const ** (1 / 3) if area else 0.0
    meas = dict(corners=used, forms={q: chosen[q] for q in used}, centre=chosen.get("centre", "none"),
                uncertified_primes=fails, collisions_inner=col_inner, collisions_edge=col_edge,
                area=area, a=a_const, c=cconst, S_triples=round(s_trip, 6), pmax=PMAX)
    if not used and not conditional:
        return "rejected", "no corner: nothing is sieved", meas
    if any(fails.values()):
        bad = [f"{q}/{chosen[q]} on {fails[q]} primes" for q in used if fails[q]]
        return "rejected", "normal form does not prove distinctness: " + ", ".join(bad), meas
    if col_inner:
        return "rejected", f"{col_inner} collisions between corners", meas
    if conditional:
        return "conditional", f"c = {cconst:.6f} if the Elliott-Halberstam conjecture holds", meas
    return "accepted", f"c = {cconst:.6f}", meas


def main():
    mem = Memory()
    mem.add_nodes(LIB)
    closed = breed(LIB)
    out = []
    for c in sorted(closed, key=lambda c: c.id):
        outcome, reason, meas = check(c)
        mem.record(c.parts, PROB, outcome, node="ES.exceptional", statement=reason, reason=json.dumps(meas),
                   survived=(outcome == "accepted"))
        out.append(dict(composite=c.id, parts=list(c.parts), outcome=outcome, reason=reason, measurement=meas))
    mem.commit()
    acc = sorted((r for r in out if r["outcome"] == "accepted"), key=lambda r: -r["measurement"]["c"])
    summary = dict(engine=f"Anaxia (ebllm) {ebllm.__version__}", library=len(LIB), closed=len(out),
                   accepted=len(acc), rejected=sum(r["outcome"] == "rejected" for r in out),
                   conditional=sum(r["outcome"] == "conditional" for r in out),
                   best=acc[0]["measurement"] if acc else None, memory=mem.stats(), pmax=PMAX)
    ladder = sorted({round(r["measurement"]["c"], 6) for r in acc})
    summary["accepted_constants"] = ladder
    json.dump(dict(summary=summary, composites=out), open("anaxia_regions.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
    best = summary["best"]
    assert best["forms"] == {"ua": "farey", "ualpha": "M", "aalpha": "I"} and best["collisions_inner"] == 0
    assert abs(best["c"] - 9 / (2 ** (8 / 3) * math.pi ** (2 / 3))) < 1e-12
    print("OK: best accepted composite is the three corners with the Farey, M and I forms, c =", best["c"])


if __name__ == "__main__":
    main()
