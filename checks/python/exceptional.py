"""Checks for "An explicit exponent for the exceptional set of the Erdos-Straus equation".

Standard library only.  Run: python3 -I checks/python/exceptional.py

1. The Type II identity 4/n = 1/(alpha u v) + 1/(a alpha n u) + 1/(a alpha n v) for every n <= 3000 in
   the class -u/a modulo p = 4 u a alpha - 1, for every triple with p prime and p < 400.
2. The normal forms: u/a = 4 u^2 alpha = (4 a^2 alpha)^(-1) modulo p for every triple and every prime
   p = 3 (mod 4) below 30000.
3. Coincidences between the corners B and C built from the parametrisation: for every
   k, m, r <= 12 and d1 e1 dividing k r, the roots tau of k tau^2 - 4 k m r tau + (k r/(d1 e1) - m^2) = 0
   give P = G d1 e1 W, d = G d1, e = G e1 and M; every split u | M with u | d gives two triples whose
   classes agree modulo p = 4P - 1, with P^2 = d e (M + r p) and P >= 1024 k^2 m^2 r^5 d1 e1.
4. The coincidence inside the centre at the prime 531359.
5. The constants: kappa, the Euler factors of the corner B, Z = 9/pi^2, the exponent
   c0 = 9/(2^(8/3) pi^(2/3)) and the ladder 0.458165, 0.577252, 0.660789, 0.727292.
6. The level (Remark VII.2): for the primes p = 3 (mod 4) below 10^6, no two triples in the corners
   at level isqrt(p)//4 share a class; at level isqrt(p)//2 exactly 17568 pairs do, all between B and C,
   the first at p = 31; the only coincidences between triples with u, a >= 2 outside the corners at
   level sqrt(p)/16 are the pair at 531359 and its mirror image.
"""
from fractions import Fraction as Fr
from math import gcd, isqrt, pi, log

FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL:", msg)


def is_prime(n):
    if n < 2:
        return False
    for q in range(2, isqrt(n) + 1):
        if n % q == 0:
            return False
    return True


def divisors(n):
    small = [d for d in range(1, isqrt(n) + 1) if n % d == 0]
    return sorted(set(small + [n // d for d in small]))


def triples(P):
    for u in divisors(P):
        for a in divisors(P // u):
            if gcd(u, a) == 1:
                yield u, a, P // (u * a)


# 1. the identity
count = 0
for p in range(3, 400, 4):
    if not is_prime(p):
        continue
    P = (p + 1) // 4
    for u, a, al in triples(P):
        r = (-u * pow(a, -1, p)) % p
        for n in range(r if r else p, 3001, p):
            v, rem = divmod(a * n + u, p)
            check(rem == 0, f"p | a n + u fails at {(p, u, a, al, n)}")
            lhs = Fr(4, n)
            rhs = Fr(1, al * u * v) + Fr(1, a * al * n * u) + Fr(1, a * al * n * v)
            check(lhs == rhs, f"identity fails at {(p, u, a, al, n)}")
            count += 1
print(f"1. Type II identity: {count} cases checked")

# 2. normal forms
count = 0
for p in range(3, 30000, 4):
    if not is_prime(p):
        continue
    P = (p + 1) // 4
    for u, a, al in triples(P):
        r = u * pow(a, -1, p) % p
        check(r == 4 * u * u * al % p, f"M-form fails at {(p, u, a, al)}")
        check(r == pow(4 * a * a * al, -1, p), f"I-form fails at {(p, u, a, al)}")
        count += 1
print(f"2. normal forms: {count} triples checked")

# 3. the parametrisation of the coincidences between B and C
count = 0
for k in range(1, 13):
    for m in range(1, 13):
        for r in range(1, 13):
            for d1 in divisors(k * r):
                for e1 in divisors(k * r // d1):
                    if gcd(d1, e1) != 1:
                        continue
                    c = Fr(k * r, d1 * e1) - m * m
                    disc = (4 * k * m * r) ** 2 - 4 * k * c
                    if disc < 0 or disc.denominator != 1 or isqrt(int(disc)) ** 2 != disc:
                        continue
                    for sgn in (1, -1):
                        tau = Fr(4 * k * m * r + sgn * isqrt(int(disc)), 2 * k)
                        if tau.denominator != 1 or tau < 1:
                            continue
                        tau = int(tau)
                        G = 4 * k * tau * r + m
                        W = 4 * G * r + tau
                        P = G * d1 * e1 * W
                        d, e = G * d1, G * e1
                        M = d1 * e1 * tau * W + r
                        if (G * G * d1 * e1) % M:
                            continue
                        p = 4 * P - 1
                        check(P * P == d * e * (M + r * p), f"P^2 = de(M + rp) fails at {(k, m, r, d1, e1)}")
                        check(tau >= 4 * m * r, f"tau < 4mr at {(k, m, r, d1, e1)}")
                        check(P >= 1024 * k * k * m * m * r ** 5 * d1 * e1, f"P bound fails at {(k, m, r)}")
                        for u in divisors(M):
                            if d % u or e % (M // u):
                                continue
                            ap = M // u
                            al, alp = d // u, e // ap
                            if al * alp != k or P % (u * al) or P % (ap * alp):
                                continue
                            a, up = P // (u * al), P // (ap * alp)
                            if gcd(u, a) != 1 or gcd(up, ap) != 1:
                                continue
                            check((u * ap - a * up) % p == 0, f"classes differ at {(P, u, a, up, ap)}")
                            count += 1
print(f"3. coincidences built from the parametrisation: {count}, all with equal classes")

# 4. the centre
p = 531359
check(is_prime(p) and (p + 1) // 4 == 81 * 328 * 5 == 1640 * 81, "centre example")
check((81 * pow(328, -1, p) - 1640 * pow(81, -1, p)) % p == 0, "centre coincidence")
X = isqrt(p) / 16
check(all(x > X for x in (81 * 328, 81 * 5, 328 * 5, 1640 * 81, 1640, 81)), "centre example not in the centre")
print("4. 81/328 = 1640/81 modulo the prime 531359; both triples lie in the centre")

# 5. constants
kappa = 9 / (2 * pi ** 2)
check(abs(kappa - 0.4559453263) < 1e-10, "kappa")
for q in (3, 5, 7, 11, 101):
    q = Fr(q)
    lhs = Fr(1, 1) / (q - 1) * q / (q - 1) + q / (q - 1) * Fr(1, 1) / (q - 1) + 1
    check(lhs == (q * q + 1) / (q - 1) ** 2, f"local factor of B at {q}")
    check(lhs / (q * q / (q - 1) ** 2) == (q * q + 1) / (q * q), f"ratio at {q}")
check(Fr(1, 2) * 1 * 2 + 2 == Fr(3, 4) * 4, "local factor at 2")
Z = 0.75
for q in range(3, 3000, 2):
    if is_prime(q):
        Z *= 1 + 1 / q ** 2
check(abs(Z - 9 / pi ** 2) < 1e-3, f"Z = {Z} against 9/pi^2")
c0 = 9 / (2 ** (8 / 3) * pi ** (2 / 3))
check(abs(c0 - 1.5 * (27 / (32 * pi ** 2)) ** (1 / 3)) < 1e-14, "c0 closed form")
check(abs(c0 ** 3 - 729 / (256 * pi ** 2)) < 1e-14, "c0 cube")
ladder = [1.5 * (j * kappa / 16) ** (1 / 3) for j in (1, 2, 3, 4)]
for got, want in zip(ladder, (0.458165, 0.577252, 0.660789, 0.727292)):
    check(abs(got - want) < 1e-6, f"ladder {got} against {want}")
print(f"5. kappa = {kappa:.10f}, Z = 9/pi^2 = {9 / pi ** 2:.10f}, c0 = {c0:.16f}")
print("   ladder: " + ", ".join(f"{x:.6f}" for x in ladder))

# 6. the level
N = 10 ** 6
sieve = bytearray([1]) * (N + 1)
sieve[0] = sieve[1] = 0
for i in range(2, isqrt(N) + 1):
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(sieve[i * i::i]))


def regions(t, X):
    u, a, al = t
    return {k for k, ok in (("A", u >= 2 and a >= 2 and u * a <= X), ("B", u * al <= X), ("C", a * al <= X)) if ok}


def outside(t, p):
    u, a, al = t
    r = isqrt(p)
    return u >= 2 and a >= 2 and 16 * u * a > r and 16 * u * al > r and 16 * a * al > r


pairs = {4: [], 2: []}
centre = []
for p in range(3, N, 4):
    if not sieve[p]:
        continue
    cls = {}
    for t in triples((p + 1) // 4):
        cls.setdefault(t[0] * pow(t[1], -1, p) % p, []).append(t)
    for v in cls.values():
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                t, s2 = v[i], v[j]
                for f in pairs:
                    X = isqrt(p) // f
                    rt, rs = regions(t, X), regions(s2, X)
                    if rt and rs:
                        pairs[f].append((p, t, s2, rt, rs))
                if outside(t, p) and outside(s2, p):
                    centre.append((p, t, s2))
check(pairs[4] == [], f"coincidences at level sqrt(p)/4: {pairs[4][:3]}")
check(len(pairs[2]) == 17568, f"{len(pairs[2])} pairs at level sqrt(p)/2")
check(all({"B", "C"} == rt | rs and len(rt) == len(rs) == 1 for _, _, _, rt, rs in pairs[2]), "B-C only")
check(pairs[2][0][:3] == (31, (1, 4, 2), (8, 1, 1)), f"first pair {pairs[2][0][:3]}")
check(sorted(centre) == sorted([(531359, (81, 328, 5), (1640, 81, 1)), (531359, (81, 1640, 1), (328, 81, 5))]),
      f"centre coincidences {centre}")
print(f"6. level: none at sqrt(p)/4, {len(pairs[2])} B-C pairs at sqrt(p)/2 from p = 31, "
      f"{len(centre)} outside the corners (both at 531359)")

if FAILS:
    raise SystemExit(f"{len(FAILS)} failures")
print("exceptional.py: all checks passed")
