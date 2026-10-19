"""check.py -- independent check in pure Python (standard library only).

1. Exact check, with fractions, of the inversion h0 = k0 * (1 (x) 1) for u, a <= 120.
2. Exact check of the local factors: 1 + 2/(p(p-1)) - (p+1)/((p-1)p^2) = 1 + 1/p^2.
3. The partial Euler product of (1 + p^-2) over odd p <= 10^6 against
   (zeta(2)/zeta(4))/(5/4) = 12/pi^2, so the constant is (1/2)(3/4)(12/pi^2) = 9/(2 pi^2).
4. W(J) for J <= 3000 against K((log J + gamma + T)^2 + M).
"""
import math
from fractions import Fraction as F


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i :: i] = bytearray(len(s[i * i :: i]))
    return [i for i in range(n + 1) if s[i]]


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def k0(d, e):
    fd, fe = factor(d), factor(e)
    r = F(1)
    for p in set(fd) | set(fe):
        i, j = fd.get(p, 0), fe.get(p, 0)
        if p == 2:
            v = {(0, 0): 1, (1, 1): -1}.get((i, j), 0)
        else:
            v = {(0, 0): F(1), (1, 0): F(1, p - 1), (0, 1): F(1, p - 1), (1, 1): F(-(p + 1), p - 1)}.get((i, j), 0)
        r *= v
        if r == 0:
            return F(0)
    return r


def main():
    ok = True
    N = 120
    divs = [[] for _ in range(N + 1)]
    for d in range(1, N + 1):
        for m in range(d, N + 1, d):
            divs[m].append(d)
    kt = {(d, e): k0(d, e) for d in range(1, N + 1) for e in range(1, N + 1)}
    for u in range(1, N + 1):
        for a in range(1, N + 1):
            s = sum(kt[d, e] for d in divs[u] for e in divs[a])
            h = F(2 * u * a, phi(4 * u * a)) if math.gcd(u, a) == 1 else F(0)
            if s != h:
                print("inversion mismatch", u, a)
                ok = False
    print("1. inversion checked for u, a <=", N)

    for p in primes_upto(10000)[1:]:
        lhs = 1 + F(2, p * (p - 1)) - F(p + 1, (p - 1) * p * p)
        if lhs != 1 + F(1, p * p):
            print("local factor mismatch", p)
            ok = False
    print("2. local factor 1 + 1/p^2 checked exactly for odd p < 10^4")

    ps = primes_upto(10 ** 6)
    prod = 1.0
    for p in ps[1:]:
        prod *= 1 + 1 / (p * p)
    target = 12 / math.pi ** 2
    print("3. prod_{2<p<=10^6} (1+p^-2) = %.10f, 12/pi^2 = %.10f" % (prod, target))
    if abs(prod - target) > 1e-6:
        ok = False

    K = 9 / (2 * math.pi ** 2)
    g = 0.57721566490153286
    T = math.log(2) / 3 + sum(math.log(p) / ((p - 1) * (p * p + 1)) for p in ps[1:])
    M = -4 / 9 * math.log(2) ** 2 - sum(
        math.log(p) ** 2 * ((p + 1) / ((p - 1) * (p * p + 1)) + 1 / ((p - 1) ** 2 * (p * p + 1) ** 2)) for p in ps[1:]
    )
    M -= (math.log(10 ** 6) + 1) / 10 ** 6  # primes above 10^6, by the prime number theorem
    print("   K = %.12f, T = %.10f, M = %.10f" % (K, T, M))
    if abs(M + 1.0429405079) > 1e-7:
        ok = False
    J = 3000
    ph = list(range(J + 1))
    for p in range(2, J + 1):
        if ph[p] == p:
            for k in range(p, J + 1, p):
                ph[k] -= ph[k] // p
    W = 0.0
    for m in range(1, J + 1):
        for x in range(1, m + 1):
            if math.gcd(m, x) != 1:
                continue
            t = 1 / ((4 if (m * x) % 2 == 0 else 2) * ph[m] * ph[x])
            W += t if x == m else 2 * t
        if m in (100, 1000, 3000):
            L = math.log(m)
            f = K * ((L + g + T) ** 2 + M)
            print("4. J = %4d  W = %.9f  formula = %.9f  error = %+.2e" % (m, W, f, W - f))
            if abs(W - f) > 3 * L / math.sqrt(m):
                ok = False
    print("Python check passed" if ok else "Python check FAILED")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
