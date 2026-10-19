/* exceptional.c -- checks for "An explicit exponent for the exceptional set of the Erdos-Straus
   equation".
   Part 1. For every prime p = 3 (mod 4) up to Y, with P = (p + 1)/4 and the level X = sqrt(p)/16,
   the Type II triples (u, a, alpha), u a alpha = P, gcd(u, a) = 1, are sorted into the corners
     A: u, a >= 2 and u a <= X,   B: u alpha <= X,   C: a alpha <= X.
   It checks that the corners are disjoint, that the class u/a (mod p) equals 4 u^2 alpha and
   (4 a^2 alpha)^(-1), that the classes are distinct inside each corner, that no class of A is a
   class of B or C, and it counts the coincidences between B and C.  It compares
   sum_p |K(p)|/p with the prediction sum_{all primes p <= Y} V_K(X_p)/p for K = A, B, C, where
   V_A(X) = sum_{ua <= X, (u,a) = 1, u,a >= 2} 1/phi(4ua) and V_B(X) = sum_{u alpha <= X} (phi(u)/u)/phi(4 u alpha).
   Part 2. For every P <= PB and p = 4P - 1 (prime or not), every pair of a triple (u, a, alpha) and
   a triple (u', a', alpha') with coprime pairs, d = u alpha, e = a' alpha', m = u d, i = a' e, m i < P^2
   and 16 m i = 1 (mod p): it checks that P^2 = d e (M + r p) with M = u a' and r >= 1, that
   t = P - 4 d e r is never negative, that t = 0 forces u = a' = 1 and P = 4 d e, and that for t > 0 the
   parametrisation of the paper holds, with tau >= 4 m r and P >= 1024 k^2 m^2 r^5 d1 e1.
   Usage: exceptional [Y] [PB]   (defaults 10^7 and 4000). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned long long u64;
typedef __uint128_t u128;
static u64 gcd(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }
static u64 pw(u64 b, u64 e, u64 m) { u64 r = 1; b %= m; while (e) { if (e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; } return r; }
static int cmp(const void *x, const void *y) { u64 a = *(const u64 *)x, b = *(const u64 *)y; return a < b ? -1 : a > b; }
static long fails = 0;
#define CHECK(c, ...) do { if (!(c)) { fails++; if (fails < 20) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

static int divisors(u64 n, u64 *dv) { int k = 0; for (u64 d = 1; d * d <= n; d++) if (n % d == 0) { dv[k++] = d; if (d * d != n) dv[k++] = n / d; } return k; }
static int in_level(u64 x, u64 p) { return (u128)256 * x * x <= p; } /* 16 x <= sqrt(p) */

int main(int argc, char **argv) {
  long Y = argc > 1 ? atol(argv[1]) : 10000000L;
  long PB = argc > 2 ? atol(argv[2]) : 4000L;
  /* ---- Part 1 ---- */
  char *comp = calloc(Y + 2, 1);
  for (long i = 2; i * i <= Y + 1; i++) if (!comp[i]) for (long j = i * i; j <= Y + 1; j += i) comp[j] = 1;
  long XM = (long)(sqrt((double)Y) / 16) + 2, FM = 4 * XM + 4;
  long *ph = malloc((FM + 1) * sizeof(long));
  for (long i = 0; i <= FM; i++) ph[i] = i;
  for (long i = 2; i <= FM; i++) if (ph[i] == i) for (long j = i; j <= FM; j += i) ph[j] -= ph[j] / i;
  double *VA = calloc(XM + 1, sizeof(double)), *VB = calloc(XM + 1, sizeof(double));
  for (long u = 1; u <= XM; u++) for (long m = u; m <= XM; m += u) { long a = m / u;
      if (u >= 2 && a >= 2 && gcd(u, a) == 1) VA[m] += 1.0 / ph[4 * m];
      VB[m] += (double)ph[u] / u / ph[4 * m]; }
  for (long m = 1; m <= XM; m++) { VA[m] += VA[m - 1]; VB[m] += VB[m - 1]; }
  static u64 dv[20000]; static u64 cls[3][200000];
  double sA = 0, sB = 0, sC = 0, sD = 0, pA = 0, pB = 0; long nBC = 0, nprimes = 0, ntrip = 0;
  for (long p = 3; p <= Y; p += 2) {
    if (comp[p]) continue;
    long X = (long)floor(sqrt((double)p) / 16); while ((u128)256 * (X + 1) * (X + 1) <= (u64)p) X++; while (X > 0 && !in_level(X, p)) X--;
    pA += VA[X] / p; pB += VB[X] / p;
    if (p % 4 != 3) continue;
    nprimes++;
    u64 P = (p + 1) / 4; int nd = divisors(P, dv); long n[3] = {0, 0, 0};
    for (int i = 0; i < nd; i++) for (int j = 0; j < nd; j++) {
      u64 u = dv[i], a = dv[j];
      if ((P / u) % a || gcd(u, a) != 1) continue;
      u64 al = P / u / a; ntrip++;
      int inA = u >= 2 && a >= 2 && in_level(u * a, p), inB = in_level(u * al, p), inC = in_level(a * al, p);
      CHECK(inA + inB + inC <= 1, "corners overlap at p=%ld (%llu,%llu,%llu)", p, u, a, al);
      u64 r = (u128)u * pw(a, p - 2, p) % p;
      u64 rM = (u128)4 * ((u128)u * u % p) % p * al % p;
      u64 rI = pw((u128)4 * ((u128)a * a % p) % p * al % p, p - 2, p);
      CHECK(r == rM && r == rI, "normal forms differ at p=%ld (%llu,%llu,%llu)", p, u, a, al);
      int k = inA ? 0 : inB ? 1 : inC ? 2 : -1;
      if (k >= 0) cls[k][n[k]++] = r;
    }
    for (int k = 0; k < 3; k++) qsort(cls[k], n[k], sizeof(u64), cmp);
    for (int k = 0; k < 3; k++) for (long t = 1; t < n[k]; t++) CHECK(cls[k][t] != cls[k][t - 1], "repeated class in corner %d at p=%ld", k, p);
    /* A against B and C; B against C, by merging sorted lists */
    long dist = n[0] + n[1] + n[2];
    for (int k = 1; k < 3; k++) { long s = 0, t = 0; while (s < n[0] && t < n[k]) { if (cls[0][s] == cls[k][t]) { CHECK(0, "A meets corner %d at p=%ld", k, p); s++; t++; } else if (cls[0][s] < cls[k][t]) s++; else t++; } }
    { long s = 0, t = 0; while (s < n[1] && t < n[2]) { if (cls[1][s] == cls[2][t]) { nBC++; dist--; s++; t++; } else if (cls[1][s] < cls[2][t]) s++; else t++; } }
    sA += (double)n[0] / p; sB += (double)n[1] / p; sC += (double)n[2] / p; sD += (double)dist / p;
  }
  printf("Part 1: primes p = 3 (mod 4) up to %ld: %ld primes, %ld Type II triples\n", Y, nprimes, ntrip);
  printf("  corners disjoint, normal forms equal, classes distinct in each corner, no class of A in B or C: %s\n", fails ? "NO" : "yes");
  printf("  coincidences between B and C: %ld\n", nBC);
  printf("  sum |A|/p = %.5f   prediction sum V_A(X_p)/p = %.5f   ratio %.4f\n", sA, pA, sA / pA);
  printf("  sum |B|/p = %.5f   prediction sum V_B(X_p)/p = %.5f   ratio %.4f\n", sB, pB, sB / pB);
  printf("  sum |C|/p = %.5f   prediction sum V_B(X_p)/p = %.5f   ratio %.4f\n", sC, pB, sC / pB);
  printf("  sum of distinct classes / p = %.5f\n", sD);
  CHECK(fabs(sA / pA - 1) < 0.05 && fabs(sB / pB - 1) < 0.05 && fabs(sC / pB - 1) < 0.05, "corner sums far from the prediction");
  /* ---- Part 2 ---- */
  long npairs = 0, nt0 = 0, ntpos = 0;
  for (u64 P = 1; P <= (u64)PB; P++) {
    u64 p = 4 * P - 1; int nd = divisors(P, dv);
    static u64 tu[20000], ta[20000], tl[20000]; int nt = 0;
    for (int i = 0; i < nd; i++) for (int j = 0; j < nd; j++) { u64 u = dv[i], a = dv[j];
      if ((P / u) % a || gcd(u, a) != 1) continue; tu[nt] = u; ta[nt] = a; tl[nt] = P / u / a; nt++; }
    for (int x = 0; x < nt; x++) for (int y = 0; y < nt; y++) {
      u64 u = tu[x], al = tl[x], ap = ta[y], alp = tl[y];
      u64 d = u * al, e = ap * alp, m = u * d, i = ap * e;
      if ((u128)m * i >= (u128)P * P) continue;
      if ((u128)16 * m % p * i % p != 1 % p) continue;
      npairs++;
      u64 M = u * ap; u128 diff = (u128)P * P - (u128)M * d * e;
      CHECK(diff % ((u128)d * e) == 0 && (diff / ((u128)d * e)) % p == 0, "P^2 - M d e not divisible by d e p at P=%llu", P);
      u64 r = (u64)(diff / ((u128)d * e) / p);
      CHECK(r >= 1, "r < 1 at P=%llu", P);
      long long t = (long long)P - 4LL * (long long)(d * e * r);
      CHECK(t >= 0, "t < 0 at P=%llu", P);
      if (t == 0) { nt0++; CHECK(u == 1 && ap == 1 && P == 4 * d * e, "t = 0 without u = a' = 1 at P=%llu", P); continue; }
      ntpos++;
      u64 G = gcd(d, e), d1 = d / G, e1 = e / G;
      CHECK(t % (long long)(G * d1 * e1) == 0, "lcm(d,e) does not divide t at P=%llu", P);
      u64 tau = (u64)t / (G * d1 * e1), W = 4 * G * r + tau;
      CHECK(P == G * d1 * e1 * W && M == d1 * e1 * tau * W + r, "parametrisation fails at P=%llu", P);
      u64 k = al * alp;
      CHECK((u128)k * M == (u128)G * G * d1 * e1 && (k * r) % (d1 * e1) == 0, "k M = G^2 d1 e1 fails at P=%llu", P);
      CHECK(G > 4 * k * tau * r, "G <= 4 k tau r at P=%llu", P);
      u64 mm = G - 4 * k * tau * r;
      CHECK(tau >= 4 * mm * r && (u128)P >= (u128)1024 * k * k * mm * mm * r * r * r * r * r * d1 * e1, "bounds fail at P=%llu", P);
    }
  }
  printf("Part 2: P <= %ld: %ld coincidences (t = 0: %ld, t > 0: %ld), parametrisation and bounds hold: %s\n",
         PB, npairs, nt0, ntpos, fails ? "NO" : "yes");
  if (fails) { printf("%ld failures\n", fails); return 1; }
  printf("exceptional.c: all checks passed\n");
  return 0;
}
