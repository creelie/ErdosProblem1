/* centre.c -- checks for the centre of the triangle in "An explicit exponent for the exceptional
   set of the Erdos-Straus equation".
   For every prime p = 3 (mod 4) up to Y, with P = (p + 1)/4, the Type II triples (u, a, alpha),
   u a alpha = P, gcd(u, a) = 1, are sorted into
     the centre:   u^2 < P, a^2 < P, alpha^2 < P,
     the corners at the level X = sqrt(p)/4 (16 X^2 <= p):
       A: u, a >= 2 and u a <= X,   B: u alpha <= X,   C: a alpha <= X.
   It checks that
     (1) the classes u/a (mod p) of the centre are pairwise distinct;
     (2) no class of the centre is a class of A;
     (3) no centre triple with u < a shares its class with B, and none with u > a with C;
     (4) every coincidence of a centre triple t = (u, a, alpha) with u < a and a triple
         t'' = (u'', a'', alpha'') of C satisfies, with s = a'' alpha'' and u'' a - u a'' = k p,
           k >= 1,  j = a - 4 k s >= 1,  D = s a'' - a alpha j >= 1,
           P j = s (u a'' - k),  u D = s k,  alpha j (j + 4 k s) = s a'' - D,
         and that t and t'' are recovered from (a'', alpha'', k, D, j); the coincidences of a centre
         triple with u > a and a triple of B are the images of these under (u, a, alpha) -> (a, u, alpha);
     (5) the number of these coincidences, nu'(p), has a small weighted sum sum nu'(p)/p.
   It also counts the triples in the three strips of depth < 1/42 used in the paper.
   Usage: centre [Y]   (default 10^7). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned long long u64;
typedef __uint128_t u128;
static u64 gcd(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }
static u64 pw(u64 b, u64 e, u64 m) { u64 r = 1; b %= m; while (e) { if (e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; } return r; }
static long fails = 0;
#define CHECK(c, ...) do { if (!(c)) { fails++; if (fails < 20) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

typedef struct { u64 c, u, a, al; int reg; } item; /* reg: 0 centre, 1 A, 2 B, 3 C, 4 other */
static int cmp(const void *x, const void *y) { u64 a = ((const item *)x)->c, b = ((const item *)y)->c; return a < b ? -1 : a > b; }
static int divisors(u64 n, u64 *dv) { int k = 0; for (u64 d = 1; d * d <= n; d++) if (n % d == 0) { dv[k++] = d; if (d * d != n) dv[k++] = n / d; } return k; }
static int in_level(u64 x, u64 p) { return (u128)16 * x * x <= p; } /* 4 x <= sqrt(p) */

/* the coincidence of a centre triple t (u < a) with a triple t2 of C */
static void family(u64 p, u64 P, const item *t, const item *t2, long *nfam, int print) {
  u64 u = t->u, a = t->a, al = t->al, u2 = t2->u, a2 = t2->a, al2 = t2->al;
  u64 s = a2 * al2;
  u128 lhs = (u128)u2 * a, rhs = (u128)u * a2;
  CHECK(lhs > rhs && (lhs - rhs) % p == 0, "u''a - ua'' not a positive multiple of p at p=%llu", p);
  u64 k = (u64)((lhs - rhs) / p);
  CHECK(k >= 1 && 4 * k * s < a, "k or j out of range at p=%llu", p);
  u64 j = a - 4 * k * s;
  u128 sa2 = (u128)s * a2, aaj = (u128)a * al * j;
  CHECK(sa2 > aaj, "D <= 0 at p=%llu", p);
  u64 D = (u64)(sa2 - aaj);
  CHECK((u128)P * j == (u128)s * ((u128)u * a2 - k), "P j = s(u a'' - k) fails at p=%llu", p);
  CHECK((u128)u * D == (u128)s * k, "u D = s k fails at p=%llu", p);
  CHECK((u128)al * j * (j + 4 * k * s) == sa2 - D, "alpha j (j + 4ks) = s a'' - D fails at p=%llu", p);
  /* recovery from (a'', alpha'', k, D, j) */
  u64 s_ = a2 * al2, u_ = s_ * k / D, a_ = j + 4 * k * s_;
  u64 al_ = (u64)(((u128)s_ * a2 - D) / ((u128)j * a_));
  u64 P_ = u_ * a_ * al_, u2_ = P_ / s_;
  CHECK(u_ == u && a_ == a && al_ == al && P_ == P && u2_ == u2, "recovery fails at p=%llu", p);
  if (print)
    printf("  p=%llu: (%llu,%llu,%llu) and (%llu,%llu,%llu), (a'',alpha'',k,D,j)=(%llu,%llu,%llu,%llu,%llu)\n",
           p, u, a, al, u2, a2, al2, a2, al2, k, D, j);
  (*nfam)++;
}

int main(int argc, char **argv) {
  long Y = argc > 1 ? atol(argv[1]) : 10000000L;
  char *comp = calloc(Y + 2, 1);
  for (long i = 2; i * i <= Y + 1; i++) if (!comp[i]) for (long j = i * i; j <= Y + 1; j += i) comp[j] = 1;
  static u64 dv[20000]; static item it[400000];
  long nprimes = 0, ncentre = 0, nstrip = 0, nfam = 0, nfam_mirror = 0, nfirst = 0;
  double sfam = 0, sstrip = 0;
  u64 first_p = 0;
  for (long p = 3; p <= Y; p += 4) {
    if (comp[p]) continue;
    nprimes++;
    u64 P = (p + 1) / 4; int nd = divisors(P, dv), n = 0;
    double lP = log((double)P);
    for (int i = 0; i < nd; i++) for (int jj = 0; jj < nd; jj++) {
      u64 u = dv[i], a = dv[jj];
      if ((P / u) % a || gcd(u, a) != 1) continue;
      u64 al = P / u / a;
      int reg;
      if (u * u < P && a * a < P && al * al < P) reg = 0;
      else if (u >= 2 && a >= 2 && in_level(u * a, p)) reg = 1;
      else if (in_level(u * al, p)) reg = 2;
      else if (in_level(a * al, p)) reg = 3;
      else reg = 4;
      if (reg == 0) {
        ncentre++;
        /* strip: the largest coordinate is at least (1/2 - 1/42) log P */
        u64 mx = u > a ? (u > al ? u : al) : (a > al ? a : al);
        if (log((double)mx) >= (0.5 - 1.0 / 42) * lP) { nstrip++; sstrip += 1.0 / p; }
      }
      item *t = &it[n++];
      t->u = u; t->a = a; t->al = al; t->reg = reg;
      t->c = (u128)u * pw(a, p - 2, p) % p;
      /* the three normal forms */
      CHECK((u128)t->c * a % p == u % p, "class at p=%ld", p);
    }
    qsort(it, n, sizeof(item), cmp);
    for (int i = 0; i < n;) {
      int k = i; while (k < n && it[k].c == it[i].c) k++;
      for (int x = i; x < k; x++) for (int y = i; y < k; y++) {
        if (x == y) continue;
        item *s = &it[x], *t = &it[y];
        if (s->reg != 0) continue;
        if (t->reg == 0 && x < y) CHECK(0, "two centre triples share a class at p=%ld", p);
        if (t->reg == 1) CHECK(0, "centre and A share a class at p=%ld", p);
        if (t->reg == 2) {
          CHECK(s->u > s->a, "centre (u<a) and B share a class at p=%ld", p);
          /* mirror: (a, u, alpha) in the centre with a < u against (u'', a'', alpha'')^iota in C */
          item ms = {0, s->a, s->u, s->al, 0}, mt = {0, t->a, t->u, t->al, 3};
          family(p, P, &ms, &mt, &nfam_mirror, 0);
          sfam += 1.0 / p;
        }
        if (t->reg == 3) {
          CHECK(s->u < s->a, "centre (u>a) and C share a class at p=%ld", p);
          family(p, P, s, t, &nfam, 1);
          sfam += 1.0 / p;
          if (!first_p) first_p = p;
          if ((u64)p == first_p) nfirst++;
        }
      }
      i = k;
    }
  }
  printf("primes p = 3 (mod 4) up to %ld: %ld; centre triples %ld; strip triples %ld, sum 1/p %.6f\n",
         Y, nprimes, ncentre, nstrip, sstrip);
  printf("coincidences centre (u<a) - C: %ld; centre (u>a) - B: %ld; first at p = %llu; sum nu'/p = %.3e\n",
         nfam, nfam_mirror, first_p, sfam);
  if (fails) { printf("%ld failures\n", fails); return 1; }
  printf("centre: classes distinct, no coincidence with A, none of u<a with B or of u>a with C, "
         "every other coincidence in the family; OK\n");
  return 0;
}
