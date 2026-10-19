/* check.c -- computes W(J) = sum over coprime u, a <= J of 1/phi(4ua) exactly in
   double precision and compares it with the three-term formula
   W(J) = K((log J + gamma + T)^2 + M) + error,  K = 9/(2 pi^2).
   T and M are summed over the primes below 10^7, with the tail of M estimated by the
   prime number theorem. Exit code 0 when every error is
   below 3 J^{-1/2} log J and the values at J = 100, 300, 1000, 3000 agree with the
   ratios 0.62795, 0.59503, 0.57114, 0.55559 reported by Dahan. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static long gcdl(long a, long b) { while (b) { long t = a % b; a = b; b = t; } return a; }

int main(void) {
  const int J = 20000, P = 10000000;
  long *ph = malloc(sizeof(long) * (J + 1));
  for (int i = 0; i <= J; i++) ph[i] = i;
  for (int p = 2; p <= J; p++) if (ph[p] == p) for (int k = p; k <= J; k += p) ph[k] -= ph[k] / p;
  char *comp = calloc(P + 1, 1);
  double T = log(2.0) / 3, M = -4.0 / 9 * log(2.0) * log(2.0);
  for (long p = 2; p <= P; p++) {
    if (comp[p]) continue;
    for (long q = p * p; q <= P; q += p) comp[q] = 1;
    if (p == 2) continue;
    double x = (double)p, l = log(x);
    T += l / ((x - 1) * (x * x + 1));
    M -= l * l * ((x + 1) / ((x - 1) * (x * x + 1)) + 1 / ((x - 1) * (x - 1) * (x * x + 1) * (x * x + 1)));
  }
  /* The terms of M decay like (log p)^2/p^2; the primes above P contribute about
     (log P + 1)/P by the prime number theorem, which is subtracted here. */
  M -= (log((double)P) + 1) / P;
  const double PI = 3.14159265358979323846, K = 9 / (2 * PI * PI), g = 0.57721566490153286;
  printf("K = %.12f  T = %.12f  M = %.12f\n", K, T, M);
  int fail0 = fabs(M + 1.0429405079) > 1e-8;
  int marks[] = {100, 300, 1000, 3000, 10000, 20000};
  double dahan[] = {0.62795, 0.59503, 0.57114, 0.55559};
  int mi = 0, fail = fail0;
  double W = 0;
  for (int m = 1; m <= J; m++) {
    for (int x = 1; x <= m; x++) {
      if (gcdl(m, x) != 1) continue;
      long u = m, a = x;
      double t = ((u * a) % 2 == 0) ? 1.0 / (4.0 * ph[u] * ph[a]) : 1.0 / (2.0 * ph[u] * ph[a]);
      W += (x == m) ? t : 2 * t;
    }
    if (m == marks[mi]) {
      double L = log((double)m), f = K * ((L + g + T) * (L + g + T) + M), e = W - f;
      printf("J = %5d  W = %.9f  formula = %.9f  error = %+.3e  W/log^2 J = %.5f\n", m, W, f, e, W / (L * L));
      if (fabs(e) > 3 * L / sqrt((double)m)) fail = 1;
      if (mi < 4 && fabs(W / (L * L) - dahan[mi]) > 6e-6) fail = 1;
      mi++;
    }
  }
  puts(fail ? "C check FAILED" : "C check passed");
  return fail;
}
