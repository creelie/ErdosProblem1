/* constant_M.c -- sums T and M of the paper over the odd primes up to P in long double,
   subtracts the prime-number-theorem estimate (log P + 1)/P for the primes above P
   from M, and prints both values. Usage: constant_M P (default 10^8). */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  long P = argc > 1 ? atol(argv[1]) : 100000000L;
  char *c = calloc(P + 1, 1);
  if (!c) return 2;
  long double T = logl(2.0L) / 3, M = -4.0L / 9 * logl(2.0L) * logl(2.0L);
  for (long p = 2; p <= P; p++) {
    if (c[p]) continue;
    if (p <= P / p)
      for (long q = p * p; q <= P; q += p) c[q] = 1;
    if (p == 2) continue;
    long double x = p, l = logl(x);
    T += l / ((x - 1) * (x * x + 1));
    M -= l * l * ((x + 1) / ((x - 1) * (x * x + 1)) + 1 / ((x - 1) * (x - 1) * (x * x + 1) * (x * x + 1)));
  }
  long double tail = (logl((long double)P) + 1) / P;
  printf("P = %ld  T = %.12Lf  M (primes up to P) = %.12Lf  M (with tail) = %.12Lf\n", P, T, M, M - tail);
  int ok = fabsl(M - tail + 1.0429405079L) < 1e-9L && fabsl(T - 0.3131549610L) < 1e-10L;
  puts(ok ? "constant check passed" : "constant check FAILED");
  return !ok;
}
