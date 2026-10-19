// chain.cpp -- brute-force check of the divisor-chain lemma of "An explicit exponent for the
// exceptional set of the Erdos-Straus equation".
//
// Lemma. Let 1 <= A < B, G = min(B/A, A), let T be the product of the prime powers of n whose
// prime is < G, and let q be the largest prime factor of n in [G, A] (q = 1 if there is none).
// If q T > A, then n has a divisor d with A < d < B.
//
// For every n <= NMAX and every window (A, B) of a grid, the program checks the lemma, checks that
// a prime factor in (A, B) is such a divisor, and records how much of the true set
// {n : n has a divisor in (A, B)} the two sufficient conditions cover.  As a control it also tests
// the rule without the step bound ("the A-smooth part of n exceeds A"), which is false; the program
// must find counterexamples to it.
//
// Usage: chain [NMAX]   (default 2000000)
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char** argv) {
  const long NMAX = argc > 1 ? std::atol(argv[1]) : 2000000L;
  std::vector<int> spf(NMAX + 1, 0);
  for (long i = 2; i <= NMAX; ++i)
    if (!spf[i])
      for (long j = i; j <= NMAX; j += i)
        if (!spf[j]) spf[j] = (int)i;
  // windows (A, B): A from 3 to 3000, B/A from 1.1 to 40
  std::vector<std::pair<double, double>> win;
  for (double A : {3.0, 7.5, 20.0, 55.0, 150.0, 400.0, 1100.0, 3000.0})
    for (double r : {1.1, 1.3, 1.7, 2.5, 4.0, 7.0, 15.0, 40.0}) win.push_back({A, A * r});
  const size_t W = win.size();
  std::vector<long> truth(W, 0), byprime(W, 0), bychain(W, 0), either(W, 0), falsehits(W, 0);
  long fails = 0, checked = 0;
  std::vector<long> divs;
  std::vector<std::pair<long, int>> fac;
  for (long n = 1; n <= NMAX; ++n) {
    fac.clear();
    for (long m = n; m > 1;) {
      int p = spf[m], e = 0;
      while (m % p == 0) { m /= p; ++e; }
      fac.push_back({p, e});
    }
    divs.assign(1, 1);
    for (auto [p, e] : fac) {
      size_t cur = divs.size();
      long pw = 1;
      for (int k = 1; k <= e; ++k) {
        pw *= p;
        for (size_t i = 0; i < cur; ++i) divs.push_back(divs[i] * pw);
      }
    }
    std::sort(divs.begin(), divs.end());
    for (size_t w = 0; w < W; ++w) {
      const double A = win[w].first, B = win[w].second, G = std::min(B / A, A);
      auto it = std::upper_bound(divs.begin(), divs.end(), (long)std::floor(A));
      const bool has = it != divs.end() && (double)*it < B;
      double T = 1, smoothA = 1;
      long q = 1;
      bool prime_in = false;
      for (auto [p, e] : fac) {
        double pe = std::pow((double)p, e);
        if (p < G) T *= pe;
        else if (p <= A) q = std::max(q, (long)p);
        if (p > A && p < B) prime_in = true;
        if (p <= A) smoothA *= pe;
      }
      const bool chain = (double)q * T > A;
      ++checked;
      if ((chain || prime_in) && !has) {
        if (++fails < 10) std::printf("FAIL: n=%ld A=%g B=%g\n", n, A, B);
      }
      truth[w] += has;
      byprime[w] += prime_in;
      bychain[w] += chain;
      either[w] += (chain || prime_in);
      if (smoothA > A && !prime_in && !has) ++falsehits[w];
    }
  }
  long ctrl = 0;
  for (size_t w = 0; w < W; ++w) ctrl += falsehits[w];
  std::printf("n <= %ld, %zu windows, %ld cases\n", NMAX, W, checked);
  std::printf("%8s %8s %10s %10s %10s\n", "A", "B/A", "true", "rule", "rule/true");
  for (size_t w = 0; w < W; w += 9)
    std::printf("%8.1f %8.1f %10ld %10ld %10.4f\n", win[w].first, win[w].second / win[w].first, truth[w],
                either[w], truth[w] ? (double)either[w] / truth[w] : 1.0);
  std::printf("control: the rule without the step bound fails %ld times\n", ctrl);
  if (fails) { std::printf("%ld failures of the lemma\n", fails); return 1; }
  if (ctrl == 0) { std::printf("the control rule did not fail; the test is not sensitive\n"); return 1; }
  std::printf("chain lemma: no counterexample; OK\n");
  return 0;
}
