// check.cpp -- exact rational check of the inversion h0 = k0 * (1 (x) 1):
// for all u, a <= 400, sum_{d|u, e|a} k0(d,e) = h0(u,a), where
//   h0(u,a) = [gcd(u,a)=1] * 2ua/phi(4ua)
//   k0(d,e) = 0 unless d = d1 m, e = e1 m with d1, e1, m squarefree, pairwise coprime,
//             d1 e1 odd; then k0 = prod_{p | d1 e1} 1/(p-1) * prod_{p | m, p odd} (-(p+1)/(p-1)) * (-1)^[2 | m].
// It also checks phi(4ua) = c phi(u) phi(a), c = 4 or 2, for coprime u, a <= 400.
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <vector>

using i128 = __int128;
struct Q {
  i128 n, d;
  Q(i128 a = 0, i128 b = 1) : n(a), d(b) { norm(); }
  static i128 g(i128 a, i128 b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { i128 t = a % b; a = b; b = t; } return a; }
  void norm() { if (d < 0) { n = -n; d = -d; } i128 h = g(n, d); if (h > 1) { n /= h; d /= h; } }
  Q operator+(const Q &o) const { return Q(n * o.d + o.n * d, d * o.d); }
  Q operator*(const Q &o) const { return Q(n * o.n, d * o.d); }
  bool operator==(const Q &o) const { return n == o.n && d == o.d; }
};

static const int N = 400;
static long phi_(long n) { long r = n; for (long p = 2; p * p <= n; p++) if (n % p == 0) { while (n % p == 0) n /= p; r -= r / p; } if (n > 1) r -= r / n; return r; }

static Q k0(long d, long e) {
  long m = std::gcd(d, e), d1 = d / m, e1 = e / m;
  if (std::gcd(d1, m) != 1 || std::gcd(e1, m) != 1) return Q(0);
  Q r(1);
  auto sqf = [](long x) { for (long p = 2; p * p <= x; p++) if (x % (p * p) == 0) return false; return true; };
  if (!sqf(d1) || !sqf(e1) || !sqf(m) || d1 % 2 == 0 || e1 % 2 == 0) return Q(0);
  for (long x : {d1 * e1}) { long y = x; for (long p = 3; p <= y; p += 2) if (y % p == 0) { r = r * Q(1, p - 1); y /= p; } }
  long y = m;
  if (y % 2 == 0) { r = r * Q(-1); y /= 2; }
  for (long p = 3; p <= y; p += 2) if (y % p == 0) { r = r * Q(-(p + 1), p - 1); y /= p; }
  return r;
}

int main() {
  int fail = 0;
  std::vector<std::vector<long>> divs(N + 1);
  for (long d = 1; d <= N; d++) for (long k = d; k <= N; k += d) divs[k].push_back(d);
  for (long u = 1; u <= N; u++)
    for (long a = 1; a <= N; a++) {
      Q s(0);
      for (long d : divs[u]) for (long e : divs[a]) s = s + k0(d, e);
      Q h = std::gcd(u, a) == 1 ? Q(2 * u * a, phi_(4 * u * a)) : Q(0);
      if (!(s == h)) { std::printf("mismatch at u=%ld a=%ld\n", u, a); fail = 1; }
      if (std::gcd(u, a) == 1) {
        long c = (u * a) % 2 == 0 ? 4 : 2;
        if (phi_(4 * u * a) != c * phi_(u) * phi_(a)) { std::printf("totient mismatch %ld %ld\n", u, a); fail = 1; }
      }
    }
  std::puts(fail ? "C++ check FAILED" : "C++ check passed: inversion and totient factorisation for u, a <= 400");
  return fail;
}
