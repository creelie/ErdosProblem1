// mean_values.cpp -- the mean values of the corners in "An explicit exponent for the exceptional set
// of the Erdos-Straus equation":
//   V_A(X) = sum_{ua <= X, gcd(u,a) = 1, u, a >= 2} 1/phi(4ua),
//   V_B(X) = sum_{u alpha <= X} (phi(u)/u) / phi(4 u alpha),
// both equal to (kappa/2) log^2 X + O(log X) with kappa = 9/(2 pi^2).  The program prints
// (V - (kappa/2) log^2 X) / log X along X = 10^3, ..., 2 * 10^6 and checks that it stays below 1 in
// absolute value and settles: it moves by less than 0.01 between X = 10^5 and X = 2 * 10^6, so that
// V / ((kappa/2) log^2 X) tends to 1.
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

int main() {
  const long XM = 2000000;
  std::vector<long> ph(4 * XM + 1);
  std::iota(ph.begin(), ph.end(), 0L);
  for (long i = 2; i <= 4 * XM; i++)
    if (ph[i] == i)
      for (long j = i; j <= 4 * XM; j += i) ph[j] -= ph[j] / i;
  std::vector<double> fA(XM + 1, 0.0), fB(XM + 1, 0.0);
  for (long u = 1; u <= XM; u++)
    for (long m = u; m <= XM; m += u) {
      long a = m / u;
      double w = 1.0 / ph[4 * m];
      if (u >= 2 && a >= 2 && std::gcd(u, a) == 1) fA[m] += w;
      fB[m] += w * double(ph[u]) / double(u);
    }
  const double kappa = 9.0 / (2.0 * M_PI * M_PI);
  const long marks[] = {1000, 10000, 100000, 1000000, 2000000};
  double cA = 0, cB = 0, dA[5], dB[5];
  int mi = 0, bad = 0;
  for (long m = 1; m <= XM; m++) {
    cA += fA[m];
    cB += fB[m];
    if (mi < 5 && m == marks[mi]) {
      double l = std::log(double(m)), main = kappa / 2 * l * l;
      std::printf("X = %7ld  V_A = %9.5f  V_B = %9.5f  (kappa/2) log^2 X = %9.5f  "
                  "(V_A - main)/log X = %.4f  (V_B - main)/log X = %.4f\n",
                  m, cA, cB, main, (cA - main) / l, (cB - main) / l);
      dA[mi] = (cA - main) / l;
      dB[mi] = (cB - main) / l;
      if (std::fabs(dA[mi]) > 1 || std::fabs(dB[mi]) > 1) bad++;
      mi++;
    }
  }
  if (std::fabs(dA[4] - dA[2]) > 0.01 || std::fabs(dB[4] - dB[2]) > 0.01) bad++;
  if (bad) {
    std::printf("mean_values.cpp: %d failures\n", bad);
    return 1;
  }
  std::printf("mean_values.cpp: all checks passed\n");
  return 0;
}
