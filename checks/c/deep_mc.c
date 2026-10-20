/* Robust rule with exact n-intervals: for each assignment of the parts >= G
   (bin width h) to three groups, each choice of leftover group L and targets
   rA, rB for the remainder (parts < G) gives the open interval of n on which
   every corner of the achievable (d2,d3) set meets (C1)-(C4); the modulus
   passes if these intervals cover T = [1/6 - th, 1/3]. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static uint64_t S = 88172645463325252ULL;
static double rnd(void){ S ^= S << 13; S ^= S >> 7; S ^= S << 17; return (S >> 11) * (1.0/9007199254740992.0); }
static double TAU = 0.003;
static int MODE = 0; /* 0 fill with reserve G between two groups (the paper), 1 one filled group, 2 no reserve, 3 reserve G per group */
static int pd(double s, double *p, double *dust){
  int m = 0; double rem = 1.0;
  while (rem * s > TAU){ double v = rem * rnd(); if (v * s > TAU) p[m++] = v * s; rem -= v; }
  *dust = rem * s; return m;
}
static int cmp(const void*a,const void*b){ double x=*(double*)a,y=*(double*)b; return (x<y)-(x>y);}
static int NS; static double sS, G, H, TLO, THI;
#define MAXI 4096
static double IL[MAXI], IR[MAXI]; static int NI; static int covered;
static int dcmp(const void*a,const void*b){ int i=*(int*)a,j=*(int*)b; return (IL[i]>IL[j])-(IL[i]<IL[j]); }
static int check_cover(void){
  static int idx[MAXI]; for (int i = 0; i < NI; i++) idx[i] = i;
  qsort(idx, NI, sizeof(int), dcmp);
  double reach = TLO; int started = 0;
  for (int k = 0; k < NI; k++){ int i = idx[k];
    if (IL[i] >= reach){ if (!(started == 0 && IL[i] < TLO)) { if (IL[i] >= reach) return 0; } }
    if (IR[i] > reach) reach = IR[i]; started = 1;
    if (reach > THI) return 1; }
  return 0;
}
static void add(double a, double b){
  if (a < TLO - 1) a = TLO - 1; if (b <= a) return;
  /* merge if contained in an existing interval */
  for (int i = 0; i < NI; i++) if (IL[i] <= a && IR[i] >= b) return;
  if (NI < MAXI){ IL[NI] = a; IR[NI] = b; NI++; }
}
static void leaf(const double *lo, const double *hi, double Rlo){
  for (int L = 0; L < 3; L++){
    int A = (L + 1) % 3, B = (L + 2) % 3;
    for (int i = 0; i <= NS; i++) for (int j = 0; j <= NS; j++){
      double rA = Rlo * i / NS, rB = Rlo * j / NS;
      if (MODE == 1 && i && j) continue;
      double used;
      if (MODE == 3) { used = (i ? rA + G : 0) + (j ? rB + G : 0); if (used > Rlo + 1e-12) continue; }
      else { used = rA + rB + ((MODE != 2 && i && j) ? G : 0); if (used >= Rlo) continue; }
      double Alo = lo[A] + rA, Ahi = hi[A] + (i ? rA + G : 0);
      double Blo = lo[B] + rB, Bhi = hi[B] + (j ? rB + G : 0);
      double nmin = -1, nmax = 2;
      for (int c = 0; c < 4; c++){
        double dA = (c & 1) ? Ahi : Alo, dB = (c & 2) ? Bhi : Blo;
        double d[3]; d[A] = dA; d[B] = dB; d[L] = sS - dA - dB;
        double d2 = d[1], d3 = d[2], t;
        t = sS - d2 - d3; if (t > nmin) nmin = t;
        t = (1 - d2 - 2*d3) / 2; if (t < nmax) nmax = t;
        t = (2 - sS - 3*d2 - 2*d3) / 2; if (t < nmax) nmax = t;
        t = 2 - sS - 4*d2 - d3; if (t < nmax) nmax = t;
      }
      if (nmax > nmin && nmax > TLO && nmin < THI) add(nmin, nmax);
    }
  }
}
static double z[64]; static int M; static double RLO;
static double cap0, cap1, cap2;
static void dfs(int i, double b0, double b1, double b2, int c0, int c1, int c2){
  if (covered) return;
  if (b0 - c0*H/2 > cap0 || b1 - c1*H/2 > cap1 || b2 - c2*H/2 > cap2) return;
  if (i == M){
    double lo[3] = {b0 - c0*H/2, b1 - c1*H/2, b2 - c2*H/2}, hi[3] = {b0 + c0*H/2, b1 + c1*H/2, b2 + c2*H/2};
    int before = NI; leaf(lo, hi, RLO);
    if (NI > before && check_cover()) covered = 1;
    return;
  }
  dfs(i+1, b0 + z[i], b1, b2, c0+1, c1, c2); dfs(i+1, b0, b1 + z[i], b2, c0, c1+1, c2); dfs(i+1, b0, b1, b2 + z[i], c0, c1, c2+1);
}
int main(int argc, char **argv){
  double th = atof(argv[1]); long N = atol(argv[2]); G = atof(argv[3]); H = atof(argv[4]);
  int MAXP = argc > 5 ? atoi(argv[5]) : 12;
  if (argc > 6) S += atol(argv[6]) * 0x9E3779B97F4A7C15ULL;
  NS = argc > 7 ? atoi(argv[7]) : 10;
  MODE = argc > 8 ? atoi(argv[8]) : 0;
  sS = 0.5 + th; TLO = 1.0/6 - th; THI = 1.0/3;
  /* each factor must stay below these at the most lenient n = TLO */
  cap0 = THI; cap2 = (1 - 2*TLO) / 2; cap1 = (2 - TLO - sS) / 4;
  long acc = 0, over = 0;
  for (long it = 0; it < N; it++){
    double s1 = 2*th + (0.5 - 3*th) * rnd(), s2 = sS - s1, p1[64], p2[64], D1, D2;
    int m1 = pd(s1, p1, &D1), m2 = pd(s2, p2, &D2);
    double all[128]; int ma = 0;
    for (int i = 0; i < m1; i++) all[ma++] = p1[i]; for (int i = 0; i < m2; i++) all[ma++] = p2[i];
    qsort(all, ma, sizeof(double), cmp);
    M = 0; double big = 0; while (M < ma && all[M] >= G){ z[M] = all[M]; big += all[M]; M++; }
    if (M > MAXP){ over++; continue; }
    RLO = sS - big - M*H/2;
    NI = 0; covered = 0; dfs(0, 0, 0, 0, 0, 0, 0); acc += covered;
  }
  printf("mode=%d th=%.4f G=%.3f h=%.4f NS=%d rule=%.4f over=%.4f (N=%ld)\n", MODE, th, G, H, NS, (double)acc/N, (double)over/N, N);
  return 0;
}
