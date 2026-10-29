/* certify_K.c -- lower bound for the factorable share of the centre on a cell of depths.
 *
 * Configurations: the prime factors of q = 4ua of size >= G (exponents of x), the largest of
 * them z1, the others in a multiset, and R, the total of the prime factors below G.  The
 * largest factor is solved for, z1 = s - S' - R, so a cell is a multiset of intervals for the
 * other factors and an interval for R.  The measure of a cell is
 *     prod_J (int_J dz/z)^{m_J} / m_J!  x  sum over labels  int Psi_c(R) / z1 dR,
 * where Psi_c(R) = int rho(t/G) rho((R - t)/G) dt over the smooth parts of u and a with the
 * strip constraints of the label; the R-integral uses exact band tables and Jensen's
 * inequality for 1/z1.  All masses use directed rounding.  The pass check uses fixed targets
 * for the greedy fill of two groups by the factors below G, and every covering interval is
 * certified by two witness points checked with directed rounding.
 *
 * Usage: certify_K tha thb G H0 RMINB KMAX [EPS] [ALLPASS] [TAU] [HT]
 */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#pragma STDC FENV_ACCESS ON

/* directed operations; valid while the rounding mode is FE_UPWARD */
static double dn_add(double a, double b){ volatile double x = -a; volatile double r = x - b; return -r; }
static double dn_sub(double a, double b){ volatile double r = b - a; return -r; }
static double dn_mul(double a, double b){ volatile double x = -a; volatile double r = x * b; return -r; }
static double dn_div(double a, double b){ volatile double x = -a; volatile double r = x / b; return -r; }
static double up_add(double a, double b){ volatile double r = a + b; return r; }
static double up_sub(double a, double b){ volatile double r = a - b; return r; }
static double up_mul(double a, double b){ volatile double r = a * b; return r; }
static double up_div(double a, double b){ volatile double r = a / b; return r; }
#define UP() fesetround(FE_UPWARD)
#define NEAR() fesetround(FE_TONEAREST)

/* ---------------- Dickman's function on the grid i / NR ---------------- */
#define NR 2048
#define UMAX 20
#define NRHO (UMAX * NR + 2)
static double rho_lo[NRHO], rho_up[NRHO];
static void dickman(void){
  static double cuu[NRHO + 1], cdu[NRHO + 1], cul[NRHO + 1], cdl[NRHO + 1];
  UP();
  cuu[0] = cdu[0] = cul[0] = cdl[0] = 0;
  for (int i = 0; i < NRHO; i++){
    if (i <= NR){ rho_lo[i] = rho_up[i] = 1; }
    else {
      /* i rho(i/NR) / NR = int_{u-1}^{u} rho; left sums from above, right sums from below */
      double su = up_sub(cuu[i], cdu[i - NR]);              /* >= sum_{j=i-NR}^{i-1} rho_up[j] */
      rho_up[i] = up_div(su, (double)i);
      double sl = dn_sub(cdl[i], cul[i - NR + 1]);          /* <= sum_{j=i-NR+1}^{i-1} rho_lo[j] */
      if (sl < 0) sl = 0;
      rho_lo[i] = dn_div(sl, (double)(i - 1));
      if (rho_up[i] > rho_up[i - 1]) rho_up[i] = rho_up[i - 1];
      if (rho_lo[i] > rho_up[i]) { fprintf(stderr, "dickman bounds cross at %d\n", i); exit(2); }
    }
    cuu[i + 1] = up_add(cuu[i], rho_up[i]); cdu[i + 1] = dn_add(cdu[i], rho_up[i]);
    cul[i + 1] = up_add(cul[i], rho_lo[i]); cdl[i + 1] = dn_add(cdl[i], rho_lo[i]);
  }
  NEAR();
}

/* ---------------- band tables for the smooth parts ---------------- */
/* sub-grid step dp = G / NR (so rho(t/G) at t = a dp is rho(a/NR)); bands of BW sub-steps,
   band i is R in [i dR, (i+1) dR], dR = BW dp.  NB bands cover R < NB dR = 16 G. */
#define BW 16
#define NB (16 * NR / BW)
static double G, dp, dR;
static int CMAX;                       /* constraint index c (t >= c dR) up to CMAX */
static double *CB1lo, *CB1up;         /* [(c)*(NB+1) + i]: cumulative over bands < i, constraint c */
static double CBUlo[NB + 1], CBUup[NB + 1], CM1lo[NB + 1], CM1up[NB + 1];
#define IDX(c, i) ((size_t)(c) * (NB + 1) + (i))
static void bands(void){
  int NA = NB * BW + 2 * BW;
  double *f = malloc(sizeof(double) * (NA + 2)), *F = malloc(sizeof(double) * (NA + 2));
  double *PL = malloc(sizeof(double) * (NA + 3)), *PU = malloc(sizeof(double) * (NA + 3));
  double *QL = malloc(sizeof(double) * (NA + 3)), *QU = malloc(sizeof(double) * (NA + 3));
  double *g = malloc(sizeof(double) * (NA + 2));
  CB1lo = calloc((size_t)(CMAX + 1) * (NB + 1), sizeof(double));
  CB1up = calloc((size_t)(CMAX + 1) * (NB + 1), sizeof(double));
  UP();
  for (int a = 0; a <= NA; a++){ f[a] = rho_lo[a + 1]; F[a] = rho_up[a]; }
  PL[0] = PU[0] = QL[0] = QU[0] = 0;
  for (int a = 0; a <= NA; a++){ PL[a + 1] = dn_add(PL[a], f[a]); PU[a + 1] = up_add(PU[a], f[a]);
                                 QL[a + 1] = dn_add(QL[a], F[a]); QU[a + 1] = up_add(QU[a], F[a]); }
  double area = dn_mul(dp, dp), half = dn_mul(0.5, area);
  double areau = up_mul(dp, dp), halfu = up_mul(0.5, areau);
  for (int c = 0; c <= CMAX; c++){ CB1lo[IDX(c, 0)] = 0; CB1up[IDX(c, 0)] = 0; }
  CBUlo[0] = CBUup[0] = CM1lo[0] = CM1up[0] = 0;
  for (int i = 0; i < NB; i++){
    int lo = BW * i, hi = BW * i + BW - 2;          /* full squares: lo <= a + b <= hi */
    int amax = BW * i + BW - 1;
    double bu = 0;                                   /* upper bound of the band (max values) */
    for (int a = 0; a <= amax; a++){
      int b0 = lo - a; if (b0 < 0) b0 = 0; int b1 = hi - a;
      double s = 0, su = 0;
      if (b1 >= b0){ s = dn_sub(PL[b1 + 1], PU[b0]); if (s < 0) s = 0; su = up_sub(QU[b1 + 1], QL[b0]); }
      double v = dn_mul(area, dn_mul(f[a], s));
      double vu = up_mul(areau, up_mul(F[a], su));
      /* half squares a + b = lo - 1 (upper triangle) and a + b = hi + 1 (lower triangle) */
      int bh = lo - 1 - a; if (bh >= 0){ v = dn_add(v, dn_mul(half, dn_mul(f[a], f[bh]))); vu = up_add(vu, up_mul(halfu, up_mul(F[a], F[bh]))); }
      bh = hi + 1 - a; if (bh >= 0){ v = dn_add(v, dn_mul(half, dn_mul(f[a], f[bh]))); vu = up_add(vu, up_mul(halfu, up_mul(F[a], F[bh]))); }
      g[a] = v; bu = up_add(bu, vu);
    }
    /* suffix sums over a, sampled at a = BW c */
    double suf = 0;
    for (int a = amax; a >= 0; a--){
      suf = dn_add(suf, g[a]);
      if (a % BW == 0){ int c = a / BW; if (c <= CMAX){
        CB1lo[IDX(c, i + 1)] = dn_add(CB1lo[IDX(c, i)], suf);
        CB1up[IDX(c, i + 1)] = up_add(CB1up[IDX(c, i)], suf); } }
    }
    for (int c = amax / BW + 1; c <= CMAX; c++){ CB1lo[IDX(c, i + 1)] = CB1lo[IDX(c, i)]; CB1up[IDX(c, i + 1)] = CB1up[IDX(c, i)]; }
    double Ri = dn_mul((double)i, dR);
    double b0 = suf;  /* band lower bound, unconstrained */
    CBUlo[i + 1] = dn_add(CBUlo[i], bu); CBUup[i + 1] = up_add(CBUup[i], bu);
    CM1lo[i + 1] = dn_add(CM1lo[i], dn_mul(Ri, b0)); CM1up[i + 1] = up_add(CM1up[i], dn_mul(Ri, b0));
  }
  NEAR();
  free(f); free(F); free(PL); free(PU); free(QL); free(QU); free(g);
}
/* lower bound of int_{i0 dR}^{i1 dR} Psi_c(R) dR, constraint t >= c dR on one side */
static double Qlo(int c, int i0, int i1){ double q = dn_sub(CB1lo[IDX(c, i1)], CB1up[IDX(c, i0)]); return q > 0 ? q : 0; }
static double Qup0(int i0, int i1){ return up_sub(CBUup[i1], CBUlo[i0]); }

/* ---------------- the depth cell ---------------- */
static double tha, thb, sa, sb, TLO, THI = 1.0 / 3, cap[3], H0, EPS, TAU, HT; static int RMINB, KMAX, ALLPASS;

/* ---------------- the pass check ---------------- */
#define MARG 1e-9
typedef struct { double a, b, c; } line_t;
static const double CAN[4] = {1, -2, -2, -1};
static const double CC[4][3] = {{-1,0,0},{0,-1,-2},{0,-3,-2},{0,-4,-1}};
static const double CE[4] = {0, 1, 2, 2};
static const double CS[4] = {0, 0, -1, -1};
/* group data of an assignment: others in each group, the group p of z1, R range */
typedef struct { double Olo[3], Ohi[3], OloT, OhiT; int p; double Rlo, Rhi; } gd_t;
typedef struct { double blo[3], bhi[3]; } bx_t;
/* boxes of the group sums (round to nearest; used for the search) */
static void boxes(const gd_t *d, bx_t *b){
  for (int j = 0; j < 3; j++){
    if (j == d->p){ b->blo[j] = sa - d->Rhi - (d->OhiT - d->Ohi[j]); b->bhi[j] = sb - d->Rlo - (d->OloT - d->Olo[j]); }
    else { b->blo[j] = d->Olo[j]; b->bhi[j] = d->Ohi[j]; }
  }
}
/* min over the cell of b_A + R (single fill of A) and of s - b_L (double fill) */
static double budA_of(const gd_t *d, int A){ return d->p == A ? sa - (d->OhiT - d->Ohi[A]) : d->Olo[A] + d->Rlo; }
static double budAB_of(const gd_t *d, int L){ return d->p == L ? (d->OloT - d->Olo[L]) + d->Rlo : sa - d->Ohi[L]; }
/* The configuration-dependent part of constraint k in mode (L, fA, fB) is
     kappa_A D_A + kappa_B D_B + sigma s,
   with D_X in [x, x + G) for a filled group and D_X = O_X + [X = p] z1, z1 = s - R - O_1 - O_2 - O_3,
   for an unfilled one.  In the independent coordinates O_1, O_2, O_3, R, s of the cell it is affine,
   and its minimum over the cell is taken at the corner chosen by the signs of the coefficients. */
typedef struct { double cO[3], cR, cs, K, cx, cy; } aff_t;
static void affine(int p, int L, int fA, int fB, int k, aff_t *a){
  int A = (L + 1) % 3, B = (L + 2) % 3;
  double kA = CC[k][A] - CC[k][L], kB = CC[k][B] - CC[k][L], sg = CC[k][L] + CS[k];
  a->cO[0] = a->cO[1] = a->cO[2] = 0; a->cR = 0; a->cs = sg; a->K = CE[k]; a->cx = a->cy = 0;
  int X[2] = {A, B}, f[2] = {fA, fB}; double kk[2] = {kA, kB};
  for (int i = 0; i < 2; i++){
    if (f[i]){ if (i == 0) a->cx = kk[i]; else a->cy = kk[i]; continue; }
    a->cO[X[i]] += kk[i];
    if (X[i] == p){ a->cs += kk[i]; a->cR -= kk[i]; for (int j = 0; j < 3; j++) a->cO[j] -= kk[i]; }
  }
}
static double Kmin(const gd_t *d, int L, int fA, int fB, int k){
  aff_t a; affine(d->p, L, fA, fB, k, &a);
  double K = a.K;
  if (fA && a.cx < 0) K += a.cx * G;
  if (fB && a.cy < 0) K += a.cy * G;
  for (int j = 0; j < 3; j++) K += a.cO[j] * (a.cO[j] >= 0 ? d->Olo[j] : d->Ohi[j]);
  K += a.cR * (a.cR >= 0 ? d->Rlo : d->Rhi) + a.cs * (a.cs >= 0 ? sa : sb);
  return K;
}
/* n-interval of a mode with witness points (x, y) at both ends */
static int mode3(const gd_t *d, const bx_t *bx, int L, int fA, int fB, double *nlo, double *nhi, double *w){
  int A = (L + 1) % 3, B = (L + 2) % 3;
  line_t reg[8]; int nr = 0;
  if (fA) reg[nr++] = (line_t){-1, 0, -(bx->bhi[A] + MARG)}; else { reg[nr++] = (line_t){1, 0, 0}; reg[nr++] = (line_t){-1, 0, 0}; }
  if (fB) reg[nr++] = (line_t){0, -1, -(bx->bhi[B] + MARG)}; else { reg[nr++] = (line_t){0, 1, 0}; reg[nr++] = (line_t){0, -1, 0}; }
  if (fA && fB) reg[nr++] = (line_t){1, 1, budAB_of(d, L) - G - MARG};
  else if (fA) reg[nr++] = (line_t){1, 0, budA_of(d, A) - MARG};
  else if (fB) reg[nr++] = (line_t){0, 1, budA_of(d, B) - MARG};
  double l0 = 0, lA = 0, lB = 0, u0[3], uA[3], uB[3]; int nu = 0;
  for (int k = 0; k < 4; k++){
    double kA = CC[k][A] - CC[k][L], kB = CC[k][B] - CC[k][L];
    double K = Kmin(d, L, fA, fB, k);
    double cA = fA ? kA : 0, cB = fB ? kB : 0;
    if (CAN[k] > 0){ l0 = -K / CAN[k] + MARG; lA = -cA / CAN[k]; lB = -cB / CAN[k]; }
    else { double ww = -CAN[k]; u0[nu] = K / ww - MARG; uA[nu] = cA / ww; uB[nu] = cB / ww; nu++; }
  }
  for (int j = 0; j < nu; j++) reg[nr++] = (line_t){lA - uA[j], lB - uB[j], u0[j] - l0 - MARG};
  line_t ln[16]; int nl = 0;
  for (int i = 0; i < nr; i++) ln[nl++] = reg[i];
  for (int i = 0; i < nu; i++) for (int j = i + 1; j < nu; j++) ln[nl++] = (line_t){uA[i] - uA[j], uB[i] - uB[j], u0[j] - u0[i]};
  double best_hi = -1e9, best_lo = 1e9; int any = 0;
  for (int i = 0; i < nl; i++) for (int j = i + 1; j < nl; j++){
    double det = ln[i].a * ln[j].b - ln[i].b * ln[j].a;
    if (fabs(det) < 1e-12) continue;
    double x = (ln[i].c * ln[j].b - ln[i].b * ln[j].c) / det, y = (ln[i].a * ln[j].c - ln[i].c * ln[j].a) / det;
    int ok = 1;
    for (int t = 0; t < nr && ok; t++) if (reg[t].a * x + reg[t].b * y > reg[t].c + 1e-12) ok = 0;
    if (!ok) continue;
    any = 1;
    double l = l0 + lA * x + lB * y; if (l < best_lo){ best_lo = l; w[0] = x; w[1] = y; }
    double m = 1e9; for (int t = 0; t < nu; t++){ double u = u0[t] + uA[t] * x + uB[t] * y; if (u < m) m = u; }
    if (m > best_hi){ best_hi = m; w[2] = x; w[3] = y; }
  }
  if (!any || best_hi < best_lo) return 0;
  *nlo = best_lo; *nhi = best_hi; return 1;
}
/* rigorous check (directed rounding) that (x, y, n) is feasible for the mode for every
   configuration of the cell; rb holds rigorous group data */
typedef struct { double blo[3], bhi[3], budA[3], budAB[3], Olo[3], Ohi[3], Rlo, Rhi; int p; } rb_t;
static int verify_pt(const rb_t *r, int L, int fA, int fB, double x, double y, double n){
  int A = (L + 1) % 3, B = (L + 2) % 3;
  if (fA && !(x > r->bhi[A])) return 0;
  if (fB && !(y > r->bhi[B])) return 0;
  if (fA && fB){ if (!(up_add(x, y) < dn_sub(r->budAB[L], G))) return 0; }
  else if (fA){ if (!(x < r->budA[A])) return 0; }
  else if (fB){ if (!(y < r->budA[B])) return 0; }
  double xu = up_add(x, G), yu = up_add(y, G);
  for (int k = 0; k < 4; k++){
    aff_t a; affine(r->p, L, fA, fB, k, &a);          /* small integer coefficients, exact */
    double t = dn_mul(CAN[k], n);
    if (fA) t = dn_add(t, dn_mul(a.cx, a.cx >= 0 ? x : xu));
    if (fB) t = dn_add(t, dn_mul(a.cy, a.cy >= 0 ? y : yu));
    for (int j = 0; j < 3; j++) if (a.cO[j] != 0) t = dn_add(t, dn_mul(a.cO[j], a.cO[j] >= 0 ? r->Olo[j] : r->Ohi[j]));
    if (a.cR != 0) t = dn_add(t, dn_mul(a.cR, a.cR >= 0 ? r->Rlo : r->Rhi));
    t = dn_add(t, dn_mul(a.cs, a.cs >= 0 ? sa : sb));
    t = dn_add(t, a.K);
    if (!(t > 0)) return 0;
  }
  return 1;
}

/* validation: random configurations of the cell, direct evaluation of C1-C4 */
static int VALID; static long nsamp, nviol;
static unsigned long long VS = 88172645463325252ULL;
static double vr(void){ VS ^= VS << 13; VS ^= VS >> 7; VS ^= VS << 17; return (VS >> 11) * (1.0 / 9007199254740992.0); }
static void sample_check(const rb_t *r, int L, int fA, int fB, double x, double y, double n){
  int A = (L + 1) % 3, B = (L + 2) % 3;
  for (int it = 0; it < 64; it++){
    double O[3], b[3], d[3];
    for (int j = 0; j < 3; j++){ double u = vr(); if (it < 8) u = (it >> (j % 3)) & 1; O[j] = r->Olo[j] + u * (r->Ohi[j] - r->Olo[j]); }
    double R = r->Rlo + vr() * (r->Rhi - r->Rlo), s_ = sa + vr() * (sb - sa);
    if (it < 8){ R = (it & 1) ? r->Rhi : r->Rlo; s_ = (it & 2) ? sb : sa; }
    double z1 = s_ - R - O[0] - O[1] - O[2];
    for (int j = 0; j < 3; j++) b[j] = O[j] + (j == r->p ? z1 : 0);
    /* the fill must be possible for this configuration */
    if (fA && x < b[A]) { nviol++; continue; }
    if (fB && y < b[B]) { nviol++; continue; }
    if (fA && fB && x + y + G > s_ - b[L]) { nviol++; continue; }
    if (fA && !fB && x > b[A] + R) { nviol++; continue; }
    if (fB && !fA && y > b[B] + R) { nviol++; continue; }
    d[A] = fA ? x + G * vr() * 0.999999 : b[A];
    d[B] = fB ? y + G * vr() * 0.999999 : b[B];
    d[L] = s_ - d[A] - d[B];
    for (int k = 0; k < 4; k++){
      double t = CAN[k] * n + CC[k][0] * d[0] + CC[k][1] * d[1] + CC[k][2] * d[2] + CS[k] * s_ + CE[k];
      if (!(t > 0)){ nviol++; if (nviol < 10) fprintf(stderr, "violation k=%d t=%g p=%d L=%d f=%d%d n=%g\n", k, t, r->p, L, fA, fB, n); break; }
    }
    nsamp++;
  }
}
/* ---------------- cells ---------------- */
#define MAXJ 48
typedef struct { int nj; double lo[MAXJ], hi[MAXJ]; int m[MAXJ]; int r0, r1; int shrunk; } cell_t;
static long nodes, npass, nfail, npart, ndrop, nverr, nshrink;
static double Vpass, Vfail, Vdrop, Vpart, Vslack; static long ntight;

/* coverage */
#define MAXI 4096
static double IL[MAXI], IR[MAXI]; static int NI, covered;
static int icmp(const void *a, const void *b){ int i = *(const int *)a, j = *(const int *)b; return (IL[i] > IL[j]) - (IL[i] < IL[j]); }
static int check_cover(void){
  static int idx[MAXI]; for (int i = 0; i < NI; i++) idx[i] = i;
  qsort(idx, NI, sizeof(int), icmp);
  double reach = TLO;
  for (int k = 0; k < NI; k++){ int i = idx[k]; if (IL[i] > reach) return 0; if (IR[i] > reach) reach = IR[i]; if (reach >= THI) return 1; }
  return 0;
}
static int dominated(double a, double b){ if (b < TLO || a > THI) return 1; for (int i = 0; i < NI; i++) if (IL[i] <= a && IR[i] >= b) return 1; return 0; }

static const cell_t *CUR; static int occ[MAXJ], nocc; static int comp[MAXJ][3]; static double z1lo_c;
static void rig_boxes2(int p, rb_t *r){
  double Olo[3] = {0,0,0}, Ohi[3] = {0,0,0};
  for (int t = 0; t < nocc; t++){ int J = occ[t];
    for (int j = 0; j < 3; j++) if (comp[J][j]){ Olo[j] = dn_add(Olo[j], dn_mul(comp[J][j], CUR->lo[J])); Ohi[j] = up_add(Ohi[j], up_mul(comp[J][j], CUR->hi[J])); } }
  double Rlo = dn_mul((double)CUR->r0, dR), Rhi = up_mul((double)CUR->r1, dR);
  for (int j = 0; j < 3; j++){
    int j1 = (j + 1) % 3, j2 = (j + 2) % 3;
    double oth_hi = up_add(Ohi[j1], Ohi[j2]), oth_lo = dn_add(Olo[j1], Olo[j2]);
    if (j == p){ r->blo[j] = dn_sub(dn_sub(sa, Rhi), oth_hi); r->bhi[j] = up_sub(up_sub(sb, Rlo), oth_lo); }
    else { r->blo[j] = Olo[j]; r->bhi[j] = Ohi[j]; }
    r->budA[j] = (p == j) ? dn_sub(sa, oth_hi) : dn_add(Olo[j], Rlo);
    r->budAB[j] = (p == j) ? dn_add(oth_lo, Rlo) : dn_sub(sa, Ohi[j]);
    r->Olo[j] = Olo[j]; r->Ohi[j] = Ohi[j];
  }
  r->Rlo = Rlo; r->Rhi = Rhi; r->p = p;
}
static void leaf(int p, const double *Olo, const double *Ohi){
  gd_t d; bx_t bx;
  for (int j = 0; j < 3; j++){ d.Olo[j] = Olo[j]; d.Ohi[j] = Ohi[j]; }
  d.OloT = Olo[0] + Olo[1] + Olo[2]; d.OhiT = Ohi[0] + Ohi[1] + Ohi[2]; d.p = p;
  d.Rlo = CUR->r0 * dR; d.Rhi = CUR->r1 * dR;
  boxes(&d, &bx);
  for (int j = 0; j < 3; j++) if (bx.blo[j] > cap[j]) return;
  int before = NI, haverb = 0; rb_t rb;
  for (int L = 0; L < 3; L++) for (int f = 0; f < 4; f++){
    double a, b, w[4]; int fA = f & 1, fB = f >> 1;
    if (!mode3(&d, &bx, L, fA, fB, &a, &b, w)) continue;
    if (dominated(a, b)) continue;
    if (!haverb){ UP(); rig_boxes2(p, &rb); NEAR(); haverb = 1; }
    UP(); int ok = verify_pt(&rb, L, fA, fB, w[0], w[1], a) && verify_pt(&rb, L, fA, fB, w[2], w[3], b); NEAR();
    if (!ok){ nverr++; continue; }
    if (VALID) for (int e = 0; e < 2; e++) sample_check(&rb, L, fA, fB, w[2*e], w[2*e+1], e ? b : a);
    if (NI < MAXI){ IL[NI] = a; IR[NI] = b; NI++; }
  }
  if (NI > before && check_cover()) covered = 1;
}
static void adfs(int t, int p, double *Olo, double *Ohi){
  if (covered) return;
  /* prune by the caps: lower bounds of the group sums */
  for (int j = 0; j < 3; j++){ double lb = Olo[j] + (j == p ? z1lo_c : 0); if (lb > cap[j]) return; }
  if (t == nocc){ leaf(p, Olo, Ohi); return; }
  int J = occ[t], m = CUR->m[J]; double lo = CUR->lo[J], hi = CUR->hi[J];
  for (int a = 0; a <= m; a++) for (int b = 0; a + b <= m; b++){
    int c = m - a - b; comp[J][0] = a; comp[J][1] = b; comp[J][2] = c;
    Olo[0] += a * lo; Ohi[0] += a * hi; Olo[1] += b * lo; Ohi[1] += b * hi; Olo[2] += c * lo; Ohi[2] += c * hi;
    adfs(t + 1, p, Olo, Ohi);
    Olo[0] -= a * lo; Ohi[0] -= a * hi; Olo[1] -= b * lo; Ohi[1] -= b * hi; Olo[2] -= c * lo; Ohi[2] -= c * hi;
    if (covered) return;
  }
}
static int passes(const cell_t *c, double z1lo){
  if (ALLPASS) return 1;
  CUR = c; z1lo_c = z1lo; nocc = 0;
  for (int J = c->nj - 1; J >= 0; J--) if (c->m[J]) occ[nocc++] = J;   /* largest first */
  NI = 0; covered = 0;
  for (int p = 0; p < 3 && !covered; p++){ double Olo[3] = {0,0,0}, Ohi[3] = {0,0,0}; adfs(0, p, Olo, Ohi); }
  return covered;
}

/* ---------------- masses ---------------- */
static double fact[64];
/* lower bound of the measure of a valid cell (all directed); also an upper bound proxy */
static double mass_lo(const cell_t *c, double z1lo, double *ub){
  UP();
  double zm = 1, ES = 0, zmu = 1; int k = 1;
  double twoth = up_mul(2, thb);
  for (int J = 0; J < c->nj; J++) if (c->m[J]){
    /* ln(hi/lo) = log1p((hi - lo)/lo); glibc log1p is within 1 ulp, two ulps of margin */
    double dlo = dn_div(dn_sub(c->hi[J], c->lo[J]), c->lo[J]), dhi = up_div(up_sub(c->hi[J], c->lo[J]), c->lo[J]);
    NEAR(); double Llo = nextafter(nextafter(log1p(dlo), 0), 0), Lup = nextafter(nextafter(log1p(dhi), 1), 1); UP();
    double mu = dn_div(dn_sub(c->hi[J], c->lo[J]), Lup);
    double t = 1, tu = 1; for (int i = 0; i < c->m[J]; i++){ t = dn_mul(t, Llo); tu = up_mul(tu, Lup); }
    zm = dn_mul(zm, dn_div(t, fact[c->m[J]])); zmu = up_mul(zmu, up_div(tu, fact[c->m[J]]));
    ES = dn_add(ES, dn_mul(c->m[J], mu)); k += c->m[J];
  }
  int i0 = c->r0, i1 = c->r1;
  double Rlo = dn_mul((double)i0, dR);
  /* expected R under Psi on the box, from below */
  double num = dn_sub(CM1lo[i1], CM1up[i0]), den = Qup0(i0, i1);
  double ER = (den > 0 && num > 0) ? dn_div(num, den) : Rlo; if (ER < Rlo) ER = Rlo;
  double inv0 = dn_div(1, up_sub(up_sub(sb, ES), ER));
  double inv1 = dn_div(1, up_sub(up_sub(sb, ES), Rlo));
  /* labels: for every interval the number j of factors on the side of u, and the side of z1 */
  int occJ[MAXJ], no = 0; for (int J = 0; J < c->nj; J++) if (c->m[J]) occJ[no++] = J;
  int jj[MAXJ]; for (int t = 0; t < no; t++) jj[t] = 0;
  double W0 = 0, M1 = 0;
  for (;;){
    double w = 1, Su = 0, Sa = 0;
    for (int t = 0; t < no; t++){ int J = occJ[t]; int m = c->m[J], j = jj[t];
      w *= fact[m] / (fact[j] * fact[m - j]);                         /* exact integers */
      Su = dn_add(Su, dn_mul(j, c->lo[J])); Sa = dn_add(Sa, dn_mul(m - j, c->lo[J])); }
    for (int side = 0; side < 2; side++){
      double su = side ? dn_add(Su, z1lo) : Su, sa_ = side ? Sa : dn_add(Sa, z1lo);
      double cu = up_sub(twoth, su), ca = up_sub(twoth, sa_);
      if (cu <= 0 && ca <= 0) W0 += w;
      else if (cu <= 0 || ca <= 0){
        double cc = cu > 0 ? cu : ca; int ci = (int)ceil(cc / dR);
        if (dn_mul(ci, dR) < cc) ci++;
        if (ci <= CMAX && ci < i1) M1 = dn_add(M1, dn_mul(w, Qlo(ci, i0, i1)));
      }
    }
    int t = 0; while (t < no && jj[t] == c->m[occJ[t]]){ jj[t] = 0; t++; }
    if (t == no) break; jj[t]++;
  }
  double m0 = dn_mul(dn_mul(W0, Qlo(0, i0, i1)), inv0);
  double m1 = dn_mul(M1, inv1);
  double res = dn_mul(zm, dn_add(m0, m1));
  /* crude upper bound: all 2^k labels, Psi from above, 1/z1 <= 1/z1lo */
  *ub = up_mul(up_mul(zmu, up_mul(ldexp(1.0, k), Qup0(i0, i1))), up_div(1, z1lo));
  NEAR();
  return res;
}


/* estimate of the measure of the valid part of a cell (round to nearest; only for refinement) */
static double mass_est(const cell_t *c, double z1lo, double z1hi, double zmin){
  double zm = 1, ES = 0; int no = 0, occJ[MAXJ];
  for (int J = 0; J < c->nj; J++) if (c->m[J]){
    double L = log(c->hi[J] / c->lo[J]); zm *= pow(L, c->m[J]) / fact[c->m[J]];
    ES += c->m[J] * (c->hi[J] - c->lo[J]) / L; occJ[no++] = J; }
  int i0 = c->r0, i1 = c->r1;
  double Q0 = CB1lo[IDX(0, i1)] - CB1lo[IDX(0, i0)]; if (Q0 <= 0) return 0;
  double ER = (CM1lo[i1] - CM1lo[i0]) / Q0;
  double sm = 0.5 * (sa + sb), thm = 0.5 * (tha + thb), z1m = sm - ES - ER;
  double frac = 1;
  if (z1lo < zmin){ frac = (z1hi - zmin) / (z1hi - z1lo); if (frac <= 0) return 0; if (frac > 1) frac = 1; if (z1m < zmin) z1m = zmin; }
  int jj[MAXJ]; for (int t = 0; t < no; t++) jj[t] = 0;
  double W0 = 0, M1 = 0;
  for (;;){
    double w = 1, Su = 0, Sa = 0;
    for (int t = 0; t < no; t++){ int J = occJ[t]; int m = c->m[J], j = jj[t]; double mid = 0.5 * (c->lo[J] + c->hi[J]);
      w *= fact[m] / (fact[j] * fact[m - j]); Su += j * mid; Sa += (m - j) * mid; }
    for (int side = 0; side < 2; side++){
      double su = side ? Su + z1m : Su, sa_ = side ? Sa : Sa + z1m;
      double cu = 2 * thm - su, ca = 2 * thm - sa_;
      if (cu <= 0 && ca <= 0) W0 += w;
      else if (cu <= 0 || ca <= 0){ double cc = cu > 0 ? cu : ca; int ci = (int)floor(cc / dR + 0.5);
        if (ci <= CMAX && ci < i1){ double q = CB1lo[IDX(ci, i1)] - CB1lo[IDX(ci, i0)]; if (q > 0) M1 += w * q; } }
    }
    int t = 0; while (t < no && jj[t] == c->m[occJ[t]]){ jj[t] = 0; t++; }
    if (t == no) break;
    jj[t]++;
  }
  return frac * zm * (W0 * Q0 + M1) / z1m;
}

/* ---------------- the tree ---------------- */
/* split the cell along its largest uncertainty; returns 0 if no side is wider than the limits */
static void explore(cell_t *c, int known);
static int split(const cell_t *c, int known, int rmin, double hmin){
  int best = -2; double bw = 0;
  if (c->r1 - c->r0 > rmin){ best = -1; bw = (c->r1 - c->r0) * dR; }
  for (int J = 0; J < c->nj; J++) if (c->m[J] && c->hi[J] - c->lo[J] > hmin * 1.0001){ double w = c->m[J] * (c->hi[J] - c->lo[J]); if (w > bw){ bw = w; best = J; } }
  if (best == -2 || c->nj >= MAXJ - 2) return 0;
  if (best == -1){
    cell_t d = *c; int mid = (c->r0 + c->r1) / 2;
    d.r1 = mid; explore(&d, known); d.r0 = mid; d.r1 = c->r1; explore(&d, known);
    return 1;
  }
  /* children: the interval J halved, its m factors shared out; empty intervals are dropped */
  int J = best, m = c->m[J]; double mid = 0.5 * (c->lo[J] + c->hi[J]);
  for (int a = 0; a <= m; a++){
    cell_t d; d.r0 = c->r0; d.r1 = c->r1; d.nj = 0; d.shrunk = c->shrunk;
    for (int i = 0; i < c->nj; i++){
      if (i == J){
        if (a){ d.lo[d.nj] = c->lo[J]; d.hi[d.nj] = mid; d.m[d.nj] = a; d.nj++; }
        if (m - a){ d.lo[d.nj] = mid; d.hi[d.nj] = c->hi[J]; d.m[d.nj] = m - a; d.nj++; }
      } else if (c->m[i]){ d.lo[d.nj] = c->lo[i]; d.hi[d.nj] = c->hi[i]; d.m[d.nj] = c->m[i]; d.nj++; }
    }
    explore(&d, known);
  }
  return 1;
}
/* known = 1: the cell lies in a cell that passed, so every configuration in it passes */
static void explore(cell_t *c, int known){
  nodes++;
  double Slo = 0, Shi = 0; int top = -1;
  for (int J = 0; J < c->nj; J++) if (c->m[J]){ Slo += c->m[J] * c->lo[J]; Shi += c->m[J] * c->hi[J]; top = J; }
  double Rlo = c->r0 * dR, Rhi = c->r1 * dR;
  double z1lo = sa - Shi - Rhi - 1e-12, z1hi = sb - Slo - Rlo + 1e-12;
  double topLo = top >= 0 ? c->lo[top] : G, topHi = top >= 0 ? c->hi[top] : G;
  if (z1hi < G || z1hi < topLo) return;                       /* empty */
  int valid = z1lo > G && z1lo > topHi;
  if (!valid){
    /* crude upper bound of the measure (z1 >= G): drop if negligible */
    { double zmu = 1; int k = 1;
      for (int J = 0; J < c->nj; J++) if (c->m[J]){ double L = log(c->hi[J] / c->lo[J]) * (1 + 1e-12); zmu *= pow(L, c->m[J]) / fact[c->m[J]]; k += c->m[J]; }
      double ubp = zmu * ldexp(1.0, k) * Qup0(c->r0, c->r1) / G * (1 + 1e-9);
      if (ubp < EPS){ ndrop++; Vdrop += ubp; return; } }
    if (split(c, known, 1, HT)) return;
    /* at the smallest size: keep the part of the cell in which z1 exceeds G and the top
       factors, by lowering the top end of R and then shrinking the top interval */
    if (!c->shrunk){
      cell_t d = *c; d.shrunk = 1;
      double Shi2 = Shi;
      if (!(sa - Shi2 - d.r1 * dR - 1e-12 > G)){
        int r1 = (int)floor((sa - Shi2 - G - 2e-12) / dR); if (r1 > d.r1) r1 = d.r1;
        while (r1 > d.r0 && !(sa - Shi2 - r1 * dR - 1e-12 > G)) r1--;
        d.r1 = r1;
      }
      if (d.r1 > d.r0 && top >= 0){
        double Rh = d.r1 * dR, Soth = Shi2 - d.m[top] * d.hi[top];
        double zt = sa - Soth - Rh - 1e-12;             /* z1 + (top factors) at the worst corner */
        if (!(zt - d.m[top] * d.hi[top] > d.hi[top])){
          double h = zt / (d.m[top] + 1), h2 = (zt - G) / d.m[top]; if (h2 < h) h = h2;
          h -= 1e-9;
          if (h > d.lo[top]) d.hi[top] = h; else d.r1 = d.r0;
        }
      }
      if (d.r1 > d.r0){
        double S2 = 0; int t2 = -1; for (int J = 0; J < d.nj; J++) if (d.m[J]){ S2 += d.m[J] * d.hi[J]; t2 = J; }
        double z1l = sa - S2 - d.r1 * dR - 1e-12;
        if (z1l > G && (t2 < 0 || z1l > d.hi[t2])){ nshrink++; explore(&d, known); }
      }
    }
    npart++; Vpart += mass_est(c, z1lo, z1hi, topHi > G ? topHi : G); return;
  }
  double ub = 0, mlo = mass_lo(c, z1lo, &ub);
  if (ub < EPS){ ndrop++; Vdrop += ub; return; }
  if (!known) known = passes(c, z1lo);
  if (known){
    double me = mass_est(c, z1lo, z1hi, G);
    if (me - mlo > TAU * me){ ntight++; if (split(c, 1, 1, HT)) return; }
    npass++; Vpass += mlo; if (me > mlo) Vslack += me - mlo; return;
  }
  if (split(c, 0, RMINB, H0)) return;
  nfail++; Vfail += mlo;
}

int main(int argc, char **argv){
  if (argc < 7){ fprintf(stderr, "usage: certify_K tha thb G H0 RMINB KMAX [EPS] [ALLPASS] [TAU] [HT]\n"); return 1; }
  tha = atof(argv[1]); thb = atof(argv[2]); G = atof(argv[3]); H0 = atof(argv[4]); RMINB = atoi(argv[5]); KMAX = atoi(argv[6]);
  EPS = argc > 7 ? atof(argv[7]) : 1e-9; ALLPASS = argc > 8 ? atoi(argv[8]) : 0;
  VALID = getenv("VALIDATE") != NULL;
  TAU = argc > 9 ? atof(argv[9]) : 0.01; HT = argc > 10 ? atof(argv[10]) : H0 / 4;
  fact[0] = 1; for (int i = 1; i < 64; i++) fact[i] = fact[i - 1] * i;
  sa = 0.5 + tha; sb = 0.5 + thb;
  if (sa - 0.5 != tha || sb - 0.5 != thb){ fprintf(stderr, "depths must be dyadic\n"); return 1; }
  /* the range of n to cover is widened by 1e-12 at both ends; the caps are raised by 1e-12 */
  TLO = 1.0 / 6 - thb - 1e-12; THI = 1.0 / 3 + 1e-12;
  cap[0] = THI + 1e-12; cap[2] = (1 - 2 * TLO) / 2 + 1e-12; cap[1] = (2 - TLO - sa) / 4 + 1e-12;
  dp = G / NR; dR = dp * BW;
  CMAX = (int)ceil(2 * thb / dR) + 1; if (CMAX > NB) CMAX = NB;
  clock_t t0 = clock();
  dickman(); bands();
  double tb = (double)(clock() - t0) / CLOCKS_PER_SEC;
  /* the interval of the other factors: [G, top] dyadic, top >= sb */
  double span = sb - G; int lev = 0; while (H0 * ldexp(1.0, lev) < span) lev++;
  double topv = G + H0 * ldexp(1.0, lev);
  for (int k = 1; k <= KMAX; k++){
    cell_t c; c.shrunk = 0; c.nj = 1; c.lo[0] = G; c.hi[0] = topv; c.m[0] = k - 1; c.r0 = 0; c.r1 = NB;
    if (k == 1){ c.nj = 0; }
    explore(&c, 0);
    fprintf(stderr, "k=%d nodes=%ld pass=%ld fail=%ld part=%ld drop=%ld Vpass=%.6f Vfail=%.6f Vdrop=%.2e Vpart~%.5f Vslack~%.5f tight=%ld verr=%ld shrink=%ld\n", k, nodes, npass, nfail, npart, ndrop, Vpass, Vfail, Vdrop, Vpart, Vslack, ntight, nverr, nshrink);
  }
  double tt = (double)(clock() - t0) / CLOCKS_PER_SEC;
  if (VALID) fprintf(stderr, "validation: %ld samples, %ld violations\n", nsamp, nviol);
  double strip = 0.5 - 3 * thb;
  printf("th=[%.5f,%.5f] G=%.4f H0=%.4f RMINB=%d KMAX=%d Vpass=%.6f Vfail=%.6f strip=%.6f U=%.4f nodes=%ld verr=%ld tables=%.1fs time=%.1fs\n",
         tha, thb, G, H0, RMINB, KMAX, Vpass, Vfail, strip, Vpass / strip, nodes, nverr, tb, tt);
  return 0;
}
