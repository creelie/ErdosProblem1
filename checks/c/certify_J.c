/* certify_J.c -- a certified upper bound for the integral J of "An explicit exponent for the
 * exceptional set of the Erdos-Straus equation" (Proposition prop:cert and Appendix app:cert).
 *
 * All arithmetic runs with the rounding mode FE_UPWARD, so every sum, product and quotient of the
 * program is an upper bound for the exact value; lower bounds are formed as -(upper bound of -x).
 * Grid steps are powers of two, so grid points are exact, and consecutive cells of the depth share
 * their end points exactly.
 *
 * For 0 < v < 1/42 put al = 2v, th = min(1/10 - 7v/5, 1/2 - 19v), xi = min(th - al, al), S = 1/2 + v.
 * The paper defines l(v) = int U_v(s, S - s) ds and J = int_0^{1/42} l(v) dv.  By Lemma lem:relax,
 *   l(v) <= int_0^al P(t) Q(t) dt,
 *   P(t) = int_0^t rho(t1/xi) rho((t - t1)/xi) dt1,
 *   Q(t) = int_{th - t}^{S - th} w_t(r) w_t(S - t - r) dr,
 * where w_t is the solution of r w(r) = 1_A(r) + int_{A, z < r} w(r - z) dz for
 * A = (xi, max(al - t, xi)] u [th, oo), vanishing for r <= xi.  On a cell a <= v <= b every quantity is
 * replaced by a majorant valid for all v in the cell (Lemma lem:mono): rho(t/xi) <= rho(t/xi_hi), w for
 * the larger set (xi_lo, max(2b - t, xi_hi)] u [th(b), oo), and the ranges of t, r and S widened; the
 * result is capped by the bound 1/2 + b - 2 th(b) of Proposition prop:var.
 *
 * Usage: certify_J [NH NB Q1 N2 N3 Q4 CLAIM]   (defaults 64 32 1.05 120 120 1.03 0.000335; about 25
 *        minutes on one core)
 *   NH     grid steps of w per xi (default 64)
 *   NB     fine t steps per coarse t step; the fine step is xi/1024 (default 32)
 *   Q1     ratio of the cells on [2^-18, 0.005] (default 1.05)
 *   N2 N3  numbers of equal cells on [0.005, 1/54] and [1/54, 1/44] (default 120, 120)
 *   Q4     ratio of the distances to 1/42 on [1/44, 1/42 - 1e-6] (default 1.03)
 *   CLAIM  the bound stated in the paper; the program exits with status 1 unless J < CLAIM
 * VERBOSE=1 in the environment prints the bound on every cell; certify_J.out is that output for the
 * default arguments.
 *
 * Compile with  cc -O2 -std=c99 -frounding-math -ffp-contract=off -o certify_J certify_J.c -lm
 * GCC ignores the FENV_ACCESS pragma; -frounding-math and the volatile temporaries keep the compiler
 * from folding or reordering the rounded operations.
 */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma STDC FENV_ACCESS ON

/* lower-rounded operations under FE_UPWARD */
static double dn_add(double a, double b) { volatile double x = -a; volatile double r = x - b; return -r; }
static double dn_sub(double a, double b) { volatile double r = b - a; return -r; }
static double dn_mul(double a, double b) { volatile double x = -a; volatile double r = x * b; return -r; }
static double dn_div(double a, double b) { volatile double x = -a; volatile double r = x / b; return -r; }
static double up_add(double a, double b) { volatile double r = a + b; return r; }
static double up_sub(double a, double b) { volatile double r = a - b; return r; }
static double up_mul(double a, double b) { volatile double r = a * b; return r; }
static double up_div(double a, double b) { volatile double r = a / b; return r; }

#define NG 1024            /* grid steps per unit for Dickman's function */
#define UMAX 26            /* rho tabulated on [0, UMAX] */
#define TCUT 24            /* t is cut at TCUT * xi_hi; the tail is bounded with rho(12) */

static double rhob[UMAX * NG + 1];   /* rhob[i] >= rho(i / NG) */
static double cconv[TCUT * NG + 1];  /* cconv[k] >= sum_{i + j = k} rhob[i] rhob[j] */

/* th(v) from below and from above */
static double th_dn(double v) {
    double a = dn_sub(dn_div(1, 10), up_div(up_mul(7, v), 5));
    double b = dn_sub(0.5, up_mul(19, v));
    return a < b ? a : b;
}
static double th_up(double v) {
    double a = up_sub(up_div(1, 10), dn_div(dn_mul(7, v), 5));
    double b = up_sub(0.5, dn_mul(19, v));
    return a < b ? a : b;
}
static double xi_dn(double v) {
    double a = dn_sub(th_dn(v), 2 * v);
    return a < 2 * v ? a : 2 * v;
}

static void dickman(void) {
    static double cu[UMAX * NG + 2], cd[UMAX * NG + 2];
    cu[0] = cd[0] = 0;
    for (int i = 0; i <= UMAX * NG; i++) {
        if (i <= NG) rhob[i] = 1;
        else {
            /* i rho(u_i) / NG = int_{u_i - 1}^{u_i} rho <= sum_{j = i - NG}^{i - 1} rhob[j] / NG */
            double s = up_sub(cu[i], cd[i - NG]);
            rhob[i] = up_div(s, (double)i);
        }
        cu[i + 1] = up_add(cu[i], rhob[i]);
        cd[i + 1] = dn_add(cd[i], rhob[i]);
    }
    for (int k = 0; k <= TCUT * NG; k++) {
        double s = 0;
        for (int i = 0; i <= k; i++) s = up_add(s, up_mul(rhob[i], rhob[k - i]));
        cconv[k] = s;
    }
}

/* Upper bound for the sup of w on the cells [jh, (j+1)h], 0 <= j < n, for A = (x0, w1] u [t0, oo). */
static void wtable(double x0, double w1, double t0, double h, long n, double *P, double *I, double *cu, double *cd) {
    cu[0] = cd[0] = 0;
    for (long j = 0; j < n; j++) {
        double rj = j * h, rj1 = (j + 1) * h;   /* exact */
        double p = 0;
        if (rj1 > x0) {
            double ind = ((rj <= w1) || (rj1 >= t0)) ? 1 : 0;
            double s = ind;
            /* v in [rj - w1, rj1 - x0] */
            double lo = dn_sub(rj, w1), hi = up_sub(rj1, x0);
            if (lo < 0) lo = 0;
            if (hi > 0) {
                long mlo = (long)floor(lo / h), mhi = (long)ceil(hi / h) - 1;
                if (mhi >= j) { fprintf(stderr, "grid not explicit\n"); exit(2); }
                if (mhi >= mlo) s = up_add(s, up_sub(cu[mhi + 1], cd[mlo]));
            }
            /* v in (0, rj1 - t0] */
            if (rj1 > t0) {
                double hi2 = up_sub(rj1, t0);
                long mhi = (long)ceil(hi2 / h) - 1;
                if (mhi >= j) { fprintf(stderr, "grid not explicit (2)\n"); exit(2); }
                if (mhi >= 0) s = up_add(s, cu[mhi + 1]);
            }
            double r0 = rj > x0 ? rj : x0;
            p = up_div(s, r0);
        }
        P[j] = p;
        I[j] = up_mul(h, p);
        cu[j + 1] = up_add(cu[j], I[j]);
        cd[j + 1] = dn_add(cd[j], I[j]);
    }
}

static double pow2_floor(double x) { int e; frexp(x, &e); return ldexp(1.0, e - 1); } /* largest 2^k <= x */

static int NH = 64, NB = 32;
#define CLAIM "0.000335"
static double *P, *I, *CU, *CD;
static long NMAX = 0;

/* upper bound for int_a^b l(v) dv */
static double cell(double a, double b, int verbose) {
    double xlo = xi_dn(a), xb = xi_dn(b);
    if (xb < xlo) xlo = xb;
    double xhi = up_mul(2, b), x2 = up_sub(th_up(a), 2 * a);
    if (x2 < xhi) xhi = x2;
    double tlo = th_dn(b);
    double Sa = dn_add(0.5, a), Sb = up_add(0.5, b);
    double al_hi = 2 * b;
    if (!(xlo > 0) || !(al_hi < tlo)) { fprintf(stderr, "bad cell %g %g\n", a, b); exit(2); }

    double h = pow2_floor(xlo / NH);
    long n = (long)ceil(Sb / h) + 2;
    if (n > NMAX) {
        NMAX = n;
        P = realloc(P, n * sizeof(double)); I = realloc(I, n * sizeof(double));
        CU = realloc(CU, (n + 1) * sizeof(double)); CD = realloc(CD, (n + 1) * sizeof(double));
    }
    double gt = up_mul(1.0 / NG, xhi);       /* fine t step, >= xhi / NG */
    double Gt = gt * NB;                      /* coarse t step (exact) */
    double tmax = al_hi;
    int tail = 0;
    if (tmax > TCUT * xhi) { tmax = TCUT * xhi; tail = 1; }
    long K = (long)ceil(tmax / Gt);
    if ((K * NB) > TCUT * NG) K = TCUT * NG / NB;
    if (dn_mul((double)K, Gt) < tmax) { fprintf(stderr, "t grid short\n"); exit(2); }

    double L = 0, wprev = -1;
    for (long k = 0; k < K; k++) {
        double tau0 = dn_mul((double)k, Gt), tau1 = up_mul((double)(k + 1), Gt);
        double w1 = up_sub(al_hi, tau0);
        if (w1 < xhi) w1 = xhi;
        if (w1 != wprev) { wtable(xlo, w1, tlo, h, n, P, I, CU, CD); wprev = w1; }
        /* band integral of P over [tau0, tau1]: gt^2/2 sum (c_k + c_{k-1}) over the fine bands */
        double band = 0;
        for (long f = k * NB; f < (k + 1) * NB; f++) {
            band = up_add(band, cconv[f]);
            if (f > 0) band = up_add(band, cconv[f - 1]);
        }
        band = up_mul(up_mul(up_mul(gt, gt), 0.5), band);
        /* Q over r in [th(b) - tau1, Sb - th(b)]: U vanishes unless th <= s <= S - th */
        double rlo = dn_sub(tlo, tau1), rhi = up_sub(Sb, tlo);
        long m0 = rlo > 0 ? (long)floor(rlo / h) : 0, m1 = (long)ceil(rhi / h) - 1;
        if (m1 > n - 1) m1 = n - 1;
        double sa = dn_sub(Sa, tau1), sb = up_sub(Sb, tau0);
        double Q = 0;
        for (long m = m0; m <= m1; m++) {
            if (I[m] == 0) continue;
            double vlo = dn_sub(sa, (m + 1) * h), vhi = up_sub(sb, m * h);
            long q0 = (long)floor(vlo / h), q1 = (long)ceil(vhi / h) - 1;
            if (q0 < 0) q0 = 0;
            if (q1 > n - 1) { fprintf(stderr, "v grid short\n"); exit(2); }
            double mx = 0;
            for (long q = q0; q <= q1; q++) if (P[q] > mx) mx = P[q];
            Q = up_add(Q, up_mul(I[m], mx));
        }
        L = up_add(L, up_mul(band, Q));
    }
    if (tail) {
        /* int_{t > 24 xhi} P <= 2 (xhi int_12^oo rho)(xhi e^gamma) <= 2 * (13/12) * 1.7811 rho(12) xhi^2,
           Q <= Sb / xlo^2 since w <= 1 / xlo */
        double r12 = rhob[12 * NG];
        double tl = up_mul(up_mul(up_mul(2 * 1.7811, up_div(13, 12)), r12), up_mul(xhi, xhi));
        tl = up_mul(tl, up_div(Sb, up_mul(xlo, xlo)));
        L = up_add(L, tl);
    }
    /* l(v) <= 1/2 + v - 2 th(v) <= 1/2 + b - 2 th(b) (Proposition var) */
    double r1 = up_sub(Sb, up_mul(2, tlo));
    if (r1 < L) L = r1;
    double res = up_mul(up_sub(b, a), L);
    if (verbose) printf("cell [%.9f, %.9f] xi in [%.6g, %.6g] h=%.3g n=%ld K=%ld  lbar=%.6f  int<=%.4e\n",
                        a, b, xlo, xhi, h, n, K, L, res);
    return res;
}

int main(int argc, char **argv) {
    /* the claim is read rounded down, so that J < claim (in binary) implies J < CLAIM (in decimal) */
    fesetround(FE_DOWNWARD);
    double claim = strtod(argc > 7 ? argv[7] : CLAIM, NULL);
    fesetround(FE_UPWARD);
    int verbose = getenv("VERBOSE") != NULL;
    NH = argc > 1 ? atoi(argv[1]) : 64;
    NB = argc > 2 ? atoi(argv[2]) : 32;
    double q1 = argc > 3 ? atof(argv[3]) : 1.05;
    int n2 = argc > 4 ? atoi(argv[4]) : 120;
    int n3 = argc > 5 ? atoi(argv[5]) : 120;
    double q4 = argc > 6 ? atof(argv[6]) : 1.03;
    if (NH < 16 || (NB & (NB - 1)) || NB < 1 || NB > 1024 || !(q1 > 1) || n2 < 1 || n3 < 1 || !(q4 > 1)) {
        fprintf(stderr, "usage: certify_J [NH>=16 NB(power of 2) Q1>1 N2 N3 Q4>1 CLAIM]\n");
        return 2;
    }
    dickman();
    printf("rho(2) <= %.10f (exact 1 - log 2 = 0.3068528194), rho(3) <= %.10f (0.0486083883), rho(12) <= %.3e\n",
           rhob[2 * NG], rhob[3 * NG], rhob[12 * NG]);

    /* break points (any doubles: the cells tile [a0, e1] exactly) */
    double a0 = ldexp(1.0, -18), v2 = 0.005, v54 = 1.0 / 54, v44 = 1.0 / 44;
    double e42 = up_div(1, 42), eps = 1e-6;
    double total = 0, part[6] = {0};
    int ncell[6] = {0};

    /* [0, a0]: l(v) <= 1/2 + v - 2 th(v) <= 1/2 + a0 - 2 th(a0) */
    part[0] = up_mul(a0, up_sub(up_add(0.5, a0), up_mul(2, th_dn(a0))));

    /* [a0, v2]: geometric cells of ratio q1 */
    {
        double a = a0;
        while (a < v2) {
            double b = a * q1;
            if (b > v2 || v2 - b < 0.1 * (b - a)) b = v2;
            part[1] = up_add(part[1], cell(a, b, verbose)); ncell[1]++;
            a = b;
        }
    }
    for (int i = 0; i < n2; i++) {
        double a = v2 + (v54 - v2) * i / n2, b = (i + 1 == n2) ? v54 : v2 + (v54 - v2) * (i + 1) / n2;
        part[2] = up_add(part[2], cell(a, b, verbose)); ncell[2]++;
    }
    for (int i = 0; i < n3; i++) {
        double a = v54 + (v44 - v54) * i / n3, b = (i + 1 == n3) ? v44 : v54 + (v44 - v54) * (i + 1) / n3;
        part[3] = up_add(part[3], cell(a, b, verbose)); ncell[3]++;
    }
    /* [v44, e1]: geometric in d = 1/42 - v */
    double e1;
    {
        double d0 = 1.0 / 42 - v44;
        double d = d0, a = v44;
        while (d > eps) {
            double dn = d / q4;
            if (dn < eps) dn = eps;
            double b = 1.0 / 42 - dn;
            part[4] = up_add(part[4], cell(a, b, verbose)); ncell[4]++;
            a = b; d = dn;
        }
        e1 = a;
    }
    /* [e1, 1/42]: l(v) <= 1/2 + v - 2 th(v) <= 1/2 + 1/42 - 2 th(1/42) */
    part[5] = up_mul(up_sub(e42, e1), up_sub(up_add(0.5, e42), up_mul(2, th_dn(e42))));
    for (int i = 0; i < 6; i++) total = up_add(total, part[i]);
    const char *name[6] = {"[0, 2^-18]", "[2^-18, 0.005]", "[0.005, 1/54]", "[1/54, 1/44]",
                           "[1/44, 1/42 - 1e-6]", "[1/42 - 1e-6, 1/42]"};
    for (int i = 0; i < 6; i++) printf("%-22s %4d cells  %.9e\n", name[i], ncell[i], part[i]);
    printf("J <= %.9e\n", total);
    if (!(total < claim)) {
        printf("the bound %.9e is not below the claim %.9e\n", total, claim);
        return 1;
    }
    printf("certified: J < %s\n", argc > 7 ? argv[7] : CLAIM);
    return 0;
}
