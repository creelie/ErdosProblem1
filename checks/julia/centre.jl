# centre.jl -- independent check, in Julia, of the constant of the centre strips in "An explicit
# exponent for the exceptional set of the Erdos-Straus equation".
#
# For a depth delta in (0, 1/42), al = 2 delta, be = min(1/10 - 7 delta/5, 1/2 - 19 delta) and
# ga = min(be - al, al), a pair (u, a) with u a ~ x^(1/2 + delta) is lost when u a has no prime
# factor in (x^al, x^be) and q T <= x^al, where T is the x^ga-smooth part of u a and q its largest
# prime factor in (x^ga, x^al].  The loss is computed here by a method different from the Python
# program: the density b_M of products of primes in (ga, M] is marched on a uniform grid, the
# be-rough part has the Buchstab density, the two are convolved explicitly, and the integrals are
# Riemann sums (trapezoid in s1, piecewise linear cubature on the triangle t1 + t2 <= al, midpoints
# in delta).  The method is of first order in the grid step, so the check compares the losses at
# five depths with the values of the Python quadrature, and the constant at a moderate resolution.
#
# Usage: julia centre.jl
# Dickman rho and Buchstab omega on a fine grid (Simpson-type cumulative integration)
const H0 = 1.0e-5; const UMAX = 40.0
const N0 = round(Int, UMAX / H0) + 1; const K1 = round(Int, 1 / H0)
function tables()
    x = [i * H0 for i in 0:N0-1]
    rho = ones(N0)
    for i in K1+2:N0
        # rho' = -rho(u-1)/u, trapezoid step
        rho[i] = rho[i-1] - H0 / 2 * (rho[i-1-K1] / x[i-1] + rho[i-K1] / x[i])
    end
    om = zeros(N0)
    for i in K1+1:min(2K1 + 1, N0); om[i] = 1 / x[i]; end
    for i in 2K1+2:N0
        # (u omega)' = omega(u-1)
        om[i] = (x[i-1] * om[i-1] + H0 / 2 * (om[i-1-K1] + om[i-K1])) / x[i]
    end
    x, rho, om
end
const X0, RHO, OM = tables()
interp(tab, u) = begin
    v = u / H0; i = floor(Int, v); f = v - i
    i >= N0 - 1 ? 0.0 : (1 - f) * tab[i+1] + f * tab[i+2]
end
rhof(u) = u <= 1 ? 1.0 : (u <= 2 ? 1 - log(u) : interp(RHO, u))
omf(u) = u < 1 ? 0.0 : (u <= 2 ? 1 / u : interp(OM, u))
tau(d) = min(0.1 - 1.4d, 0.5 - 19d)

# density of the X^ga-smooth part, uncovered mass for one delta, as a function of s1
function loss_delta(d, M_PER)
    al = 2d; be = tau(d); S = 0.5 + d
    ga = min(be - al, al)
    h = ga / M_PER
    n = ceil(Int, S / h) + 2
    r = [i * h for i in 0:n]
    K = [r[i] >= be ? omf(r[i] / be) / be : 0.0 for i in 1:n+1]
    D = [rhof(r[i] / ga) for i in 1:n+1]
    tmax = min(al, 24ga); nt = floor(Int, tmax / h + 1e-9)
    # phi tables for M = max(al - t, ga), t = k h, k = 0..nt
    phis = Vector{Vector{Float64}}(undef, nt + 1)
    for k in 0:nt
        M = max(al - k * h, ga)
        b = zeros(n + 1); C = zeros(n + 1)   # C = cumulative integral of b
        jg = round(Int, ga / h); jM = round(Int, M / h)
        for i in 1:n
            y = r[i+1]
            ind = (i > jg && i <= jM) ? 1.0 : 0.0
            # integral of b(y - z) for z in (ga, min(M, y)) = C(y - ga) - C(y - min(M,y))
            lo = i - min(jM, i); hi = i - jg
            q = hi > 0 ? C[hi+1] - C[lo+1] : 0.0
            b[i+1] = (ind + q) / y
            C[i+1] = C[i] + h / 2 * (b[i] + b[i+1])
        end
        # phi = b + K + b*K
        phi = similar(b)
        for i in 0:n
            conv = 0.0
            for j in 1:i-1
                bj = b[j+1]; kk = K[i-j+1]
                (bj == 0.0 || kk == 0.0) && continue
                conv += bj * kk
            end
            phi[i+1] = b[i+1] + K[i+1] + conv * h
        end
        phis[k+1] = phi
    end
    # integrate over s1 in (be, S - be) on a uniform grid; the (t1, t2) triangle uses the
    # piecewise linear cubature (weights h^2 inside, h^2/2 on the edges, h^2/6 and h^2/3 at the corners)
    ns = max(200, ceil(Int, (S - 2be) / h))
    hs = (S - 2be) / ns
    edge = abs(nt * h - tmax) < 1e-12 * max(1.0, tmax)
    wt(k1, k2) = begin
        if !edge
            (k1 == 0 ? 0.5 : 1.0) * (k2 == 0 ? 0.5 : 1.0)
        else
            n0 = (k1 == 0) + (k2 == 0) + (k1 + k2 == nt)
            n0 == 0 ? 1.0 : (n0 == 1 ? 0.5 : ((k1 == 0 && k2 == 0) ? 1 / 6 : 1 / 3))
        end
    end
    tot = 0.0
    for i in 0:ns
        s1 = be + i * hs; s2 = S - s1
        g = 0.0
        for k1 in 0:nt, k2 in 0:(nt - k1)
            ph = phis[k1 + k2 + 1]
            a1 = s1 - k1 * h; a2 = s2 - k2 * h
            v1 = a1 / h; j1 = floor(Int, v1); f1 = v1 - j1
            v2 = a2 / h; j2 = floor(Int, v2); f2 = v2 - j2
            p1 = (1 - f1) * ph[j1+1] + f1 * ph[j1+2]
            p2 = (1 - f2) * ph[j2+1] + f2 * ph[j2+2]
            g += wt(k1, k2) * D[k1+1] * D[k2+1] * p1 * p2
        end
        g *= h * h
        tot += (i == 0 || i == ns ? 0.5 : 1.0) * g
    end
    tot * hs
end


const REF = [(0.01, 3.097365e-3), (0.02, 1.853518e-2), (0.022, 2.520337e-2),
             (0.023, 4.280047e-2), (0.0235, 1.716346e-1)]   # Python quadrature, 8 nodes, 32 steps

function main()
    ok = true
    println("loss at five depths, grid step ga/64, against the Python quadrature")
    for (d, ref) in REF
        v = loss_delta(d, 64)
        rel = abs(v - ref) / ref
        println("  delta = ", d, ": Julia ", round(v, sigdigits = 6), ", Python ", ref, ", relative difference ",
                round(rel, sigdigits = 2))
        ok &= rel < 0.015
    end
    kap = 9 / (2pi^2); Iceil = 1 / 84 - 3 / 3528
    J = 0.0; ND = 10
    for (a, b) in ((0.0, 1 / 54), (1 / 54, 1 / 44), (1 / 44, 1 / 42))
        for i in 1:ND
            J += loss_delta(a + (i - 0.5) * (b - a) / ND, 16) * (b - a) / ND
        end
    end
    c = 1.5 * ((kap / 4) * (0.75 + 6 * (Iceil - J)))^(1 / 3)
    println("J = ", round(J, sigdigits = 6), " (grid step ga/16, ten midpoints per window), c1 = ", round(c, digits = 6))
    ok &= 0.6791 < c < 0.6794
    ok || error("centre check failed")
    println("OK: the Julia losses agree with the Python quadrature, c1 = 0.6792")
end

main()
