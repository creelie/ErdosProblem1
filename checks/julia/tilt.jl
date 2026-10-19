# tilt.jl -- the two bounds of the tilt lemma in "An explicit exponent for the exceptional set of the
# Erdos-Straus equation", in the model N(lambda) = a lambda^2.  Then
#   F(s) = 2a/s^2,  mu(s) = 4a/s^3,  V(s) <= 12a/s^4.
# Upper bound (Rankin):      log L(Q) <= min_s (s T + F(s)),            T = log Q.
# Lower bound (tilted measure): log L(Q) >= F(s) + s (T - 2 sqrt(2V)) - log 2  at mu(s) + sqrt(2V(s)) = T.
# Both, divided by (3/2)(4a)^(1/3) T^(2/3), tend to 1; the program checks this for the four values of a
# of the paper (one, two and three corners, and the whole triangle) and prints the constants
# c = (3/2) a^(1/3) of the bound N exp(-(c - eps)(log N)^(2/3)).
using Printf

function upper(a, T)
    lo, hi = 1e-12, 10.0
    f(s) = s * T + 2a / s^2
    for _ in 1:300                       # golden-section search on a convex function
        m1 = lo + (hi - lo) * 0.381966; m2 = hi - (hi - lo) * 0.381966
        f(m1) < f(m2) ? (hi = m2) : (lo = m1)
    end
    f((lo + hi) / 2)
end

function lower(a, T)
    g(s) = 4a / s^3 + sqrt(24a / s^4) - T   # decreasing in s
    lo, hi = 1e-12, 1e3
    for _ in 1:400
        mid = (lo + hi) / 2
        g(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    s = (lo + hi) / 2
    2a / s^2 + s * (T - 2sqrt(24a / s^4)) - log(2)
end

kappa = 9 / (2pi^2)
fails = 0
for (name, a) in (("one corner", kappa / 16), ("two corners", kappa / 8),
                  ("three corners", 3kappa / 16), ("whole triangle", kappa / 4))
    exact(T) = 1.5 * (4a)^(1 / 3) * T^(2 / 3)
    ru = Float64[]; rl = Float64[]
    for T in (1e2, 1e4, 1e6, 1e9, 1e12, 1e15)
        push!(ru, upper(a, T) / exact(T)); push!(rl, lower(a, T) / exact(T))
    end
    @printf("%-15s a = %.10f  c = %.6f  upper/exact %s  lower/exact %s\n", name, a, 1.5 * a^(1 / 3),
            join([@sprintf("%.6f", x) for x in ru], " "), join([@sprintf("%.4f", x) for x in rl], " "))
    global fails += !all(abs.(ru .- 1) .< 1e-6)
    global fails += !(issorted(rl) && rl[end] > 0.99 && rl[end] < 1)
end
c0 = 9 / (2^(8 / 3) * pi^(2 / 3))
fails += abs(c0 - 1.5 * (3kappa / 16)^(1 / 3)) > 1e-14
fails += abs(c0^3 - 729 / (256pi^2)) > 1e-14
@printf("c0 = 9/(2^(8/3) pi^(2/3)) = %.16f,  c0^3 = 729/(256 pi^2) = %.16f\n", c0, 729 / (256pi^2))
fails == 0 || error("tilt.jl: $fails failures")
println("tilt.jl: all checks passed")
