# check.jl -- independent check in Julia (no packages).
# 1. Exact rational check of the local factor 1 + 2/(p(p-1)) - (p+1)/((p-1)p^2) = 1 + 1/p^2.
# 2. W(J) for J <= 2000 against K((log J + gamma + T)^2 + M), K = 9/(2 pi^2).
# 3. The ratio K / (3/(2 pi^2)) = 3 and the exponent (1/2)(K/(25e))^(1/3) = 0.0943...

function sieve(n)
    s = trues(n); s[1] = false
    for i in 2:isqrt(n)
        if s[i]
            for j in i*i:i:n
                s[j] = false
            end
        end
    end
    findall(s)
end

function main()
    ok = true
    ps = sieve(10^6)
    for p in ps[2:end]
        p > 5000 && break
        q = big(p)
        lhs = 1 + 2 // (q * (q - 1)) - (q + 1) // ((q - 1) * q^2)
        lhs == 1 + 1 // q^2 || (println("local factor mismatch at ", p); ok = false)
    end
    println("1. local factors exact for odd p < 5000")
    K = 9 / (2 * pi^2)
    γ = Base.MathConstants.eulergamma
    T = log(2) / 3 + sum(log(p) / ((p - 1) * (p^2 + 1)) for p in Float64.(ps[2:end]))
    M = -4 / 9 * log(2)^2 - sum(log(p)^2 * ((p + 1) / ((p - 1) * (p^2 + 1)) + 1 / ((p - 1)^2 * (p^2 + 1)^2)) for p in Float64.(ps[2:end]))
    M -= (log(1e6) + 1) / 1e6  # primes above 10^6, by the prime number theorem
    abs(M + 1.0429405079) < 1e-7 || (println("M = ", M); ok = false)
    J = 2000
    ph = collect(0:J)
    for p in 2:J
        if ph[p+1] == p
            for k in p:p:J
                ph[k+1] -= ph[k+1] ÷ p
            end
        end
    end
    W = 0.0
    for m in 1:J
        for x in 1:m
            gcd(m, x) == 1 || continue
            t = 1 / ((iseven(m * x) ? 4 : 2) * ph[m+1] * ph[x+1])
            W += x == m ? t : 2t
        end
        if m in (100, 500, 2000)
            L = log(m)
            f = K * ((L + γ + T)^2 + M)
            println("2. J = $m  W = $(round(W, digits=9))  formula = $(round(f, digits=9))  error = $(W - f)")
            abs(W - f) < 3L / sqrt(m) || (ok = false)
        end
    end
    r = K / (3 / (2 * pi^2))
    c = 0.5 * (K / (25 * exp(1)))^(1 / 3)
    println("3. K/(3/(2pi^2)) = $r, exponent = $c")
    (abs(r - 3) < 1e-12 && abs(c - 0.0943) < 1e-4) || (ok = false)
    println(ok ? "Julia check passed" : "Julia check FAILED")
    exit(ok ? 0 : 1)
end

main()
