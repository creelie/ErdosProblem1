**An explicit exponent for the exceptional set of the Erdős–Straus equation**, by Deep Bhattacharjee (corresponding author) and Ushashi Bhattacharya.

Let E(N) be the set of n ≤ N for which 4/n = 1/x + 1/y + 1/z has no solution in positive integers. Vaughan proved in 1970 that |E(N)| ≪ N exp(−c (log N)^{2/3}) for some c > 0, and his argument gives no value of c. For the primes in E(N) a recent sieve of Dahan gives c = 0.065, and c = 0.0943 with the constant of its sieve sum.

**Theorem A.** For every ε > 0 and all large N, |E(N)| ≤ N exp(−(c₁ − ε)(log N)^{2/3}) for all n, with

c₁ = (3/2) ((9/(8π²)) (40/49 − 6𝒥 + 6𝒦))^{1/3},

where 𝒥 is an explicit double integral with 0 ≤ 𝒥 < 0.000335 and 𝒦 is the measure of the factorable pairs of the deep strips, with 𝒦 > 0.00859, so that c₁ > 0.6931. In closed form c₁ ≥ (286983/(98560π²))^{1/3} = 0.6657102497…; a Monte Carlo simulation gives 𝒦 ≈ 0.0138 and c₁ ≈ 0.7015.

**Theorem B.** With the residue classes from the factorisations (p + 1)/4 = uaα in which one product of two of the three factors is at most √p (the three corners of the triangle of factorisations), the large sieve gives no better exponent than c₀ = 9/(2^{8/3} π^{2/3}) = 0.6607886973…

The step past c₀ uses the centre of the triangle, where u, a and α are all below √P: its classes are distinct, they agree with no class of the corner A, and with the corners B and C only along one explicit family whose contribution converges. Three strips of the centre are counted with Maynard's theorem on primes in progressions to large moduli, and 𝒥 measures the pairs for which the rule that finds a divisor in Maynard's range fails. The Erdős–Straus conjecture itself remains open.

New in this release: the deep strips. Behind the strips of depth 1/42, Maynard's second theorem (primes in progressions to moduli up to p^{3/5}, Mem. AMS 306 (2025), no. 1543) counts the pairs whose modulus splits into three factors of suitable sizes; the large prime factors are placed in three groups and the small ones fill two of them to prescribed sizes. No splitting exists beyond depth 1/10 (Lemma VII.4), which is the level p^{3/5}. The bound 𝒦 > 0.00859 is proved by `checks/c/certify_K.c`, which bounds the density of the factorable pairs from below on 47 cells of depth with every rounding directed against the bound (about five hours on one core; output `checks/c/certify_K.out`), and `checks/python/sum_K.py` adds up the cells in exact rational arithmetic. This takes c₁ past 0.6797, the value with every triple of the strips of v3.1.0. The paper also cites T. Newton's preprint on the exact divisor tests for the same primes (doi:10.5281/zenodo.23760689), whose families are classes to fixed moduli and leave the exponent unchanged.

The constructions were selected by the engine Anaxia of the BLLM project in three runs: with version 3.0.5 (doi:10.5281/zenodo.23469389) the corners among 128 constructions and the strips among 37, and with version 3.1.0 (doi:10.5281/zenodo.23484160) the splitting of the moduli in the deep strips among 15. The records and the scripts that reproduce them are in `engine/`.

Lean 4 with Mathlib checks 59 statements, among them the bound behind the depth 1/10, the areas of the deep strips and the final inequality c₁ > 0.6931 (axioms propext, Classical.choice, Quot.sound only); C, C++, Python and Julia programs check the corners and the centre over the 332,398 primes p ≡ 3 (mod 4) below 10⁷, the divisor chain in 128,000,000 cases, the integral 𝒥 by two independent methods, the mean values and the two bounds of the large sieve; the release workflow reruns the C and C++ checks, `verify.yml` reruns the certified bound for 𝒥 and the exact sum of the cells of 𝒦 on every pull request.

The release also contains the paper *A power-saving asymptotic for the two-parameter Type II sieve sum of the Erdős–Straus equation*, unchanged since v2.0.0.

Files:
- `ErdosStrausExceptionalSet_final2.pdf`: the paper
- `ErdosStrausExceptionalSet_final2-tex.zip`: LaTeX source (`main.tex`, figures as PNG)
- `ErdosStrausExceptionalSet_final2-arxiv.tar.gz`: the same source, ready for arXiv
- `ErdosStrausSieveConstant_final2.pdf`, `ErdosStrausSieveConstant_final2-tex.zip`, `ErdosStrausSieveConstant_final2-arxiv.tar.gz`: the sieve-sum paper
