**An explicit exponent for the exceptional set of the Erdős–Straus equation**, by Deep Bhattacharjee (corresponding author) and Ushashi Bhattacharya.

Let E(N) be the set of n ≤ N for which 4/n = 1/x + 1/y + 1/z has no solution in positive integers. Vaughan proved in 1970 that |E(N)| ≪ N exp(−c (log N)^{2/3}) for some c > 0, and his argument gives no value of c. For the primes in E(N) a recent sieve of Dahan gives c = 0.065, and c = 0.0943 with the constant of its sieve sum.

**Theorem A.** For every ε > 0 and all large N, |E(N)| ≤ N exp(−(c₁ − ε)(log N)^{2/3}) for all n, with

c₁ = (3/2) ((9/(8π²)) (40/49 − 6𝒥))^{1/3},

where 𝒥 is an explicit double integral with 0 ≤ 𝒥 < 0.000335, so that c₁ > 0.67916. In closed form 𝒥 ≤ 533/64680, so that c₁ ≥ (286983/(98560π²))^{1/3} = 0.6657102497…; numerically c₁ = 0.6792.

**Theorem B.** With the residue classes from the factorisations (p + 1)/4 = uaα in which one product of two of the three factors is at most √p (the three corners of the triangle of factorisations), the large sieve gives no better exponent than c₀ = 9/(2^{8/3} π^{2/3}) = 0.6607886973…

The step past c₀ uses the centre of the triangle, where u, a and α are all below √P: its classes are distinct, they agree with no class of the corner A, and with the corners B and C only along one explicit family whose contribution converges. Three strips of the centre are counted with Maynard's theorem on primes in progressions to large moduli, and 𝒥 measures the pairs for which the rule that finds a divisor in Maynard's range fails. The Erdős–Straus conjecture itself remains open.

New in this release: the bound 𝒥 < 0.000335 (Proposition IX.2) is proved by a computation in which every rounding error is directed upward. The program `checks/c/certify_J.c` cuts the depths into 625 cells, bounds the loss on each cell from above and adds the bounds, giving 𝒥 ≤ 0.0003341 in about 25 minutes on one core; its output is `checks/c/certify_J.out`. Appendix A.8 of the paper proves the method, and Table III and Figure 3 give the result by ranges of the depth.

The constructions were selected by the engine Anaxia of the BLLM project, version 3.0.5 (doi:10.5281/zenodo.23469389), in two runs: the corners among 128 constructions and the strips among 37, recorded in version 3.1.0 (doi:10.5281/zenodo.23484160). The records and the scripts that reproduce them are in `engine/`.

Lean 4 with Mathlib checks 52 statements, among them the closed-form steps of the certified bound and the final inequality c₁ > 0.67916 (axioms propext, Classical.choice, Quot.sound only); C, C++, Python and Julia programs check the corners and the centre over the 332,398 primes p ≡ 3 (mod 4) below 10⁷, the divisor chain in 128,000,000 cases, the integral 𝒥 by two independent methods, the mean values and the two bounds of the large sieve; the release workflow reruns the C and C++ checks, and `verify.yml` reruns the certified bound on every pull request.

The release also contains the paper *A power-saving asymptotic for the two-parameter Type II sieve sum of the Erdős–Straus equation*, unchanged since v2.0.0.

Files:
- `ErdosStrausExceptionalSet_final2.pdf`: the paper
- `ErdosStrausExceptionalSet_final2-tex.zip`: LaTeX source (`main.tex`, figures as PNG)
- `ErdosStrausExceptionalSet_final2-arxiv.tar.gz`: the same source, ready for arXiv
- `ErdosStrausSieveConstant_final2.pdf`, `ErdosStrausSieveConstant_final2-tex.zip`, `ErdosStrausSieveConstant_final2-arxiv.tar.gz`: the sieve-sum paper
