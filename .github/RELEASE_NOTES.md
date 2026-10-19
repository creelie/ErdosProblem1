**An explicit exponent for the exceptional set of the Erdős–Straus equation**, by Deep Bhattacharjee (corresponding author) and Ushashi Bhattacharya.

Let E(N) be the set of n ≤ N for which 4/n = 1/x + 1/y + 1/z has no solution in positive integers. Vaughan proved in 1970 that |E(N)| ≪ N exp(−c (log N)^{2/3}) for some c > 0, and his argument gives no value of c. For the primes in E(N) a recent sieve of Dahan gives c = 0.065, and c = 0.0943 with the constant of its sieve sum.

**Theorem A.** For every ε > 0 and all large N, |E(N)| ≤ N exp(−(c₀ − ε)(log N)^{2/3}) with c₀ = 9/(2^{8/3} π^{2/3}) = 0.6607886973…, for all n and not only the primes.

**Theorem B.** With the residue classes of the Type II identities used for Theorem A, the large sieve gives no better exponent than c₀.

The three corners of the triangle of factorisations and their normal forms were selected among 128 candidate constructions by the engine Anaxia of the BLLM project, version 3.0.5 (doi:10.5281/zenodo.23469389); the record of the run and a script that reproduces it are in `engine/`. The Erdős–Straus conjecture itself remains open.

Lean 4 with Mathlib checks the algebraic steps (axioms propext, Classical.choice, Quot.sound only); C, C++, Python and Julia programs check the corners over the 332,398 primes p ≡ 3 (mod 4) below 10⁷, the six-parameter family of coincidences, the mean values and the two bounds of the large sieve.

The release also contains the revised paper of v1.0.0, now titled *A power-saving asymptotic for the two-parameter Type II sieve sum of the Erdős–Straus equation*. Its leading constant 9/(2π²) and linear term were found earlier by Z. Liu (doi:10.5281/zenodo.22650516); the revision credits Liu, states that the constant term and the error O(J^{-1/2} log J) are what the paper adds, and records the archive DOI of v1.0.0.

Files:
- `ErdosStrausExceptionalSet_final2.pdf`: the paper
- `ErdosStrausExceptionalSet_final2-tex.zip`: LaTeX source (`main.tex`, figures as PNG)
- `ErdosStrausExceptionalSet_final2-arxiv.tar.gz`: the same source, ready for arXiv
- `ErdosStrausSieveConstant_final2.pdf`, `ErdosStrausSieveConstant_final2-tex.zip`, `ErdosStrausSieveConstant_final2-arxiv.tar.gz`: the revised paper of v1.0.0
