# Erdős–Straus: an explicit exponent for the exceptional set and the Type II sieve sum

Deep Bhattacharjee (corresponding author, [ORCID 0000-0003-0466-750X](https://orcid.org/0000-0003-0466-750X)) and Ushashi Bhattacharya ([ORCID 0009-0002-2254-3914](https://orcid.org/0009-0002-2254-3914)).

Contact: d.bhattacharjee@erl-forschung.de (corresponding), itsdeep@live.com; for questions about this repository, egspl@proton.me (GitHub [creelie](https://github.com/creelie)).

The Erdős–Straus conjecture says that 4/n = 1/x + 1/y + 1/z has a solution in positive integers for every n ≥ 2. It is open, and neither paper here proves it.

| Paper | PDF |
|---|---|
| An explicit exponent for the exceptional set of the Erdős–Straus equation | [`preprintErdosStrausExceptionalSet/ErdosStrausExceptionalSet.pdf`](preprintErdosStrausExceptionalSet/ErdosStrausExceptionalSet.pdf) |
| A power-saving asymptotic for the two-parameter Type II sieve sum of the Erdős–Straus equation | [`preprintErdosStrausSieveConstant/ErdosStrausSieveConstant.pdf`](preprintErdosStrausSieveConstant/ErdosStrausSieveConstant.pdf) |

## An explicit exponent for the exceptional set

Let E(N) be the set of n ≤ N for which the equation has no solution. Vaughan proved in 1970 that |E(N)| ≪ N exp(−c (log N)^{2/3}) for some c > 0, and his argument gives no value of c. For the primes in E(N), which all satisfy n ≡ 1 (mod 24), Dahan's sieve ([arXiv:2608.24035](https://arxiv.org/abs/2608.24035)) gives c = 0.065, and c = 0.0943 with the constant 9/(2π²) of its sieve sum found by Liu.

For the same primes T. Newton ([doi:10.5281/zenodo.23760689](https://doi.org/10.5281/zenodo.23760689)) studies the exact divisor tests attached to a prescribed first denominator: his finite menu of Type II congruences covers 25651/26496 of their reduced residue classes. His families are classes to fixed moduli, so they leave the exponent unchanged; the paper here sieves with classes to moduli that grow with N.

**Theorem A.** For every ε > 0 and all large N,

|E(N)| ≤ N exp(−(c₁ − ε)(log N)^{2/3}),  c₁ = (3/2) ((9/(8π²)) (40/49 − 6𝒥 + 6𝒦))^{1/3},

where 𝒥 is an explicit double integral and 𝒦 the measure of the factorable pairs of the deep strips (Section 7 of the paper). Moreover 0 ≤ 𝒥 < 0.000335 and 𝒦 > 0.00859, so that c₁ > 0.6931, and in closed form 𝒥 ≤ 533/64680 and 𝒦 ≥ 0, so that

c₁ ≥ (286983/(98560 π²))^{1/3} = 0.6657102497…

Numerically 𝒥 = 0.0003, and a Monte Carlo simulation gives 𝒦 ≈ 0.0138 and c₁ ≈ 0.7015. The bound 𝒥 < 0.000335 is proved by `checks/c/certify_J.c`, which bounds the integral cell by cell with every rounding directed upward, and the bound 𝒦 > 0.00859 by `checks/c/certify_K.c`, which bounds the density of the factorable pairs from below on 47 cells of depth with every rounding directed against the bound; `checks/python/sum_K.py` adds up the cells in exact rational arithmetic.

**Theorem B.** With the residue classes from the factorisations (p + 1)/4 = uaα in which one product of two of the three factors is at most √p, the large sieve gives no better exponent than c₀ = 9/(2^{8/3} π^{2/3}) = 0.6607886973…

The factorisations (p + 1)/4 = uaα fill a triangle. Its three corners, where one product of two variables is at most √p, give the classes that the Bombieri–Vinogradov theorem can count; each corner has a normal form of the class (u/a, 4u²α or (4a²α)⁻¹) that separates its classes, and the corners alone give exactly c₀ (Theorem B). The gain to c₁ comes from the centre of the triangle, where u, a and α are all below √P. Its classes are distinct and agree with no class of the corner A, and with B or C only along one explicit family whose contribution converges. Three strips of the centre along the corners, of depth 1/42, are counted with Maynard's theorem on primes in progressions to large moduli (Mem. AMS 306 (2025), no. 1542), using a rule that finds a divisor of the modulus in Maynard's range; 𝒥 measures the pairs where the rule fails. With every triple of the strips the constant would be 0.6797. Behind the strips, from depth 1/42 to depth 1/10, Maynard's second theorem (primes in progressions to moduli up to p^{3/5}, Mem. AMS 306 (2025), no. 1543) counts the pairs whose modulus splits into three factors of suitable sizes; the large prime factors are placed in three groups and the small ones fill two of them to prescribed sizes. No splitting exists beyond depth 1/10, which is the level p^{3/5}. These deep strips add 6𝒦 to the share and take c₁ past 0.6797. With the whole centre (which needs the Elliott–Halberstam conjecture) the constant would be 0.7273.

The construction was selected by the engine Anaxia of the BLLM project in three runs: with version 3.0.5 ([doi:10.5281/zenodo.23469389](https://doi.org/10.5281/zenodo.23469389)) the corners among 128 constructions and the strips among 37, and with version 3.1.0 ([doi:10.5281/zenodo.23484160](https://doi.org/10.5281/zenodo.23484160)) the splitting of the moduli in the deep strips among 15; the records are in `engine/`.

## The Type II sieve sum

B. Dahan ([arXiv:2608.24035](https://arxiv.org/abs/2608.24035)) bounds the primes n ≡ 1 (mod 24) whose Type II solutions are not found by a search of depth J through a sieve governed by

W(J) = Σ_{u, a ≤ J, gcd(u, a) = 1} 1/φ(4ua),

proved W(J) ≥ (3/(2π²) + o(1)) (log J)², and asked for W(J) ≥ (0.45 + o(1)) (log J)².

**Theorem.** W(J) = κ((log J + γ + T)² + M) + O(J^{-1/2} log J), with κ = 9/(2π²) = 0.4559453…, T = 0.3131549610… and M = −1.0429405079…, both given by explicit sums over primes.

The first two terms, κ(log J)² and 2κ(γ + T) log J with an error O(1), were found earlier by Z. Liu ([doi:10.5281/zenodo.22650516](https://doi.org/10.5281/zenodo.22650516)), who also derived the exponent 0.0943 for Dahan's bound. The theorem adds the constant term κ((γ + T)² + M) = −0.1140687011… and the error O(J^{-1/2} log J).

## Contents

| Path | Content |
|---|---|
| `exceptional/paper/` | LaTeX source of the exceptional-set paper (REVTeX) and its eight figures as PNG |
| `paper/` | LaTeX source of the sieve-sum paper (amsart) and its three figures as PNG |
| `preprintErdosStrausExceptionalSet/`, `preprintErdosStrausSieveConstant/` | the two papers as PDF |
| `dist/` | for each paper the PDF, the LaTeX source zip and the arXiv tarball |
| `lean/` | Lean 4 and Mathlib, 59 theorems in `Axioms.lean`: the factorisation of φ(4ua), the local factors, the value 9/π², the Type II identity, the normal forms of the three corners, the identities and inequalities of the six-parameter family, the optimisation of the large sieve, c₀³ = 729/(256π²); the Farey bound and the parametrisation of the coincidences of the centre, the family of coincidences with the primes 531359 and 1442783; the range of Maynard's theorem, the divisor chain, the area of the strips and the share 40/49; the closed-form steps of the certified bound for 𝒥 and the inequality c₁ > 0.67916; for the deep strips the bound 5s < 3 behind the depth 1/10, the monotonicity of Maynard's conditions, the areas 3/98 and 88/3675, the share 24/25 and the inequality c₁ > 0.6931 |
| `checks/` | C, C++, Python and Julia programs: W(J) to J = 20000; the corners and the centre over the 332,398 primes p ≡ 3 (mod 4) below 10⁷; the six-parameter family against all coincidences with (p + 1)/4 ≤ 4000; the mean values to X = 2·10⁶; the two bounds of the large sieve to T = 10¹⁵; the level of the distinctness argument to 10⁶; the divisor chain in 128,000,000 cases; the integral 𝒥 by quadrature (Python) and by the midpoint rule (Julia); a certified upper bound for 𝒥 in directed rounding (C, `certify_J.c`, with its cell-by-cell output `certify_J.out`); for the deep strips the Monte Carlo simulation `deep_mc.c` with a larger-sample check (`deep_mc_check.sh`, `deep_mc_check.txt`), the certified lower bound for 𝒦 (C, `certify_K.c`, output `certify_K.out`) and its exact sum (Python, `sum_K.py`) |
| `engine/` | the three Anaxia runs (`anaxia_regions.json` for the corners, `anaxia_centre.json` for the strips, `anaxia_deep.json` for the deep strips), the scripts that reproduce them and the comparison of a fresh run with each record |
| `figures/` | the scripts that draw the figures |
| `scripts/run_all.sh` | builds the Lean project, checks the axioms, runs the programs and the engine reproduction |
| `scripts/build_paper.sh` | builds both papers and the release files |
| `MasterErdosStrausExceptionalSet_masterzipfinal2.zip`, `MasterErdosStrausSieveConstant_masterzipfinal2.zip` | every file of the repository with the files of each paper |

Every Lean theorem depends at most on the axioms propext, Classical.choice and Quot.sound, with no `sorry`. The analytic estimates of the proofs are not formalised.

To rerun everything: `bash scripts/run_all.sh` (Lean through elan; the engine step installs Anaxia from [creelie/eBLLM](https://github.com/creelie/eBLLM) at tag v3.0.5 in a virtual environment).

## Releases

Every release is archived on Zenodo. The concept DOI [10.5281/zenodo.23480292](https://doi.org/10.5281/zenodo.23480292) always resolves to the latest version.

| Version | Content | DOI |
|---|---|---|
| 4.0.1 | two corrections in the appendix: the error term of the hyperbola sum (Appendix A.2) and the bound for the small variables in Case (ii) of Appendix A.5; theorems and constants unchanged | [10.5281/zenodo.23770866](https://doi.org/10.5281/zenodo.23770866) |
| 4.0.0 | the deep strips through Maynard's second theorem: certified c₁ > 0.6931, Monte Carlo c₁ ≈ 0.7015 (Section 7, Appendix A.11, final inequality in Lean) | [10.5281/zenodo.23764189](https://doi.org/10.5281/zenodo.23764189) |
| 3.1.0 | the certified bound c₁ > 0.67916 (Appendix A.8, final inequality in Lean) | [10.5281/zenodo.23486200](https://doi.org/10.5281/zenodo.23486200) |
| 3.0.0 | Theorem A with the centre and Maynard's strips: c₁ = 0.6792, at least 0.66571 in closed form | [10.5281/zenodo.23483694](https://doi.org/10.5281/zenodo.23483694) |
| 2.0.0 | the exceptional-set paper with c₀ = 0.6607886973 from the three corners | [10.5281/zenodo.23480673](https://doi.org/10.5281/zenodo.23480673) |
| 1.0.0 | the sieve-sum paper | [10.5281/zenodo.23480293](https://doi.org/10.5281/zenodo.23480293) |

## License

MIT, see `LICENSE`.
