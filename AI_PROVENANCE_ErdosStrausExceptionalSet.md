# AI provenance: An explicit exponent for the exceptional set of the Erdős–Straus equation

Authors: Deep Bhattacharjee (corresponding author, ORCID 0000-0003-0466-750X)
and Ushashi Bhattacharya (ORCID 0009-0002-2254-3914).

- Paper: https://github.com/creelie/ErdosProblem1/blob/main/preprintErdosStrausExceptionalSet/ErdosStrausExceptionalSet.pdf
- Repository: https://github.com/creelie/ErdosProblem1
- Archive of release v4.0.0: https://doi.org/10.5281/zenodo.23764189

## Result

Let E(N) be the set of n ≤ N for which 4/n = 1/x + 1/y + 1/z has no solution in
positive integers. Vaughan (1970) proved |E(N)| ≪ N exp(−c (log N)^{2/3}) for some
c > 0, with no value of c. For the primes in E(N), which satisfy n ≡ 1 (mod 24),
the sieve of Dahan (arXiv:2608.24035) gives c = 0.065, and c = 0.0943 with the
constant 9/(2π²) of its sieve sum (Liu, doi:10.5281/zenodo.22650516). The paper
proves the bound for the whole set E(N) with every exponent below
c₁ = (3/2)((9/(8π²))(40/49 − 6𝒥 + 6𝒦))^{1/3}, where 𝒥 is an explicit double
integral with 0 ≤ 𝒥 < 0.000335 and 𝒦 is the measure of the factorable pairs of
the deep strips, with 𝒦 > 0.00859, so that c₁ > 0.6931; in closed form
𝒥 ≤ 533/64680 and 𝒦 ≥ 0, so that c₁ ≥ (286983/(98560π²))^{1/3} = 0.6657102497…;
a Monte Carlo simulation gives 𝒦 ≈ 0.0138 and c₁ ≈ 0.7015 (Theorem A). It also proves that the residue classes of
the three corners of the triangle of factorisations give no better exponent
than c₀ = 9/(2^{8/3} π^{2/3}) = 0.6607886973… (Theorem B); the gain from c₀ to
c₁ comes from classes of the centre of the triangle, counted in strips of depth
1/42 with Maynard's theorem on primes in progressions to large moduli
(Mem. AMS 306 (2025), no. 1542) and, behind them to the depth 1/10, with his
second theorem on well-factorable moduli up to p^{3/5} (Mem. AMS 306 (2025),
no. 1543), for the pairs whose modulus splits into three factors of suitable
sizes.

## Who did what

Deep Bhattacharjee chose the problem, directed the work, set the format and the
target journals, and ran the computations. The constructions of the residue
classes were selected by the engine Anaxia of the BLLM project, which Deep
Bhattacharjee develops, in three runs, the first two with version 3.0.5
(doi:10.5281/zenodo.23469389):
the three corners of the triangle of factorisations with their normal forms,
among 128 candidate constructions (`engine/anaxia_regions.json`, reproduced by
`engine/anaxia_regions.py`), and the strips of the centre with the rule that
finds a divisor in Maynard's range, among 37 (`engine/anaxia_centre.json`,
reproduced by `engine/anaxia_centre.py`); and, with version 3.1.0
(doi:10.5281/zenodo.23484160), the splitting of the moduli in the deep strips,
among 15 (`engine/anaxia_deep.json`, reproduced by `engine/anaxia_deep.py`).
The later versions, up to 3.4.0, add no construction technique that applies to this problem.
Claude (Anthropic), used through
Claude Code under his direction, produced the proofs (the distinctness of the
classes within and between the corners and the centre, the six-parameter family
of coincidences between two corners, the family of coincidences between the
centre and a corner, the divisor chain, the density of the pairs that fail the
rule, the form of Maynard's second theorem with sizes that depend on the long
variable, the factorable profiles and the bound at the depth 1/10, the mean
values, the tilted lower bound and the matching upper bound),
the Lean 4 formalisation, the checking programs, the figures and the LaTeX
text. A review of version 4.0.0 prepared with OpenAI Codex pointed out
two corrections in the appendix, the intermediate error term of the hyperbola
sum (Appendix A.2) and the bound for the small variables in Case (ii) of the
proof of the adapted well-factorable estimate (Appendix A.5); they are made in
version 4.0.1 and change no theorem or constant. No other AI system took part.

## How the result is checked

- Lean 4.34.1 with Mathlib (59 theorems): the three normal forms, the
  identities and the three inequalities of the six-parameter family, the
  relation φ(qn) = qφ(n) for q | n, the local factors, the optimisation of the
  large sieve, c₀³ = 729/(256π²); the Farey bound of the centre, the
  parametrisation of its coincidences with the corner C and the family of such
  coincidences with the primes 531359 and 1442783; the range of Maynard's
  theorem, the divisor chain, the area of the strips and the share 40/49; the
  bound 5s < 3 behind the depth 1/10, the monotonicity of Maynard's
  conditions, the areas 3/98 and 88/3675, the share 24/25 and the final
  inequalities c₁ > 0.67916 and c₁ > 0.6931. All
  depend at most on propext, Classical.choice and Quot.sound, with no sorry or
  added axiom. The analytic estimates are not formalised.
- Programs in C, C++, Python and Julia: over the 332,398 primes p ≡ 3 (mod 4)
  below 10⁷ (37,093,810 Type II triples) the corners are disjoint, the normal
  forms agree, the classes in each corner are distinct and no two corners share
  a class; the 6,943,856 triples of the centre have distinct classes, and their
  only coincidences with the corners are the four pairs of the explicit family;
  the divisor chain holds in 128,000,000 cases; the six-parameter family is
  checked against all coincidences with (p + 1)/4 ≤ 4000; the mean values are
  computed to X = 2·10⁶; the two bounds of the large sieve are evaluated to
  T = 10¹⁵; the integral 𝒥 is computed by Gauss–Legendre quadrature (Python)
  and independently by the midpoint rule (Julia), which agree on c₁ to four
  decimals, and checked by Monte Carlo at four points. The value 0.6792 is a
  numerical evaluation. The bound 𝒥 < 0.000335, hence c₁ > 0.67916, is proved
  by the C program `checks/c/certify_J.c`, which bounds the integral cell by
  cell with every rounding directed upward (the proposition on the certified bound
  and its appendix in the paper; the closed-form steps and the final
  inequality are in Lean). The bound 𝒦 > 0.00859, hence c₁ > 0.6931, is proved
  by the C program `checks/c/certify_K.c`, which bounds the density of the
  factorable pairs from below on 47 cells of depth with every rounding directed
  against the bound, and `checks/python/sum_K.py`, which adds the cells in
  exact rational arithmetic. The value 𝒦 ≈ 0.0138 is a Monte Carlo estimate
  (`checks/c/deep_mc.c`, with a 20-fold larger sample at four depths in
  `checks/c/deep_mc_check.txt`). The bound 0.66571 is proved in closed form.
- A fresh run of each engine run is compared with its record by
  `engine/compare.py`.

## What is not claimed

The paper does not prove the Erdős–Straus conjecture. Theorem B concerns the
residue classes used in the paper; other classes, or other methods, may give a
larger exponent. On 19 October 2026 the authors searched the literature for
explicit values of Vaughan's constant: Elsholtz–Tao (2013) and
Pomerance–Weingartner (Ramanujan J. 69 (2026), Paper No. 31) leave their
constants unspecified, and the explicit values 0.065 (Dahan) and 0.0943 (Liu,
from Dahan's sieve) concern the primes in E(N) only. No explicit value for the
whole set E(N), and none as large as c₀, was found. The certified value 0.6931
is not the limit of the method: the simulation gives 0.7015, every pair of the
strips and deep strips would give 0.7175 (no splitting exists beyond the depth
1/10, the level p^{3/5} of Maynard's second theorem), and the whole centre,
which needs the Elliott–Halberstam conjecture, would give 0.7273.
