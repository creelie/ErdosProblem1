# AI provenance: An explicit exponent for the exceptional set of the Erdős–Straus equation

Authors: Deep Bhattacharjee (corresponding author, ORCID 0000-0003-0466-750X)
and Ushashi Bhattacharya (ORCID 0009-0002-2254-3914).

- Paper: https://github.com/creelie/ErdosProblem1/blob/main/preprintErdosStrausExceptionalSet/ErdosStrausExceptionalSet.pdf
- Repository: https://github.com/creelie/ErdosProblem1

## Result

Let E(N) be the set of n ≤ N for which 4/n = 1/x + 1/y + 1/z has no solution in
positive integers. Vaughan (1970) proved |E(N)| ≪ N exp(−c (log N)^{2/3}) for some
c > 0, with no value of c. For the primes in E(N), which satisfy n ≡ 1 (mod 24),
the sieve of Dahan (arXiv:2608.24035) gives c = 0.065, and c = 0.0943 with the
constant 9/(2π²) of its sieve sum (Liu, doi:10.5281/zenodo.22650516). The paper
proves the bound for the whole set E(N) with every exponent below
c₀ = 9/(2^{8/3} π^{2/3}) = 0.6607886973… (Theorem A), and proves that, with the
residue classes it uses, the large sieve gives no better exponent (Theorem B).

## Who did what

Deep Bhattacharjee chose the problem, directed the work, set the format and the
target journals, and ran the computations. The construction of the residue
classes (three corners of the triangle of factorisations, each with its normal
form) was selected among 128 candidate constructions by the engine Anaxia of the
BLLM project, version 3.0.5 (doi:10.5281/zenodo.23469389), which Deep
Bhattacharjee develops; the run is recorded in `engine/anaxia_regions.json` and
reproduced by `engine/anaxia_regions.py`. Claude (Anthropic), used through
Claude Code under his direction, produced the proofs (the distinctness of the
classes within and between the corners, the six-parameter family of
coincidences between two corners, the mean values, the tilted lower bound and
the matching upper bound), the Lean 4 formalisation, the checking programs, the
figures and the LaTeX text. No other AI system took part.

## How the result is checked

- Lean 4.34.1 with Mathlib (19 theorems for this paper): the three normal forms,
  the identities and the three inequalities of the six-parameter family, the
  relation φ(qn) = qφ(n) for q | n, the local factors, the optimisation of the
  large sieve, c₀³ = 729/(256π²), and the coincidence at the prime 531359 with
  its primality. All depend at most on propext, Classical.choice and
  Quot.sound, with no sorry or added axiom. The analytic estimates are not
  formalised.
- Programs in C, C++, Python and Julia: over the 332,398 primes p ≡ 3 (mod 4)
  below 10⁷ (37,093,810 Type II triples) the corners are disjoint, the normal
  forms agree, the classes in each corner are distinct and no two corners share
  a class; the six-parameter family is checked against all coincidences with
  (p + 1)/4 ≤ 4000; the mean values are computed to X = 2·10⁶; the two bounds of
  the large sieve are evaluated to T = 10¹⁵; the level of the distinctness
  argument is checked to 10⁶.
- A fresh run of the engine is compared with its record by `engine/compare.py`.

## What is not claimed

The paper does not prove the Erdős–Straus conjecture. Theorem B concerns the
residue classes used in the paper; other classes, or other methods, may give a
larger exponent. On 19 October 2026 the authors searched the literature for
explicit values of Vaughan's constant: Elsholtz–Tao (2013) and
Pomerance–Weingartner (Ramanujan J. 69 (2026), Paper No. 31) leave their
constants unspecified, and the explicit values 0.065 (Dahan) and 0.0943 (Liu,
from Dahan's sieve) concern the primes in E(N) only. No explicit value for the
whole set E(N), and none as large as c₀, was found.
