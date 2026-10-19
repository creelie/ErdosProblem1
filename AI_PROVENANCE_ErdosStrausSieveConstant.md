# AI provenance: A power-saving asymptotic for the two-parameter Type II sieve sum of the Erdős–Straus equation

Authors: Deep Bhattacharjee (corresponding author, ORCID 0000-0003-0466-750X)
and Ushashi Bhattacharya (ORCID 0009-0002-2254-3914).

- Paper: https://github.com/creelie/ErdosProblem1/blob/main/preprintErdosStrausSieveConstant/ErdosStrausSieveConstant.pdf
- Repository: https://github.com/creelie/ErdosProblem1

## Result

B. Dahan (arXiv:2608.24035) bounds the primes n ≡ 1 (mod 24) at which a depth-J
search finds no Type II solution of 4/n = 1/x + 1/y + 1/z through the sum
W(J) of 1/φ(4ua) over coprime u, a ≤ J. He proved
W(J) ≥ (3/(2π²) + o(1))(log J)² and asked for (0.45 + o(1))(log J)². The paper
proves W(J) = (9/(2π²))((log J + γ + T)² + M) + O(J^{-1/2} log J) with explicit
constants T and M, and records the resulting constants in Dahan's bounds.

The first two terms, (9/(2π²))(log J)² and 2κ(γ + T) log J with an error O(1),
were found earlier by Z. Liu, *The Exact Reciprocal-Totient Constant in the
Two-Parameter Erdős–Straus Sieve* (Zenodo, 7 September 2026,
doi:10.5281/zenodo.22650516), who also derived the exponent 0.0943 for Dahan's
bound. The paper credits Liu; what it adds is the constant term
κ((γ + T)² + M) = −0.1140687011… and the error O(J^{-1/2} log J), with the
Lean checks.

## Who did what

Deep Bhattacharjee chose the problem, directed the work, set the format and the
target journal, and ran the computations. Claude (Anthropic), used through
Claude Code under his direction, produced the proof (the multiplicativity of
ua/φ(4ua), the two-variable inversion and the Euler product), the Lean 4
formalisation, the checking programs, the figures and the LaTeX text. No other
AI system took part.

## How the result is checked

- Lean 4.34.1 with Mathlib: the factorisation of φ(4ua), the local values of the
  inverted function, the local factors and their derivatives, the value
  Z(0,0) = 9/π², and the Type II identity, all depending only on propext,
  Classical.choice and Quot.sound, with no sorry or added axiom. The analytic
  estimates (the error terms) are not formalised.
- Programs in C, C++, Python and Julia compute W(J) exactly for J ≤ 20000 and
  compare it with the formula, check the inversion in exact rational arithmetic,
  check the local factors exactly, and compute T and M.

## What is not claimed

The paper does not prove the Erdős–Straus conjecture. The consequences for the
exceptional set rest on Dahan's arguments, which are unrefereed and are used as
stated. On 19 October 2026 the authors found version 1 of Dahan's preprint
(25 August 2026) to be its only version, and found Liu's preprint, which
answers Dahan's first question with the leading constant and the linear term.
The paper claims neither as new.
