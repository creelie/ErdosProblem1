# The exact constant of the two-parameter Type II sieve for the Erdős–Straus equation

Deep Bhattacharjee (corresponding author, [ORCID 0000-0003-0466-750X](https://orcid.org/0000-0003-0466-750X)) and Ushashi Bhattacharya ([ORCID 0009-0002-2254-3914](https://orcid.org/0009-0002-2254-3914)).

The paper is [`preprintErdosStrausSieveConstant/ErdosStrausSieveConstant.pdf`](preprintErdosStrausSieveConstant/ErdosStrausSieveConstant.pdf).

## Result

The Erdős–Straus conjecture says that 4/n = 1/x + 1/y + 1/z has a solution in positive integers for every n ≥ 2. It is open. B. Dahan ([arXiv:2608.24035](https://arxiv.org/abs/2608.24035)) bounds the primes n ≡ 1 (mod 24) whose Type II solutions are not found by a search of depth J through a sieve governed by

W(J) = Σ_{u, a ≤ J, gcd(u, a) = 1} 1/φ(4ua),

proved W(J) ≥ (3/(2π²) + o(1)) (log J)², and asked, as the first open question of the paper, for W(J) ≥ (0.45 + o(1)) (log J)².

**Theorem.** W(J) = κ((log J + γ + T)² + M) + O(J^{-1/2} log J), with κ = 9/(2π²) = 0.4559453…, T = 0.3131549610… and M = −1.0429405079…, both given by explicit sums over primes.

So W(J) ~ (9/(2π²)) (log J)², three times the earlier constant. Substituted into Dahan's arguments this raises the exponent of his bound N exp(−c (log N)^{2/3}) for the exceptional primes from c = 0.065 to c = 0.0943, and triples the constants of his other two bounds. The paper does not prove the Erdős–Straus conjecture.

## Contents

| Path | Content |
|---|---|
| `paper/` | LaTeX source (amsart) and the three figures as PNG |
| `preprintErdosStrausSieveConstant/` | the paper as a PDF |
| `dist/` | the PDF, the LaTeX source zip and the arXiv tarball |
| `lean/` | Lean 4 and Mathlib: the factorisation of φ(4ua), the local factors and derivatives, the value 9/π², the Type II identity |
| `checks/` | C, C++, Python and Julia programs: W(J) against the formula up to J = 20000, the two-variable inversion in exact arithmetic, the local factors, the constants T and M |
| `figures/` | the script that draws the figures |
| `scripts/run_all.sh` | builds the Lean project, checks the axioms and runs the programs |
| `scripts/build_paper.sh` | builds the paper and the release files |
| `MasterErdosStrausSieveConstant_masterzipfinal2.zip` | every file of the repository |

Every Lean theorem depends only on the axioms propext, Classical.choice and Quot.sound, with no `sorry`. The analytic estimates of the proof are not formalised.

## License

MIT, see `LICENSE`.
