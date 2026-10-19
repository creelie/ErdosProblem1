**The exact constant of the two-parameter Type II sieve for the Erdős–Straus equation**, by Deep Bhattacharjee (corresponding author) and Ushashi Bhattacharya.

B. Dahan ([arXiv:2608.24035](https://arxiv.org/abs/2608.24035)) bounds the primes for which a depth-J search finds no Type II solution of 4/n = 1/x + 1/y + 1/z through the sum W(J) = Σ 1/φ(4ua) over coprime u, a ≤ J, and asked for W(J) ≥ (0.45 + o(1))(log J)².

**Theorem.** W(J) = κ((log J + γ + T)² + M) + O(J^{-1/2} log J) with κ = 9/(2π²) = 0.4559453…, T = 0.3131549610…, M = −1.0429405079….

The constant is three times the one proved before; in Dahan's framework the exponent 0.065 of his bound N exp(−c(log N)^{2/3}) becomes 0.0943. The Erdős–Straus conjecture itself remains open.

Lean 4 with Mathlib checks the algebraic steps (axioms propext, Classical.choice, Quot.sound only); C, C++, Python and Julia programs check W(J) to J = 20000, the inversion in exact arithmetic and the constants.

Files:
- `ErdosStrausSieveConstant_final2.pdf`: the paper
- `ErdosStrausSieveConstant_final2-tex.zip`: LaTeX source (`main.tex`, figures as PNG)
- `ErdosStrausSieveConstant_final2-arxiv.tar.gz`: the same source, ready for arXiv
