# Journal submission

Target: *Integers: Electronic Journal of Combinatorial Number Theory* (open access, no fees, accepts short papers in elementary and analytic number theory).

Submission file: `dist/ErdosStrausSieveConstant_final2.pdf`; source `dist/ErdosStrausSieveConstant_final2-tex.zip` (amsart, `main.tex` with `figures/*.png`, no BibTeX pass).

Title: A power-saving asymptotic for the two-parameter Type II sieve sum of the Erdős–Straus equation
Authors: Deep Bhattacharjee (corresponding; itsdeep@live.com, d.bhattacharjee@erl-forschung.de; ORCID 0000-0003-0466-750X) and Ushashi Bhattacharya (rai.bhattacharya.in@gmail.com; ORCID 0009-0002-2254-3914)
MSC2020: 11D68 (primary); 11N36, 11N37, 11A25
Keywords: Erdős–Straus conjecture, Egyptian fractions, sieve methods, Euler totient function

Recommendation: do not submit this paper on its own. Its leading constant 9/(2π²) and linear term were found earlier by Z. Liu (Zenodo, 7 September 2026, doi:10.5281/zenodo.22650516), so what is new is the constant term and the power-saving error. That refinement is better used as a lemma of the exceptional-set paper (`JOURNAL_ErdosStrausExceptionalSet.md`), which cites it.

Cover letter points, if it is submitted:
- The paper determines W(J) = Σ 1/φ(4ua) over coprime u, a ≤ J up to O(J^{-1/2} log J), the question of Dahan, arXiv:2608.24035.
- The first two terms were found earlier by Liu; the paper adds the constant term κ((γ + T)² + M) and the power-saving error, with a Lean check of the algebraic steps. Corollary 1.2 records the consequences in Dahan's framework and states that they rest on his unrefereed arguments.
- The paper does not claim the Erdős–Straus conjecture.

Declarations (already in the paper): funding none; competing interests none; data and code at https://github.com/creelie/ErdosProblem1 and in the Zenodo archive (concept doi:10.5281/zenodo.23480292); AI assistance stated in the data availability paragraph.
Suggested referees: none proposed.
