# Hexagon submission fields

Title: An explicit exponent for the exceptional set of the Erdős–Straus equation

Abstract: The Erdős–Straus conjecture says that four divided by any integer greater than one is a sum of three unit fractions. In 1970 Vaughan showed that the integers up to N for which this fails are few: their number is at most N times the exponential of minus a positive constant times the two-thirds power of log N, and his argument gives no value of the constant. We prove this bound with a constant that a computation with directed rounding shows to exceed 0.6931. The residue classes we sieve with include classes to moduli up to the three-fifths power of each prime, which two recent theorems of Maynard let us count.

Primary subject: math.NT (Number Theory)
Cross-list subjects: none
Tags: Erdős–Straus conjecture, Egyptian fractions, large sieve, primes in arithmetic progressions
MSC2020: 11D68 Rational numbers as sums of fractions; 11N36 Applications of sieve methods; 11N35 Sieves; 11N13 Primes in congruence classes

Contributor 1 (submitting): Deep Bhattacharjee; ORCID 0000-0003-0466-750X; affiliation: Formerly, Electro-Gravitational Space Propulsion Laboratory (EGSPL), Bhubaneswar, Odisha 751030, India; formerly, University of Calcutta (Heramba Chandra College), 23/49 Gariahat Road, Golpark, Kolkata 700029, West Bengal, India; d.bhattacharjee@erl-forschung.de (corresponding; further addresses itsdeep@live.com and egspl@proton.me, the last for the GitHub repository)
Contributor 2: Ushashi Bhattacharya; ORCID 0009-0002-2254-3914; affiliation: Formerly, National Taiwan University, Taipei 106319, Taiwan; rai.bhattacharya.in@gmail.com

Text of the contribution: A mix of human-written and AI-generated.
Human understanding: "All parts" only after both authors have read every proof; otherwise "Some parts, but not all".

Provenance statement: The constructions of the residue classes (the three corners, the strips of the centre with the rule that finds a divisor in Maynard's range, and the splitting of the moduli in the deep strips) were selected in three runs of the engine Anaxia of the BLLM project, developed by Deep Bhattacharjee: two with version 3.0.5 (doi:10.5281/zenodo.23469389) and one with version 3.1.0 (doi:10.5281/zenodo.23484160). The proofs, the text, the figures, the Lean 4 formalisation and the C, C++, Python, Julia and shell programs were carried out by Deep Bhattacharjee with AI assistance (Claude, Anthropic) in the Claude Code environment. The algebraic steps are checked in Lean 4 with Mathlib (axioms propext, Classical.choice and Quot.sound only); the analytic estimates are not formalised; the counts of classes are checked over all primes p ≡ 3 (mod 4) below 10⁷, the integral that gives the constant is computed by two independent programs, and two C programs with directed rounding prove the bounds 𝒥 < 0.000335 and 𝒦 > 0.00859, hence c₁ > 0.6931. The paper does not prove the Erdős–Straus conjecture. Both authors approved the final version.

Tools used: Anaxia (BLLM project) 3.0.5 and 3.1.0 (doi:10.5281/zenodo.23484160); Anthropic - Claude (Claude Code); Lean community - Lean 4 (v4.34.1) with Mathlib; C; C++; Python; Julia; shell.
arXiv identifier: none (not posted on arXiv).
Formalization: https://github.com/creelie/ErdosProblem1/tree/main/lean

Source files or ZIP: ErdosStrausExceptionalSet_final2-tex.zip (main.tex and figures/ as PNG; compiles with pdflatex, no BibTeX pass). Do not upload the master zip or the Lean files here.
AI provenance files: AI_PROVENANCE_ErdosStrausExceptionalSet.md
License: CC BY 4.0
Funding: none. Competing interests: none.
