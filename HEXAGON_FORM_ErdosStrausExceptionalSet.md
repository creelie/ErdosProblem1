# Hexagon submission fields

Title: An explicit exponent for the exceptional set of the Erdős–Straus equation

Abstract: The Erdős–Straus conjecture says that four divided by any integer greater than one is a sum of three unit fractions. In 1970 Vaughan showed that the integers up to N for which this fails are few: their number is at most N times the exponential of minus a positive constant times the two-thirds power of log N. His argument gives no value of the constant, and for the primes among these integers a recent sieve gives the value 0.0943. We prove the bound for all of them with the constant 0.6607, which has a closed form in π, and we show that it is the best constant the large sieve can give with the residue classes we use. These classes come from three regions of the factorisations attached to each prime, and we prove that classes from different regions agree only at a sparse set of primes.

Primary subject: math.NT (Number Theory)
Cross-list subjects: none
Tags: Erdős–Straus conjecture, large sieve, exceptional set, Egyptian fractions
MSC2020: 11D68 Rational numbers as sums of fractions; 11N36 Applications of sieve methods; 11N35 Sieves; 11N13 Primes in congruence classes

Contributor 1 (submitting): Deep Bhattacharjee; ORCID 0000-0003-0466-750X; affiliation: Formerly, Electro-Gravitational Space Propulsion Laboratory (EGSPL), Bhubaneswar, Odisha 751030, India; formerly, University of Calcutta, Senate House, 87/1 College Street, Kolkata 700073, West Bengal, India; itsdeep@live.com (second address d.bhattacharjee@erl-forschung.de)
Contributor 2: Ushashi Bhattacharya; ORCID 0009-0002-2254-3914; affiliation: Formerly, National Taiwan University, Taipei 106319, Taiwan; rai.bhattacharya.in@gmail.com

Text of the contribution: A mix of human-written and AI-generated.
Human understanding: "All parts" only after both authors have read every proof; otherwise "Some parts, but not all".

Provenance statement: The construction of the residue classes was selected with the engine Anaxia of the BLLM project, version 3.0.5 (doi:10.5281/zenodo.23469389), developed by Deep Bhattacharjee. The proofs, the text, the figures, the Lean 4 formalisation and the C, C++, Python, Julia and shell programs were carried out by Deep Bhattacharjee with AI assistance (Claude, Anthropic) in the Claude Code environment. The algebraic steps are checked in Lean 4 with Mathlib (axioms propext, Classical.choice and Quot.sound only); the analytic estimates are not formalised; the counts of classes are checked over all primes p ≡ 3 (mod 4) below 10⁷. The paper does not prove the Erdős–Straus conjecture. Both authors approved the final version.

Tools used: Anaxia (BLLM project) 3.0.5; Anthropic - Claude (Claude Code); Lean community - Lean 4 (v4.34.1) with Mathlib; C; C++; Python; Julia; shell.
arXiv identifier: none (not posted on arXiv).
Formalization: https://github.com/creelie/ErdosProblem1/tree/main/lean

Source files or ZIP: ErdosStrausExceptionalSet_final2-tex.zip (main.tex and figures/ as PNG; compiles with pdflatex, no BibTeX pass). Do not upload the master zip or the Lean files here.
AI provenance files: AI_PROVENANCE_ErdosStrausExceptionalSet.md
License: CC BY 4.0
Funding: none. Competing interests: none.
