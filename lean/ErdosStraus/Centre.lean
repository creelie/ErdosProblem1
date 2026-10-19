import Mathlib

/-!
# The centre of the triangle

The identities and inequalities behind the centre strips of the paper
"An explicit exponent for the exceptional set of the Erdős–Straus equation".
A Type II triple for a prime `p = 3 (mod 4)` is `(u, a, α)` with `u a α = P = (p + 1)/4`; it lies in
the centre when `u`, `a` and `α` are all below `√P`.

* `centre_farey`, `centre_B_form`: the classes of the centre are distinct, and a centre triple with
  `u < a` has its `M`-form below `p`.
* `centre_C_j`, `centre_C_D`, `centre_C_alpha`, `centre_C_j_le`, `centre_C_j_zero`: the
  parametrisation of a coincidence between the centre and the corner `C`.
* `centre_family`, `centre_family_coprime`, `centre_prime_two`: the one-parameter family of
  coincidences and its second member.
* `maynard_window`, `maynard_switch_tau`, `maynard_switch_gamma`: the ranges of Maynard's theorem.
* `chain_step`: the divisor chain.
* `strip_area`, `strip_share`: the area of the strips and the largest share they can add.
* `cert_r1_mono`, `cert_xi_lower`, `cert_dickman_tail`, `cert_constant`: the closed-form steps of the
  certified bound for the loss `𝒥`.
-/

namespace ErdosStraus

/-- Farey bound in the centre: if `u² < P` and `a'² < P`, then `u a' < P`. -/
theorem centre_farey (u a' P : ℕ) (hu : u * u < P) (ha : a' * a' < P) : u * a' < P := by
  by_contra h
  have h' : P ≤ u * a' := Nat.le_of_not_lt h
  have : P * P ≤ (u * a') * (u * a') := Nat.mul_le_mul h' h'
  nlinarith

/-- A centre triple with `u < a` has `4 u² α < 4 u a α = p + 1`, so its `M`-form lies below `p`. -/
theorem centre_B_form (u a α : ℕ) (hu : 0 < u) (hα : 0 < α) (h : u < a) :
    4 * (u ^ 2 * α) < 4 * (u * a * α) := by
  have : u * u * α < u * a * α := by
    apply Nat.mul_lt_mul_of_pos_right _ hα
    exact Nat.mul_lt_mul_of_pos_left h hu
  nlinarith

/-- Centre against `C`: from `u'' s = P = u a α` and `u'' a - u a'' = k p` with `p = 4P - 1` and
`s = a'' α''`, one gets `P (a - 4 k s) = s (u a'' - k)`. -/
theorem centre_C_j (u a α u'' a'' α'' k : ℤ) (hP : u'' * (a'' * α'') = u * a * α)
    (hk : u'' * a - u * a'' = k * (4 * (u * a * α) - 1)) :
    (u * a * α) * (a - 4 * k * (a'' * α'')) = (a'' * α'') * (u * a'' - k) := by
  linear_combination (a'' * α'') * hk - a * hP

/-- Centre against `C`: with `j = a - 4 k s` and `D = s a'' - a α j`, one has `u D = s k`. -/
theorem centre_C_D (u a α s a'' k j : ℤ) (h : (u * a * α) * j = s * (u * a'' - k)) :
    u * (s * a'' - a * α * j) = s * k := by
  linear_combination -h

/-- Centre against `C`: with `a = j + 4 k s`, `α j (j + 4 k s) = s a'' - D`. -/
theorem centre_C_alpha (a α s a'' k j D : ℤ) (ha : a = j + 4 * k * s)
    (hD : D = s * a'' - a * α * j) : α * j * (j + 4 * k * s) = s * a'' - D := by
  subst ha hD
  ring

/-- Centre against `C`: `4 k s ≤ a`, so `j = a - 4 k s ≥ 0`. The hypotheses are
`u'' s = P`, `u'' a - u a'' = k p`, `u a'' > 0`, `a < p`, `s > 0` and `P > 0`, with `p = 4P - 1`. -/
theorem centre_C_j_le (u a u'' a'' s k P : ℤ) (hP : u'' * s = P)
    (hk : u'' * a - u * a'' = k * (4 * P - 1)) (hpos : 0 < u * a'') (hap : a < 4 * P - 1)
    (hs : 0 < s) (hP0 : 0 < P) : 4 * k * s ≤ a := by
  have h0 : k * (4 * P - 1) < u'' * a := by linarith
  have h1 : k * s * (4 * P - 1) < P * a := by
    have := mul_lt_mul_of_pos_right h0 hs
    have e : u'' * a * s = P * a := by rw [← hP]; ring
    linarith [this, e, show k * (4 * P - 1) * s = k * s * (4 * P - 1) by ring]
  have h4 : P * a < P * (4 * P - 1) := mul_lt_mul_of_pos_left hap hP0
  have hks : k * s < P := lt_of_mul_lt_mul_right (by linarith : k * s * (4 * P - 1) < P * (4 * P - 1))
    (by linarith : (0 : ℤ) ≤ 4 * P - 1)
  by_contra h
  have h2 : a + 1 ≤ 4 * k * s := by omega
  have h3 : (a + 1) * P ≤ 4 * k * s * P := mul_le_mul_of_nonneg_right h2 hP0.le
  nlinarith

/-- Centre against `C`: if `j = 0`, then `u` divides `a`, which is impossible in the centre. -/
theorem centre_C_j_zero (u a α a'' s k : ℤ) (hs : s ≠ 0) (ha : a = 4 * k * s)
    (h : (u * a * α) * (a - 4 * k * s) = s * (u * a'' - k)) : u ∣ a := by
  have h0 : s * (u * a'' - k) = 0 := by rw [← h, ha]; ring
  have hk : k = u * a'' := by
    rcases mul_eq_zero.mp h0 with h | h
    · exact absurd h hs
    · linarith
  exact ⟨4 * a'' * s, by rw [ha, hk]; ring⟩

/-- The family of coincidences: for `s = 16 α + 1`, `a = 8 (8 α + 1)` and `P = s a α`, the triples
`(s, a, α)` of the centre and `(a α, s, 1)` of `C` satisfy `(a α) a - s · s = 4P - 1`, so they have
the same class modulo `p = 4P - 1`. -/
theorem centre_family (α : ℤ) :
    (8 * (8 * α + 1) * α) * (8 * (8 * α + 1)) - (16 * α + 1) * (16 * α + 1)
      = 4 * ((16 * α + 1) * (8 * (8 * α + 1)) * α) - 1 := by
  ring

/-- In the family, `2 (8 α + 1) - (16 α + 1) = 1`, so `s = 16 α + 1` is prime to `a = 8 (8 α + 1)`
(and `s` is odd). -/
theorem centre_family_coprime (α : ℤ) : 2 * (8 * α + 1) - (16 * α + 1) = 1 := by ring

/-- The second member of the family below `10^7`: `α = 7`, `p = 1442783`. -/
theorem centre_prime_two : Nat.Prime 1442783 ∧ 4 * (113 * 456 * 7) - 1 = 1442783 := by
  constructor
  · norm_num
  · norm_num

/-- Maynard's window is not empty: for `0 < δ < 1/42`, `2δ` lies below `1/10 - 7δ/5` and
`1/2 - 19δ`. -/
theorem maynard_window (δ : ℝ) (h0 : 0 < δ) (h1 : δ < 1 / 42) :
    2 * δ < 1 / 10 - 7 / 5 * δ ∧ 2 * δ < 1 / 2 - 19 * δ := by
  constructor <;> linarith

/-- The upper end of the window is `1/10 - 7δ/5` exactly when `δ ≤ 1/44`. -/
theorem maynard_switch_tau (δ : ℝ) : 1 / 10 - 7 / 5 * δ ≤ 1 / 2 - 19 * δ ↔ δ ≤ 1 / 44 := by
  constructor <;> intro h <;> linarith

/-- For `δ ≤ 1/44` the width `β - α = 1/10 - 17δ/5` of the window is at least `α = 2δ` exactly
when `δ ≤ 1/54`. -/
theorem maynard_switch_gamma (δ : ℝ) : 2 * δ ≤ (1 / 10 - 7 / 5 * δ) - 2 * δ ↔ δ ≤ 1 / 54 := by
  constructor <;> intro h <;> linarith

/-- The divisor chain: if `d₀ ≤ A < d₀ q₁ ⋯ q_m` with `d₀ > 0` and every `qᵢ ≥ 1` satisfying
`qᵢ A < B`, then some partial product `d₀ q₁ ⋯ qᵢ` lies strictly between `A` and `B`. -/
theorem chain_step (A B : ℝ) :
    ∀ (qs : List ℝ) (d₀ : ℝ), 0 < d₀ → (∀ q ∈ qs, 1 ≤ q ∧ q * A < B) → d₀ ≤ A →
      A < d₀ * qs.prod → ∃ i ≤ qs.length, A < d₀ * (qs.take i).prod ∧ d₀ * (qs.take i).prod < B
  | [], d₀, _, _, h0, h1 => by simp at h1; linarith
  | q :: qs, d₀, hd, hq, h0, h1 => by
    have hq1 := hq q (List.mem_cons_self ..)
    by_cases hc : A < d₀ * q
    · refine ⟨1, by simp, ?_, ?_⟩
      · simpa using hc
      · have : d₀ * q ≤ A * q := by nlinarith
        simpa using lt_of_le_of_lt this (by linarith [hq1.2])
    · have hc : d₀ * q ≤ A := not_lt.mp hc
      have hd' : 0 < d₀ * q := by nlinarith
      have hrest : ∀ r ∈ qs, 1 ≤ r ∧ r * A < B := fun r hr => hq r (List.mem_cons_of_mem _ hr)
      have h1' : A < d₀ * q * qs.prod := by simpa [mul_assoc] using h1
      obtain ⟨i, hi, hlo, hhi⟩ := chain_step A B qs (d₀ * q) hd' hrest hc h1'
      refine ⟨i + 1, by simp; omega, ?_, ?_⟩
      · simpa [List.take_succ_cons, mul_assoc] using hlo
      · simpa [List.take_succ_cons, mul_assoc] using hhi

/-- The area of one strip: `∫₀^{1/42} (1/2 - 3δ) dδ = 1/84 - 3/3528`. -/
theorem strip_area : ∫ δ in (0 : ℝ)..(1 / 42), (1 / 2 - 3 * δ) = 1 / 84 - 3 / 3528 := by
  rw [intervalIntegral.integral_sub intervalIntegrable_const
      ((intervalIntegral.intervalIntegrable_id).const_mul 3), intervalIntegral.integral_const,
    intervalIntegral.integral_const_mul, integral_id]
  norm_num

/-- With the three strips at full density the share of the triangle would be
`3/4 + 6 (1/84 - 3/3528) = 40/49`. -/
theorem strip_share : (3 / 4 : ℝ) + 6 * (1 / 84 - 3 / 3528) = 40 / 49 := by norm_num

/-- The bound of `(R1)` used at the two ends of the certified computation:
`1/2 + δ - 2 θ(δ)` is nondecreasing in `δ`. -/
theorem cert_r1_mono (x y : ℝ) (hxy : x ≤ y) :
    1 / 2 + x - 2 * min (1 / 10 - 7 / 5 * x) (1 / 2 - 19 * x)
      ≤ 1 / 2 + y - 2 * min (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y) := by
  have h : min (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y) ≤ min (1 / 10 - 7 / 5 * x) (1 / 2 - 19 * x) :=
    min_le_min (by linarith) (by linarith)
  linarith

/-- On a cell `a ≤ δ ≤ b` the number `ξ(δ) = min(θ(δ) - 2δ, 2δ)` is at least the smaller of its
values at the ends, because it is a minimum of linear functions. -/
theorem cert_xi_lower (a b x : ℝ) (hax : a ≤ x) (hxb : x ≤ b) :
    min (min (min (1 / 10 - 7 / 5 * a) (1 / 2 - 19 * a) - 2 * a) (2 * a))
        (min (min (1 / 10 - 7 / 5 * b) (1 / 2 - 19 * b) - 2 * b) (2 * b))
      ≤ min (min (1 / 10 - 7 / 5 * x) (1 / 2 - 19 * x) - 2 * x) (2 * x) := by
  have key : ∀ c m : ℝ, min (c + m * a) (c + m * b) ≤ c + m * x := by
    intro c m
    rcases le_total 0 m with h | h
    · exact min_le_of_left_le (by nlinarith)
    · exact min_le_of_right_le (by nlinarith)
  have hA : ∀ y : ℝ, min (min (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y) - 2 * y) (2 * y)
      ≤ 1 / 10 + (-17 / 5) * y := by
    intro y
    have h1 := min_le_left (min (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y) - 2 * y) (2 * y)
    have h2 := min_le_left (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y)
    linarith
  have hB : ∀ y : ℝ, min (min (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y) - 2 * y) (2 * y)
      ≤ 1 / 2 + (-21) * y := by
    intro y
    have h1 := min_le_left (min (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y) - 2 * y) (2 * y)
    have h2 := min_le_right (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y)
    linarith
  have hC : ∀ y : ℝ, min (min (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y) - 2 * y) (2 * y)
      ≤ 0 + 2 * y := by
    intro y
    have h1 := min_le_right (min (1 / 10 - 7 / 5 * y) (1 / 2 - 19 * y) - 2 * y) (2 * y)
    linarith
  rw [← min_sub_sub_right (1 / 10 - 7 / 5 * x) (1 / 2 - 19 * x) (2 * x)]
  refine le_min (le_min ?_ ?_) ?_
  · calc _ ≤ min (1 / 10 + (-17 / 5) * a) (1 / 10 + (-17 / 5) * b) := min_le_min (hA a) (hA b)
      _ ≤ 1 / 10 + (-17 / 5) * x := key _ _
      _ = 1 / 10 - 7 / 5 * x - 2 * x := by ring
  · calc _ ≤ min (1 / 2 + (-21) * a) (1 / 2 + (-21) * b) := min_le_min (hB a) (hB b)
      _ ≤ 1 / 2 + (-21) * x := key _ _
      _ = 1 / 2 - 19 * x - 2 * x := by ring
  · calc _ ≤ min (0 + 2 * a) (0 + 2 * b) := min_le_min (hC a) (hC b)
      _ ≤ 0 + 2 * x := key _ _
      _ = 2 * x := by ring

/-- The tail of Dickman's function: from `ρ(v) ≤ ρ(v - 1)/13` for `v ≥ 13`, the integral of `ρ` over
`[12, ∞)` is at most `ρ(12)` times `∑ 13⁻ᵏ ≤ 13/12`. -/
theorem cert_dickman_tail (n : ℕ) : ∑ k ∈ Finset.range n, (1 / 13 : ℝ) ^ k ≤ 13 / 12 := by
  rw [geom_sum_eq (by norm_num)]
  have h : (0 : ℝ) ≤ (1 / 13) ^ n := by positivity
  rw [div_le_iff_of_neg (by norm_num)]
  nlinarith

/-- The certified constant: if `𝒥 ≤ 0.000335` then `c₁ > 0.67916`, where
`c₁³ = (27/8) (9/(8π²)) (40/49 - 6𝒥)`; it uses `π < 3.141593`. -/
theorem cert_constant (J : ℝ) (hJ : J ≤ 0.000335) :
    (0.67916 : ℝ) ^ 3 < 27 / 8 * (9 / (8 * Real.pi ^ 2)) * (40 / 49 - 6 * J) := by
  have hpi0 : 0 < Real.pi := Real.pi_pos
  have hpi : Real.pi < 3.141593 := Real.pi_lt_d6
  have hsq : Real.pi ^ 2 < 3.141593 ^ 2 := by nlinarith
  have hfrac : 9 / (8 * (3.141593 : ℝ) ^ 2) < 9 / (8 * Real.pi ^ 2) := by
    apply div_lt_div_of_pos_left (by norm_num) (by positivity) (by nlinarith)
  have hJ' : (40 / 49 - 6 * 0.000335 : ℝ) ≤ 40 / 49 - 6 * J := by linarith
  have hpos : (0 : ℝ) < 40 / 49 - 6 * 0.000335 := by norm_num
  calc (0.67916 : ℝ) ^ 3 < 27 / 8 * (9 / (8 * (3.141593 : ℝ) ^ 2)) * (40 / 49 - 6 * 0.000335) := by
        norm_num
    _ ≤ 27 / 8 * (9 / (8 * Real.pi ^ 2)) * (40 / 49 - 6 * J) := by
        apply mul_le_mul (by linarith) hJ' hpos.le (by positivity)

end ErdosStraus
