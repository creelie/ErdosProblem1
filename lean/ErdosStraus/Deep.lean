import Mathlib

/-!
# The deep strips

The linear arithmetic behind the deep strips: the conditions of Maynard's second theorem leave no
admissible triple of sum `1/2 + δ` once `δ ≥ 1/10`, and the areas and shares of the parts of the centre
that this allows, and the final inequality of the certified bound with the deep strips.
-/

namespace ErdosStraus

/-- The four conditions of Maynard's second theorem with margin `0` (the first three are used).  If a
triple `(d₁, d₂, d₃)` of sum `s` meets them for some `n`, then `5 s < 3`. -/
theorem deep_sum_bound (s n d₁ d₂ d₃ : ℝ) (hs : d₁ + d₂ + d₃ = s) (h1 : d₁ < n)
    (h2 : 2 * n + d₂ + 2 * d₃ < 1) (h3 : 2 * n + d₁ + 4 * d₂ + 3 * d₃ < 2) : 5 * s < 3 := by
  linarith

/-- Hence no triple of sum `1/2 + δ` is admissible once `δ ≥ 1/10`. -/
theorem deep_tenth (δ n d₁ d₂ d₃ : ℝ) (hδ : 1 / 10 ≤ δ) (hs : d₁ + d₂ + d₃ = 1 / 2 + δ)
    (h1 : d₁ < n) (h2 : 2 * n + d₂ + 2 * d₃ < 1) (h3 : 2 * n + d₁ + 4 * d₂ + 3 * d₃ < 2) : False := by
  have := deep_sum_bound (1 / 2 + δ) n d₁ d₂ d₃ hs h1 h2 h3
  linarith

/-- The conditions of Maynard's second theorem become weaker when the sizes decrease. -/
theorem deep_adm_mono (n μ d₁ d₂ d₃ e₁ e₂ e₃ : ℝ) (h₁ : e₁ ≤ d₁) (h₂ : e₂ ≤ d₂) (h₃ : e₃ ≤ d₃)
    (a1 : d₁ < n - μ) (a2 : 2 * n + d₂ + 2 * d₃ < 1 - μ)
    (a3 : 2 * n + d₁ + 4 * d₂ + 3 * d₃ < 2 - μ) (a4 : n + d₁ + 5 * d₂ + 2 * d₃ < 2 - μ) :
    e₁ < n - μ ∧ 2 * n + e₂ + 2 * e₃ < 1 - μ ∧ 2 * n + e₁ + 4 * e₂ + 3 * e₃ < 2 - μ ∧
      n + e₁ + 5 * e₂ + 2 * e₃ < 2 - μ := by
  refine ⟨by linarith, by linarith, by linarith, by linarith⟩

/-- The area of the deep strips: `∫_{1/42}^{1/6} (1/2 - 3δ) dδ = 3/98`. -/
theorem deep_area : ∫ δ in (1 / 42 : ℝ)..(1 / 6), (1 / 2 - 3 * δ) = 3 / 98 := by
  rw [intervalIntegral.integral_sub intervalIntegrable_const
      ((intervalIntegral.intervalIntegrable_id).const_mul 3), intervalIntegral.integral_const,
    intervalIntegral.integral_const_mul, integral_id]
  norm_num

/-- The part of the deep strips below the depth `1/10`: `∫_{1/42}^{1/10} (1/2 - 3δ) dδ = 88/3675`. -/
theorem deep_area_tenth : ∫ δ in (1 / 42 : ℝ)..(1 / 10), (1 / 2 - 3 * δ) = 88 / 3675 := by
  rw [intervalIntegral.integral_sub intervalIntegrable_const
      ((intervalIntegral.intervalIntegrable_id).const_mul 3), intervalIntegral.integral_const,
    intervalIntegral.integral_const_mul, integral_id]
  norm_num

/-- With every pair of the strips and of the deep strips below the depth `1/10` counted, the share
of the triangle would be `40/49 + 6 · 88/3675 = 24/25`. -/
theorem deep_share_tenth : (40 / 49 : ℝ) + 6 * (88 / 3675) = 24 / 25 := by norm_num

/-- The certified constant with the deep strips: if `𝒥 ≤ 0.000335` and `𝒦 ≥ 0.00859` then
`c₁ > 0.6931`, where `c₁³ = (27/8) (9/(8π²)) (40/49 - 6𝒥 + 6𝒦)`; it uses `π < 3.141593`. -/
theorem deep_cert_constant (J K : ℝ) (hJ : J ≤ 0.000335) (hK : 0.00859 ≤ K) :
    (0.6931 : ℝ) ^ 3 < 27 / 8 * (9 / (8 * Real.pi ^ 2)) * (40 / 49 - 6 * J + 6 * K) := by
  have hpi0 : 0 < Real.pi := Real.pi_pos
  have hpi : Real.pi < 3.141593 := Real.pi_lt_d6
  have hsq : Real.pi ^ 2 < 3.141593 ^ 2 := by nlinarith
  have hfrac : 9 / (8 * (3.141593 : ℝ) ^ 2) < 9 / (8 * Real.pi ^ 2) := by
    apply div_lt_div_of_pos_left (by norm_num) (by positivity) (by nlinarith)
  have hJK : (40 / 49 - 6 * 0.000335 + 6 * 0.00859 : ℝ) ≤ 40 / 49 - 6 * J + 6 * K := by linarith
  have hpos : (0 : ℝ) < 40 / 49 - 6 * 0.000335 + 6 * 0.00859 := by norm_num
  calc (0.6931 : ℝ) ^ 3
        < 27 / 8 * (9 / (8 * (3.141593 : ℝ) ^ 2)) * (40 / 49 - 6 * 0.000335 + 6 * 0.00859) := by
        norm_num
    _ ≤ 27 / 8 * (9 / (8 * Real.pi ^ 2)) * (40 / 49 - 6 * J + 6 * K) := by
        apply mul_le_mul (by linarith) hJK hpos.le (by positivity)

end ErdosStraus
