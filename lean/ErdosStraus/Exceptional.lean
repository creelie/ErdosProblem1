import Mathlib

/-!
# The algebra behind the explicit exponent for the exceptional set

These are the identities and inequalities used in the paper
"An explicit exponent for the exceptional set of the Erdős–Straus equation".
A Type II triple for a prime `p = 3 (mod 4)` is `(u, a, α)` with `u a α = P = (p + 1)/4`.

* `class_M_form`, `class_I_form`: the class `u/a (mod p)` equals `4 u² α` and `(4 a² α)⁻¹`.
* `corner_AB_shift`, `corner_AB_divisor`: the two identities of the comparison of the corners
  `A` and `B`.
* `corner_BC_main`, `corner_BC_t`, `corner_BC_M`, `corner_BC_quadratic`: the parametrisation of a
  coincidence between the corners `B` and `C`.
* `corner_BC_tau`, `corner_BC_G`, `corner_BC_P`: the inequalities that make the coincidences rare.
* `totient_mul_of_dvd`, `local_factor_B`: the mean value of the corner `B`.
* `tilt_value`, `exponent_cube`, `exponent_cube_rpow`: the value of the exponent.
-/

namespace ErdosStraus

open Nat

/-- The M-form: `4 (u² α) a - u = u p` with `p = 4 u a α - 1`, so `u/a ≡ 4 u² α (mod p)`. -/
theorem class_M_form (u a α : ℤ) : 4 * (u ^ 2 * α) * a - u = u * (4 * u * a * α - 1) := by
  ring

/-- The I-form: `4 (a² α) u - a = a p`, so `u/a ≡ (4 a² α)⁻¹ (mod p)`. -/
theorem class_I_form (u a α : ℤ) : 4 * (a ^ 2 * α) * u - a = a * (4 * u * a * α - 1) := by
  ring

/-- Corners `A` and `B`: if `4 a m - u = w p` with `w = u - 4 a j`, then `m = j + w u α`. -/
theorem corner_AB_shift (u a α j : ℤ) :
    (u - 4 * a * j) * (4 * u * a * α - 1) + u = 4 * a * (j + (u - 4 * a * j) * u * α) := by
  ring

/-- Corners `A` and `B`: with `m = j + w u α`, `P = u a α = a' d` and `m = u' d`,
one has `d (a u' - a' w) = a j`. -/
theorem corner_AB_divisor (u a α j w m d u' a' : ℤ) (hm : m = j + w * u * α)
    (hP : u * a * α = a' * d) (hm' : m = u' * d) : d * (a * u' - a' * w) = a * j := by
  linear_combination a * hm - a * hm' + w * hP

/-- Corners `B` and `C`: from `P² = d e (M + r p)` with `p = 4P - 1`,
`P (P - 4 d e r) = d e (M - r)`. -/
theorem corner_BC_main (P d e M r : ℤ) (h : P ^ 2 = d * e * (M + r * (4 * P - 1))) :
    P * (P - 4 * d * e * r) = d * e * (M - r) := by
  linear_combination h

/-- Corners `B` and `C`: with `d = G d₁`, `e = G e₁`, `W = 4 G r + τ` and `P = G d₁ e₁ W`,
`P - 4 d e r = G d₁ e₁ τ`. -/
theorem corner_BC_t (G d₁ e₁ r τ : ℤ) :
    G * d₁ * e₁ * (4 * G * r + τ) - 4 * (G * d₁) * (G * e₁) * r = G * d₁ * e₁ * τ := by
  ring

/-- Corners `B` and `C`: the value of `M = u a'`. -/
theorem corner_BC_M (G d₁ e₁ r τ M : ℤ) (hG : G ≠ 0) (hd : d₁ ≠ 0) (he : e₁ ≠ 0)
    (h : (G * d₁ * e₁ * (4 * G * r + τ)) * (G * d₁ * e₁ * τ) = (G * d₁) * (G * e₁) * (M - r)) :
    M = d₁ * e₁ * τ * (4 * G * r + τ) + r := by
  have hne : G * G * d₁ * e₁ ≠ 0 := by
    apply mul_ne_zero (mul_ne_zero (mul_ne_zero hG hG) hd) he
  have h2 : (G * G * d₁ * e₁) * (d₁ * e₁ * τ * (4 * G * r + τ)) = (G * G * d₁ * e₁) * (M - r) := by
    linear_combination h
  have h3 := mul_left_cancel₀ hne h2
  linear_combination -h3

/-- Corners `B` and `C`: the quadratic equation for `τ`, with `k = α α'`, `k M = G² d₁ e₁`
and `G = 4 k τ r + m`. -/
theorem corner_BC_quadratic (G d₁ e₁ r τ M k m : ℤ) (hk : k * M = G ^ 2 * d₁ * e₁)
    (hM : M = d₁ * e₁ * τ * (4 * G * r + τ) + r) (hG : G = 4 * k * τ * r + m) :
    d₁ * e₁ * (m * (4 * k * τ * r + m)) = d₁ * e₁ * (k * τ ^ 2) + k * r := by
  subst hG
  rw [hM] at hk
  linear_combination -hk

/-- Corners `B` and `C`: `τ ≥ 4 m r`. -/
theorem corner_BC_tau (k m r τ : ℤ) (hk : 1 ≤ k) (hm : 1 ≤ m) (hr : 1 ≤ r) (hτ : 1 ≤ τ)
    (h : 4 * k * τ * r * m + m ^ 2 ≤ k * τ ^ 2 + k * r) : 4 * m * r ≤ τ := by
  by_contra hc
  have hc : τ < 4 * m * r := not_le.mp hc
  have h1 : 0 ≤ (τ - 1) * (4 * m * r - 1 - τ) := by
    apply mul_nonneg <;> linarith
  have h2 : 4 * m * r - 1 ≤ τ * (4 * m * r - τ) := by nlinarith
  have h3 : 3 * r ≤ 4 * m * r - 1 := by nlinarith
  have h4 : r < τ * (4 * m * r - τ) + 1 := by linarith
  have h5 : k * (τ * (4 * m * r - τ)) < k * r := by
    have : 0 < m ^ 2 := by positivity
    nlinarith
  have h6 : τ * (4 * m * r - τ) < r := by
    by_contra h7
    have h7 : r ≤ τ * (4 * m * r - τ) := not_lt.mp h7
    have : k * r ≤ k * (τ * (4 * m * r - τ)) := by
      apply mul_le_mul_of_nonneg_left h7; linarith
    linarith
  linarith

/-- Corners `B` and `C`: `G ≥ 16 k m r²`. -/
theorem corner_BC_G (k m r τ G : ℤ) (hk : 1 ≤ k) (hm : 1 ≤ m) (hr : 1 ≤ r)
    (hτ : 4 * m * r ≤ τ) (hG : G = 4 * k * τ * r + m) : 16 * k * m * r ^ 2 ≤ G := by
  subst hG
  have h1 : 0 ≤ 4 * k * r := by positivity
  nlinarith

/-- Corners `B` and `C`: `P ≥ 1024 k² m² r⁵ d₁ e₁`. -/
theorem corner_BC_P (k m r τ G d₁ e₁ : ℤ) (hk : 1 ≤ k) (hm : 1 ≤ m) (hr : 1 ≤ r) (hτ : 1 ≤ τ)
    (hd : 1 ≤ d₁) (he : 1 ≤ e₁) (hG : 16 * k * m * r ^ 2 ≤ G) :
    1024 * k ^ 2 * m ^ 2 * r ^ 5 * (d₁ * e₁) ≤ G * d₁ * e₁ * (4 * G * r + τ) := by
  have hG0 : 0 ≤ 16 * k * m * r ^ 2 := by positivity
  have hde : 0 ≤ d₁ * e₁ := by positivity
  have h1 : (16 * k * m * r ^ 2) ^ 2 ≤ G ^ 2 := by nlinarith
  have h2 : 4 * r * (16 * k * m * r ^ 2) ^ 2 ≤ 4 * r * G ^ 2 := by
    apply mul_le_mul_of_nonneg_left h1; positivity
  have h3 : 4 * r * G ^ 2 ≤ G * (4 * G * r + τ) := by nlinarith
  have h4 : 1024 * k ^ 2 * m ^ 2 * r ^ 5 = 4 * r * (16 * k * m * r ^ 2) ^ 2 := by ring
  have h5 : 1024 * k ^ 2 * m ^ 2 * r ^ 5 ≤ G * (4 * G * r + τ) := by linarith
  have h6 := mul_le_mul_of_nonneg_right h5 hde
  linarith [show G * d₁ * e₁ * (4 * G * r + τ) = G * (4 * G * r + τ) * (d₁ * e₁) by ring]

/-- `φ(q n) = q φ(n)` when `q` divides `n`. -/
theorem totient_mul_of_dvd (q n : ℕ) (hq : q ∣ n) : φ (q * n) = q * φ n := by
  rcases Nat.eq_zero_or_pos q with h0 | hpos
  · subst h0; simp
  have h := Nat.totient_gcd_mul_totient_mul q n
  rw [Nat.gcd_eq_left hq] at h
  have hφ : 0 < φ q := Nat.totient_pos.2 hpos
  have : φ q * φ (q * n) = φ q * (q * φ n) := by rw [h]; ring
  exact Nat.eq_of_mul_eq_mul_left hφ this

/-- The local factor of the corner `B` at an odd prime: the sum over `p^i ∥ u`, `p^j ∥ α`
equals `(p² + 1)/(p - 1)²`, the same factor as for the corner `A`. -/
theorem local_factor_B (p : ℝ) (hp : 1 < p) :
    (1 / (p - 1)) * (p / (p - 1)) + (p / (p - 1)) * (1 / (p - 1)) + 1
      = (p ^ 2 + 1) / (p - 1) ^ 2 := by
  have h : p - 1 ≠ 0 := by linarith
  field_simp
  ring

/-- The local factor of the corner `B` at the prime `2`: `(1/2) · 1 · 2 + 2 = 3`, which is
`3/4` of the factor `4` of `ζ(1)²`. -/
theorem local_two_B : (1 / 2 : ℝ) * 1 * 2 + 2 = 3 / 4 * 4 := by norm_num

/-- The optimum of the tilt: if `σ³ T = 4a` with `σ > 0`, then `σ T + 2a/σ² = (3/2) σ T`. -/
theorem tilt_value (σ T a : ℝ) (hσ : 0 < σ) (h : σ ^ 3 * T = 4 * a) :
    σ * T + 2 * a / σ ^ 2 = 3 / 2 * (σ * T) := by
  have hσ2 : σ ^ 2 ≠ 0 := by positivity
  field_simp
  nlinarith [h]

/-- The cube of the exponent: `(27/8) · (3κ/16) = 729/(256 π²)` with `κ = 9/(2π²)`. -/
theorem exponent_cube :
    (27 / 8 : ℝ) * (3 / 16 * (9 / (2 * Real.pi ^ 2))) = 729 / (256 * Real.pi ^ 2) := by
  have : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp
  ring

/-- `c₀ = (3/2) (27/(32π²))^{1/3}` satisfies `c₀³ = 729/(256 π²)`. -/
theorem exponent_cube_rpow :
    ((3 / 2 : ℝ) * (27 / (32 * Real.pi ^ 2)) ^ ((1 : ℝ) / 3)) ^ 3
      = 729 / (256 * Real.pi ^ 2) := by
  have hpos : (0 : ℝ) ≤ 27 / (32 * Real.pi ^ 2) := by positivity
  have h1 : ((27 / (32 * Real.pi ^ 2) : ℝ) ^ ((1 : ℝ) / 3)) ^ 3 = 27 / (32 * Real.pi ^ 2) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hpos]
    norm_num
  rw [mul_pow, h1]
  have : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp
  ring

/-- The coincidence inside the centre: `81² - 328 · 1640 = -531359`, so `81/328 ≡ 1640/81`
modulo the prime `531359`, with `81 · 328 · 5 = 1640 · 81 · 1 = 132840 = (531359 + 1)/4`. -/
theorem centre_coincidence :
    (81 : ℤ) ^ 2 - 328 * 1640 = -531359 ∧ (81 : ℤ) * 328 * 5 = 132840 ∧
      (1640 : ℤ) * 81 * 1 = 132840 ∧ 4 * (132840 : ℤ) - 1 = 531359 := by
  norm_num

/-- `531359` is prime. -/
theorem centre_prime : Nat.Prime 531359 := by norm_num

end ErdosStraus
