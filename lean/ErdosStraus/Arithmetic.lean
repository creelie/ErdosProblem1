import Mathlib

/-!
# The arithmetic behind the exact constant of the sum `W(J)`

For coprime positive integers `u` and `a`, the totient of `4ua` splits as
`φ(4ua) = c · φ(u) · φ(a)` with `c = 4` when `ua` is even and `c = 2` when `ua` is odd.
This is the factorisation that makes `ua/φ(4ua)` multiplicative in the pair `(u, a)`.
-/

namespace ErdosStraus

open Nat

/-- `φ(4n) = 2 φ(n)` for odd `n`. -/
theorem totient_four_mul_of_odd (n : ℕ) (hn : Odd n) : φ (4 * n) = 2 * φ n := by
  have hc : Nat.Coprime 4 n := by
    have h2 : Nat.Coprime 2 n := (Nat.coprime_two_left).2 hn
    have : (4 : ℕ) = 2 ^ 2 := by norm_num
    rw [this]
    exact Nat.Coprime.pow_left 2 h2
  rw [Nat.totient_mul hc]
  rfl

/-- `φ(4n) = 4 φ(n)` for even `n`. -/
theorem totient_four_mul_of_even (n : ℕ) (hn : 2 ∣ n) : φ (4 * n) = 4 * φ n := by
  have h1 : φ (2 * n) = 2 * φ n := Nat.totient_mul_of_prime_of_dvd Nat.prime_two hn
  have h2 : 2 ∣ 2 * n := dvd_mul_right 2 n
  have h3 : φ (2 * (2 * n)) = 2 * φ (2 * n) := Nat.totient_mul_of_prime_of_dvd Nat.prime_two h2
  have : 4 * n = 2 * (2 * n) := by ring
  rw [this, h3, h1]
  ring

/-- For coprime `u` and `a`, `φ(4ua) = c · φ(u) · φ(a)`, with `c = 4` if `ua` is even
and `c = 2` if `ua` is odd. -/
theorem totient_four_mul_mul (u a : ℕ) (h : Nat.Coprime u a) :
    φ (4 * (u * a)) = (if 2 ∣ u * a then 4 else 2) * (φ u * φ a) := by
  rw [← Nat.totient_mul h]
  split_ifs with h2
  · exact totient_four_mul_of_even _ h2
  · exact totient_four_mul_of_odd _ (Nat.odd_iff.2 (Nat.two_dvd_ne_zero.1 h2))

/-! ## Local factors at an odd prime `p` -/

/-- Inverting `1 ⊗ 1` at an odd prime: the local values `k(p,1) = k(1,p) = 1/(p-1)`,
`k(p,p) = -(p+1)/(p-1)` sum to `h₀(p^i,1) = p/(p-1)` and to `h₀(p^i,p^j) = 0`. -/
theorem local_inversion_axis (p : ℝ) (hp : 1 < p) :
    1 + 1 / (p - 1) = p / (p - 1) := by
  have : p - 1 ≠ 0 := by linarith
  field_simp
  ring

theorem local_inversion_diag (p : ℝ) (hp : 1 < p) :
    1 + 1 / (p - 1) + 1 / (p - 1) - (p + 1) / (p - 1) = 0 := by
  have : p - 1 ≠ 0 := by linarith
  field_simp
  ring

/-- The local factor of `Z(0,0) = Σ k(d,e)/(de)` at an odd prime is `1 + 1/p²`. -/
theorem local_factor_odd (p : ℝ) (hp : 1 < p) :
    1 + 1 / (p - 1) * (1 / p) + 1 / (p - 1) * (1 / p) - (p + 1) / (p - 1) * (1 / p ^ 2)
      = 1 + 1 / p ^ 2 := by
  have h1 : p - 1 ≠ 0 := by linarith
  have h2 : p ≠ 0 := by linarith
  field_simp
  ring

/-- The same local factor seen from the other side: `(1 - 1/p)² · Σ_{(i,j)} h₀/(p^{i+j})`. -/
theorem local_factor_odd' (p : ℝ) (hp : 1 < p) :
    (1 - 1 / p) ^ 2 * (1 + 2 * (p / (p - 1)) * (1 / (p - 1))) = 1 + 1 / p ^ 2 := by
  have h1 : p - 1 ≠ 0 := by linarith
  have h2 : p ≠ 0 := by linarith
  field_simp
  ring

/-- Logarithmic derivative data at an odd prime: with `F = 1 + 1/p²`,
`F_x = -1/((p-1)p²)` and `F_{xy} = -(p+1)/((p-1)p²)`. -/
theorem local_Fx (p : ℝ) (hp : 1 < p) :
    1 / (p * (p - 1)) - (p + 1) / ((p - 1) * p ^ 2) = -1 / ((p - 1) * p ^ 2) := by
  have h1 : p - 1 ≠ 0 := by linarith
  have h2 : p ≠ 0 := by linarith
  field_simp
  ring

theorem local_T (p : ℝ) (hp : 1 < p) :
    (1 / ((p - 1) * p ^ 2)) / (1 + 1 / p ^ 2) = 1 / ((p - 1) * (p ^ 2 + 1)) := by
  have h1 : p - 1 ≠ 0 := by linarith
  have h2 : p ≠ 0 := by linarith
  have h3 : p ^ 2 + 1 ≠ 0 := by positivity
  field_simp

theorem local_M (p : ℝ) (hp : 1 < p) :
    (-(p + 1) / ((p - 1) * p ^ 2)) / (1 + 1 / p ^ 2)
        - ((-1 / ((p - 1) * p ^ 2)) / (1 + 1 / p ^ 2)) ^ 2
      = -((p + 1) / ((p - 1) * (p ^ 2 + 1)) + 1 / ((p - 1) ^ 2 * (p ^ 2 + 1) ^ 2)) := by
  have h1 : p - 1 ≠ 0 := by linarith
  have h2 : p ≠ 0 := by linarith
  have h3 : p ^ 2 + 1 ≠ 0 := by positivity
  field_simp
  ring

/-! ## The prime `2` -/

/-- At `p = 2` the local values are `k(2,2) = -1` and zero otherwise; they sum to
`h₀(2^i,1) = 1` and `h₀(2^i,2^j) = 0`, and the local factor is `3/4`. -/
theorem local_two : ((1 : ℝ) + 0 = 1) ∧ ((1 : ℝ) + 0 + 0 - 1 = 0) ∧ ((1 : ℝ) - 1 / 4 = 3 / 4) := by
  norm_num

/-- At `p = 2`: `T₂ = (log 2)/3` and `M₂ = -(4/9)(log 2)²` come from
`F = 3/4`, `F_x = F_{xy} = -1/4`. -/
theorem local_two_derivatives :
    -((-1 / 4 : ℝ) / (3 / 4)) = 1 / 3 ∧ (-1 / 4 : ℝ) / (3 / 4) - ((-1 / 4) / (3 / 4)) ^ 2 = -4 / 9 := by
  norm_num

/-! ## The constant -/

/-- `Π_{p>2} (1 + p^{-2}) = (ζ(2)/ζ(4)) / (1 + 1/4)`, and
`(1/2) · (3/4) · (ζ(2)/ζ(4)) / (5/4) = 9/(2π²)` with `ζ(2) = π²/6`, `ζ(4) = π⁴/90`. -/
theorem constant_value :
    (1 / 2 : ℝ) * (3 / 4) * ((Real.pi ^ 2 / 6) / (Real.pi ^ 4 / 90)) / (1 + 1 / 2 ^ 2)
      = 9 / (2 * Real.pi ^ 2) := by
  have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp
  ring

/-- `ζ(2)` and `ζ(4)` as used above. -/
theorem zeta_values :
    riemannZeta 2 = (Real.pi : ℂ) ^ 2 / 6 ∧ riemannZeta 4 = (Real.pi : ℂ) ^ 4 / 90 :=
  ⟨riemannZeta_two, riemannZeta_four⟩

/-- The leading constant is three times Dahan's constant `3/(2π²)`. -/
theorem three_times : 9 / (2 * Real.pi ^ 2) = 3 * (3 / (2 * Real.pi ^ 2)) := by ring

/-- The Type II identity: if `a n + u = v (4 a α u - 1)`, then
`4/n = 1/(α u v) + 1/(a α n u) + 1/(a α n v)`. -/
theorem typeII_identity (n u a v α : ℚ) (hn : n ≠ 0) (hu : u ≠ 0) (ha : a ≠ 0) (hv : v ≠ 0)
    (hα : α ≠ 0) (h : a * n + u = v * (4 * a * α * u - 1)) :
    4 / n = 1 / (α * u * v) + 1 / (a * α * n * u) + 1 / (a * α * n * v) := by
  field_simp
  linear_combination (-1 : ℚ) * h

/-- The first example: `n = 409`, `(u, a) = (1, 2)`, `s = 7`, `v = 117`, `α = 1`. -/
theorem example_409 : (4 : ℚ) / 409 = 1 / 117 + 1 / 818 + 1 / 95706 := by norm_num

end ErdosStraus
