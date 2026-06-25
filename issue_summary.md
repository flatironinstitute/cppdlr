# Symmetrized DLR conditioning → tprf Eliashberg failure

## Summary
The symmetrized DLR construction selected real-frequency nodes only in exact
mirror pairs. With generically small spectral weight near zero, this forced a
near-degenerate pole pair `±δ` straddling `ω=0`, ill-conditioning the DLR fit
matrices (`cond ≈ 9e9` vs `≈ 2e8` non-sym at `λ=10, ε=1e-8`). Harmless on
spectral inputs, but tprf's Eliashberg eigensolver applies a DLR-built operator
to **non-spectral** Krylov vectors, where the ill-conditioning makes the
operator lose linearity in floating point — so ARPACK returns spurious,
non-reproducible eigenvalues.

**Fix:** always include each grid's self-symmetric fixed point as a single
self-paired node (`ω=0` real-freq, `τ=β/2` imtime, `n=0` bosonic Matsubara).
This removes the near-degeneracy and restores conditioning ≈ non-symmetric.

## Mechanism
- `vals2coefs` is a solve `c = K⁻¹v`; components of `v` along the bad singular
  vector are amplified by `κ ≈ 10¹⁰`.
- A spectral input has ≈0 weight there (why DLR is only valid for
  Lehmann-representable objects); a generic/Krylov vector does not, so its
  coefficients blow up as `~κ‖v‖` (large `±δ`-pole weights that cancel *at the
  DLR nodes*).
- The operator then evaluates those coefficients in a *different* representation
  (imtime, the bubble, another frequency set) where the cancellation no longer
  holds → catastrophic cancellation, losing `~log₁₀κ ≈ 10` of ~16 digits → O(1)
  error. So although `M` is linear on paper, in float `M(a+b) − Ma − Mb` is O(1).
- ARPACK assumes a faithful linear operator; O(1) inconsistency seeded by the
  random start gives spurious, run-to-run non-reproducible eigenvalues (the
  leading value jumping instead of sitting at 0.5639).

## Consequences (cppdlr)
- Symmetrized rank `r` is now odd.
- Symmetrized **fermionic** Matsubara has no self-symmetric frequency → decouples
  to even `niom = r+1` (mirror pairs, least-squares). Non-symmetrized unchanged.

## Failing reproducer (this branch)
`test/c++/imfreq_ops.cpp::TEST(imfreq_ops, symmetrized_operator_linearity)` is a
**failing** reproducer on this branch (pre-fix). It replays tprf's Eliashberg
matvec (`eliashberg_product_fft`) reduced to one k-point/scalar, applies it to a
random non-spectral vector pair at `β=2, λ=10, ε=1e-8`, and checks
`M(a+b) = M(a) + M(b)`: the symmetric operator's residual is **~10⁴×** the
non-symmetric one (~7e-5 vs ~6e-9) — a manifestly non-linear operator — so the
test fails. With the self-symmetric-fixed-point fix it drops to ~6e-9 (≈
non-symmetric) and the test passes.

The blow-up factor comes from the `convolve`, but the effect does not depend on
it: even a single `vals2coefs`(imfreq) → `coefs2vals`(imtime) conversion already
loses linearity (~30×); the rest of the matvec only amplifies it.
