# Modulating the SMI fODF by anisotropy

Measurement report for `SMI.fODF_ModulationWeight` and `SMI.modulate_fODF`.
Everything below is measured, not assumed; the simulation is
`example_fODF_modulation.m`.

---

## 1. Why the fODF needs modulating at all

`out.plm` is stored in the normalized convention `p_00 = 1`, so

```
A(u) = 1/(4*pi) + sum_{l>=2,m} f_lm Y_lm(u)
```

and the fODF integrates to **1 in every voxel**. An isotropic CSF voxel and a
tightly coherent white matter voxel carry the same total mass. The `1/(4*pi)`
isotropic floor is a fixed offset, and

```
1/(4*pi) = 0.0796  >  0.05 = MRtrix's default iFOD2 -cutoff
```

**An unmodulated SMI fODF passes MRtrix's termination test in every voxel of the
brain, CSF included.** Measured across all four conditions below, 100% of GM and
100% of CSF voxels survive the default cutoff. The tracker has no amplitude
information left to terminate on. That is the gap modulation closes.

Two distinct operations, which must not be conflated:

| | operation | effect |
|---|---|---|
| **density** | scale all coefficients including `l=0` | `fODF -> w*fODF`. Shape and peak orientation unchanged, mass becomes `w`. AFD-like; **this is what tractography needs.** |
| **shape** | scale only `l>0` | mass stays 1, fODF sharpens or flattens. This is what the Tikhonov term already does implicitly. |

`SMI.modulate_fODF` defaults to `density`.

---

## 2. Simulation

Seven voxel classes x 200 noise realizations, 6 shells
(b = 0,1,2,3,4.5,6 ms/um^2, `Lmax = [0 2 2 2 4 4]`), Rician noise, fitted through
the real `SMI.fit`, at SNR 30 and 15, with and without the regularized
deconvolution (`lambda_nonneg = 10`, `lambda_tikhonov = 0.3`).

| class | kernel `[f Da Depar Deperp fw]` | fODF |
|---|---|---|
| WM single | `0.60 2.0 2.0 0.50 0.05` | Watson, kappa 16 |
| WM crossing 60 deg | same | two Watsons at 60 deg |
| **WM in edema** | `0.30 2.0 2.2 1.00 0.50` | Watson, kappa 12 |
| GM | `0.15 1.2 1.0 0.80 0.10` | Watson, kappa 0.8 |
| CSF | `0.02 2.0 3.0 3.00 0.95` | isotropic |
| WM/CSF interface | 50/50 signal mix | independent orientations |
| WM/GM interface | 50/50 signal mix | independent orientations |

The edema class is the control that matters. Its `f` has collapsed and its `fw`
is 0.50, but its fibres are still coherent. **Any weight that suppresses it along
with GM and CSF is a tissue-type criterion in disguise and must be rejected** —
it would delete exactly the peritumoral region the tractography is for.

The forward model is the same construction `SMI.get_plm_from_S_and_kernel`
inverts, so it round-trips to `max|error| = 1.2e-15` with the true kernel and no
noise. Ground-truth `p2`/`p4` reproduce `out.pl` to 4 decimals.

---

## 3. The four candidate anisotropy maps

Medians per class. `p2_t`/`p4_t` are ground truth; `k_p2`/`k_p4` are
`out.kernel`; `pl2`/`pl4` are `out.pl`.

**SNR 30, regularized**

| class | p2_t | p4_t | k_p2 | k_p4 | pl2 | pl4 | fODF peak |
|---|---|---|---|---|---|---|---|
| WM single | 0.90 | 0.71 | 0.725 | 0.505 | 0.654 | 0.222 | 0.496 |
| WM crossing 60 | 0.60 | 0.42 | 0.498 | 0.329 | 0.470 | 0.135 | 0.287 |
| WM in edema | 0.87 | 0.63 | 0.605 | 0.402 | 0.537 | 0.276 | 0.489 |
| GM | 0.11 | 0.01 | 0.098 | 0.157 | 0.088 | 0.149 | 0.168 |
| **CSF** | **0.00** | **0.00** | **0.306** | **0.196** | **0.076** | 0.149 | 0.169 |
| WM/CSF interface | 0.90 | 0.71 | 0.618 | 0.419 | 0.564 | 0.270 | 0.494 |
| WM/GM interface | 0.90 | 0.71 | 0.576 | 0.461 | 0.551 | 0.275 | 0.493 |

Three things to read off this table.

**The kernel `p2` has a floor in CSF.** True value 0, estimated 0.306. The
polynomial regression is trained on a prior of tissue kernels; where the `l=2`
signal carries no information it returns roughly the prior mean rather than 0.
`pl2` has no such floor (0.076) because it is a direct projection of the measured
signal. This is the opposite of what one might guess — the kernel `p2` is the
*stabler* estimate but the *biased* one.

**`p4` is worse than `p2` on both counts.** The kernel `p4` carries the same
upward bias (0.196 in CSF against a true 0), and it is bounded above by the
training prior, which draws `p4 = rand * p2 * 0.9`. The deconvolved `pl4` is
destroyed by the Tikhonov term, which damps `l=4` far harder than `l=2`: 0.222 in
single-fibre WM against a true 0.71. **Weighting by `pl4` after a regularized fit
measures the regularizer, not the tissue.**

**The edema class keeps a high `p2`.** 0.605 (kernel) and 0.537 (deconvolved),
against 0.725/0.654 for healthy single-fibre WM — a modest drop, while its `f`
has halved. This is the empirical answer to the open question in the handoff: `p2`
does hold up in edematous voxels.

---

## 4. Where the 1e13 peaks come from

At SNR 15 **without** regularization, the CSF class produced
`peak p99 = 1.9e12` and `peak max = 2.2e14` — reproducing the 1e13 amplitudes
seen in real data. The mechanism is the deconvolution's noise-amplification gain

```
g_2 = 1 / || K_2(b) ||
```

measured at 3.6 in CSF against 0.6 in WM. `p_2m` is recovered by dividing the
`l=2` signal harmonic by `K_2`, and in a nearly isotropic voxel `K_2 -> 0`.

**The regularized deconvolution removes this completely**: same voxels, same
noise, `peak max = 0.449`. That is independent confirmation that the existing
regularization work is load-bearing, and modulation should be applied on top of
it, not instead of it.

---

## 5. Scoring

### 5.1 Fraction of voxels left above the MRtrix default cutoff of 0.05

Weights clipped at 1 (a physical `p_l` cannot exceed 1; a larger value means the
deconvolution failed and must not amplify the voxel).

**SNR 30, regularized**

| class | none | kernel p2 | kernel p2^2 | pl2 | pl2^2 | **p2product** |
|---|---|---|---|---|---|---|
| WM single | 100% | 100% | 100% | 100% | 100% | 100% |
| WM crossing 60 | 100% | 100% | 100% | 100% | 100% | 100% |
| WM in edema | 100% | 100% | 100% | 100% | 100% | 100% |
| GM | 100% | 0% | 0% | 0% | 0% | **0%** |
| CSF | 100% | 65% | 0% | 0% | 0% | **0%** |
| WM/CSF interface | 100% | 100% | 100% | 100% | 100% | 100% |
| WM/GM interface | 100% | 100% | 100% | 100% | 100% | 100% |

**SNR 15, regularized**

| class | none | kernel p2 | kernel p2^2 | pl2 | pl2^2 | **p2product** |
|---|---|---|---|---|---|---|
| WM single | 100% | 100% | 100% | 100% | 100% | 100% |
| WM crossing 60 | 100% | 100% | 94% | 100% | 86% | **92%** |
| WM in edema | 100% | 100% | 100% | 100% | 100% | 100% |
| GM | 100% | 0% | 0% | 0% | 0% | **0%** |
| CSF | 100% | 54% | 0% | 20% | 1% | **0%** |
| WM/CSF interface | 100% | 100% | 100% | 100% | 100% | 100% |
| WM/GM interface | 100% | 100% | 100% | 100% | 100% | 100% |

`kernel_p2` alone is not enough — it leaves 54-80% of CSF above cutoff across the
four conditions, because of its CSF floor. `pl2` alone leaves 20-41% at SNR 15,
because of its noise tail. The product removes both, reaching 0% in three of the
four conditions and 4.5% in the fourth (SNR 15 unregularized).

`p4` is not shown because it is disqualified: `pl4` at SNR 15 regularized leaves
**84% of CSF** above cutoff, which is worse than not weighting at all.

### 5.2 Adaptive cutoff test

Threshold `T` chosen to retain 95% of all white matter **including the edema
class**, then the fraction of GM+CSF still above `T`. A weight that only wins by
discarding edema cannot win here, because `T` is forced down to keep it.

**SNR 15, regularized**

| weight | FP GM+CSF | keep edema | separation |
|---|---|---|---|
| none | 10.5% | 100% | 1.39 |
| kernel p2 | 0.0% | 100% | 2.45 |
| kernel p2^2 | 0.0% | 100% | 4.25 |
| pl2 | 0.8% | 100% | 2.25 |
| pl4 | 60.0% | 100% | 1.52 |
| **p2product** | **0.0%** | **100%** | **4.15** |
| `f` *(tissue-type reference)* | 0.0% | **94%** | 2.58 |
| `1-fw` *(tissue-type reference)* | 31.5% | **85%** | 0.86 |

The last two rows are the quantitative case against tissue-fraction weights.
`f` and `1-fw` discard 6% and 15% of the edema class respectively, while every
`p2`-based weight keeps 100% of it. `1-fw` additionally fails at its own job,
leaving 63% of GM above threshold.

Separation across all four conditions for the default weight: 27.3 / 24.1
(SNR 30 unreg/reg) and 2.0 / 4.2 (SNR 15 unreg/reg), against an unweighted
baseline of 2.6 / 2.3 / 0.9 / 1.4.

---

## 6. Recommendation

**Default: `source = 'p2product'`, `mode = 'density'`, `exponent = 1`** — the
product of the two independent `p2` estimates, `out.kernel(p2) .* out.pl(:,:,:,1)`.
A voxel is kept only if both estimates agree it is coherent, which cancels the
kernel estimate's CSF floor against the deconvolved estimate's noise tail.

Two equivalent ways to apply it. As a flag on the fit, so the same script can be
run with and without:

```matlab
options.fODF_modulation.flag_modulate = 1;         % default 0, off
out = SMI.fit(dwi, options);
out.fODF_modulated                                  % the weighted SH coefficients
out.fODF_modulation.weight                          % the per voxel weight map
```

or post hoc, on an already fitted `out`:

```matlab
[sh, w, info] = SMI.modulate_fODF(out);            % defaults
[sh, w, info] = SMI.modulate_fODF(out, struct('source','kernel_p2','exponent',2));
```

The two paths produce bit-identical output. **`out.plm`, `out.pl` and
`out.kernel` are identical whether the flag is on or off** — the modulated fODF
is returned separately in `out.fODF_modulated`, so switching the flag cannot
perturb anything upstream, and a modulated and an unmodulated run can be
compared directly. Settings and the degenerate voxel count are written to the
log file, as they are for the regularization.

- `kernel_p2` with `exponent = 2` is the near-equal alternative, and the only
  option available when `flag_fit_fODF = 0`.
- **Do not use `p4`**, in either flavour.
- **Run the regularized deconvolution.** Modulation does not repair a blown-up
  voxel: the weight is clipped at 1, so `w * 1e13` is still `1e13`. The
  regularizer is what prevents the blow-up.

### Caveats

- These are simulations. `p2` in real edema should be checked against a known
  ROI before this is trusted clinically; the simulated edema class assumes fibres
  stay coherent, which is the assumption the whole approach rests on.
- Squaring the weight costs crossing fibres at low SNR (94% retained at SNR 15,
  vs 100% at exponent 1). Crossings have genuinely lower `p2`, so any exponent
  above 1 penalizes them. Keep `exponent = 1` unless CSF suppression is
  insufficient.
- The fixed 0.05 comparisons assume the MRtrix cutoff is left at its default.
  Modulation lowers the overall amplitude scale, so the cutoff may want retuning;
  the adaptive test in 5.2 is the scale-free comparison.

### Order of operations

1. `SMI.fit` with the regularized deconvolution
2. peak truncation, if used — its thresholds are **absolute amplitudes**, so it
   must run **before** modulation or be re-derived after
3. `SMI.modulate_fODF`
4. write NIfTI

### Convention warning

In `density` mode `p_00` is no longer 1 — that is the point, but it means the
result is **not** in the convention the rest of the toolbox assumes. `sh`
includes the `l=0` term, which `out.plm` does not store; without it the density
weighting would be lost. Document which convention a saved NIfTI is in.

`sh` is in SMI's own SH basis (`SMI.get_even_SH` with `out.CS_phase`). **Verify
the coefficient ordering and the Condon-Shortley convention against MRtrix's
basis before feeding it to `tckgen`** — this has not been checked here.

---

## 7. Files

| file | contents |
|---|---|
| `SMI.m` | `grab_pl`, `grab_kernel_pl`, `fODF_ModulationWeight`, `modulate_fODF` (additions only; existing code paths untouched, so fits without modulation are bit-identical) |
| `example_fODF_modulation.m` | the simulation and scoring above |
| `REPORT_fODF_modulation.md` | this file |

`out.kernel` layout is `[f Da Depar Deperp fw (T2a T2e) p2 (p4) (p6)]`, so `p2`
sits at index 6 without T2 fitting and 8 with it. The map count alone is
ambiguous (7 maps is either `[f Da Depar Deperp fw p2 p4]` or
`[f Da Depar Deperp fw T2a T2e]`), so variable TE is detected from `out.shells`
row 4. Override with `options.kernel_p2_index` if that is ever wrong.
