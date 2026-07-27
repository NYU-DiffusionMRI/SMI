# fODF deconvolution regularization: parameter sweep report

Measurement of the regularization parameters of the SMI fODF deconvolution against a
synthetic ground truth, and the resulting change of defaults.

Everything below is produced by `example_fODF_regularization_sweep.m` (scoring in
`fODF_regularization_score.m`). It needs no data and reruns in a few minutes.

---

## 1. Why

The fODF is obtained by deconvolving the SM kernel from the DWI. That problem is
ill-conditioned: the kernel rotational invariants `Kl` decay quickly with `l`, so the high
order `plm` are noise-dominated and the fODF develops large negative lobes. Two optional
regularizers exist for it:

- **Non-negativity** (constrained spherical deconvolution, Tournier et al., NeuroImage
  2007), weight `lambda_nonneg`, threshold `tau`, initialized from a `Lmax_init` order fit.
- **Tikhonov** damping `lambda_tikhonov^2*||Gamma*plm||^2`, with `Gamma` the identity or the
  Laplace-Beltrami matrix `diag(l(l+1))/max(l(l+1))`.

Their weights were documented as "around 1 is a reasonable starting point". That was a
plausible round number, never measured. This sweep measures them.

---

## 2. Simulation setup

### 2.1 Ground truth

Each voxel is a **two-fiber crossing**. The fODF is a pair of smoothed delta functions,

```
plm_SH  = (weights * Ylm(fibers)) .* exp(-smoothing*l*(l+1))
plm     = plm_SH ./ sqrt((2l+1)/(4*pi))          % SMI normalized convention, p_00 = 1
```

with `weights = [0.6 0.4]` and `smoothing = 0.05`. The smoothing matters: it makes the
ground truth fODF genuinely non-negative (its own negative mass is 0.0034), so that a
non-negativity constraint is not being scored against a target that violates it.

The `plm` are therefore known **exactly**, not estimated. There is no reference-method bias
in the scoring.

### 2.2 Signal

Generated with the SMI forward model `SM_wFW_b_beta_TE_RealSphHarm`, so the simulation is
self-consistent with the code under test:

| | |
|---|---|
| Kernel | `f = 0.6`, `Da = 2`, `Depar = 2`, `Deperp = 0.5` um^2/ms, `fw = 0` |
| Shells | b = 1, 2, 3 ms/um^2, 64 directions each (golden-angle hemisphere) |
| b0 | 4 volumes |
| Total | 196 measurements, fixed TE |
| `Lmax` | 6 (27 `plm` coefficients for l = 2, 4, 6) |
| CS phase | 1 |
| Noise | Gaussian, `randn/SNR` on the normalized signal |

A noise-free sanity check runs first: the unregularized deconvolution recovers the ground
truth `plm` to `max|error| = 5e-16`, i.e. the forward and inverse models agree to machine
precision. Any error measured afterwards is caused by noise and by the regularization, not
by a modelling mismatch.

### 2.3 Conditions

Tuning on one geometry at one noise level produces parameters tuned to that case. The sweep
averages over a grid of conditions:

| | values |
|---|---|
| Crossing angle | 40, 60, 90 degrees |
| SNR | 20, 30, 50 |
| Noise realizations | 30 per condition |

= 9 conditions x 30 = **270 voxels scored per parameter setting**, with a fixed RNG seed so
the comparison between settings is paired (every setting sees the same noise).

### 2.4 Scores

| score | definition | role |
|---|---|---|
| **rel err fODF** | relative L2 error of the fODF amplitude over 500 hemisphere directions | **primary**, dimensionless, this is what "reconstructs the fODF like the ground truth" means |
| RMSE(plm) | RMS error of the coefficients | secondary |
| negative mass | fraction of the absolute fODF mass that is negative | physical plausibility |
| peak error | mean angle from each true fiber to the closest peak of the estimate | what tractography actually consumes |
| **peak amplitude ratio** | `peak(estimate)/peak(truth)` per voxel, averaged. `peak` is the maximum of the fODF over directions | what the regularizer *costs* in fODF height. 1 means the height is preserved, below 1 means the fODF is flattened |
| peaks/voxel | mean number of detected local maxima | fiber detection, 2 is correct above the resolution limit |

Peaks are local maxima on the 500 directions (a direction beating its 8 nearest neighbours,
neighbours by `|dot|` since even-order SH are antipodally symmetric), kept above 25% of the
voxel maximum. That threshold is necessary: a band-limited fODF has a local maximum on every
truncation ripple, and an Lmax=6 two-fiber fODF has 8 of them, of which only 2 are real.

### 2.5 Performance note

All 270 voxels of one parameter setting are deconvolved in a **single** call to
`get_plm_from_S_and_kernel`, by stacking the noise realizations along the first image
dimension. Measured 3.1x faster than one call per realization, and verified bit-identical
(difference exactly 0), because the per-call cost is dominated by building `Kell`, which
does not depend on the number of voxels.

---

## 3. Parameters tested

Staged rather than a full product, since the parameters are close to separable and the joint
grid is where the interaction actually lives.

| stage | parameter | values tested |
|---|---|---|
| 1 | `lambda_nonneg` x `lambda_tikhonov` (joint grid) | `[0 0.3 1 3 10 30 100]` x `[0 0.1 0.3 1 3 10]` |
| 2 | `tau` | `[0.02 0.05 0.1 0.2 0.4]` |
| 3 | `Lmax_init` | `[2 4 6]` |
| 4 | `TikhonovMatrix` (re-optimizing its weight) | `identity`, `laplacebeltrami` |

`lambda_nonneg = 0` means the non-negativity constraint is disabled, `lambda_tikhonov = 0`
means Tikhonov is disabled, so the grid contains the unregularized fit at its corner.

`Ndirs` (300) and `Niter` (50) were **not** swept: they control cost and the iteration limit,
not the quality of the solution. Convergence is reached in 4-6 iterations in practice, well
inside the limit.

---

## 4. Results

### 4.1 Stage 1, joint grid (rel err fODF, lower is better)

Rows `lambda_nonneg`, columns `lambda_tikhonov`:

| | 0 | 0.1 | 0.3 | 1 | 3 | 10 |
|---|---|---|---|---|---|---|
| **0** (off) | 0.6148 | 0.6044 | 0.5337 | 0.2860 | 0.3789 | 0.6881 |
| **0.3** | 0.4047 | 0.4008 | 0.3730 | 0.2574 | 0.3789 | 0.6886 |
| **1** | 0.2041 | 0.2034 | 0.1991 | 0.1964 | 0.3788 | 0.6886 |
| **3** | 0.1382 | 0.1383 | 0.1388 | 0.1659 | 0.3827 | 0.6886 |
| **10** | 0.1054 | 0.1055 | **0.1052** | 0.1203 | 0.3924 | 0.6886 |
| **30** | 0.1349 | 0.1349 | 0.1332 | 0.1138 | 0.3446 | 0.6886 |
| **100** | 0.4123 | 0.4138 | 0.4167 | 0.4202 | 0.3561 | 0.6886 |

The optimum is `lambda_nonneg = 10`, `lambda_tikhonov = 0.3`, at 0.1052.

**This optimum is interior, not a grid edge.** The first version of this sweep stopped at
`lambda_nonneg = 10` with the error still falling, which would have been an unusable result;
the grid was extended to 100 and the turnover is now explicit: 0.105 at 10, 0.133 at 30,
0.417 at 100. Past ~10 the penalty rows begin to dominate the data rows and the solution is
pulled away from the signal.

Note also that non-negativity is much the stronger of the two regularizers. Along the
`lambda_tikhonov` axis alone (top row) the best achievable is 0.286; along the
`lambda_nonneg` axis alone (first column) it is 0.105.

### 4.1b Peak amplitude: does regularization flatten the fODF?

Same grid, peak amplitude ratio `peak(estimate)/peak(truth)`:

| | 0 | 0.1 | 0.3 | 1 | 3 | 10 |
|---|---|---|---|---|---|---|
| **0** (off) | 1.116 | 1.110 | 1.065 | 0.865 | 0.577 | 0.272 |
| **0.3** | 1.063 | 1.060 | 1.032 | 0.867 | 0.577 | 0.270 |
| **1** | 1.019 | 1.017 | 1.001 | 0.873 | 0.576 | 0.270 |
| **3** | 0.988 | 0.987 | 0.976 | 0.873 | 0.571 | 0.270 |
| **10** | 1.005 | 1.004 | **0.998** | 0.924 | 0.550 | 0.270 |
| **30** | 1.058 | 1.057 | 1.055 | 1.008 | 0.613 | 0.270 |
| **100** | 1.302 | 1.304 | 1.305 | 1.292 | 1.080 | 0.270 |

Four things, none of them the naive expectation that "more regularization = smaller peaks":

1. **The unregularized fODF is inflated, not neutral.** Its ratio is 1.116, so the reference point is 12% too tall. The peak is a maximum over directions, and taking a max over a noisy field biases upward. Any comparison that treats the unregularized peak height as ground truth is starting from a biased baseline.

2. **Tikhonov is the flattening mechanism**, and it is strong and monotone: 1.065 at `lambda_tikhonov` = 0.3, 0.865 at 1, 0.577 at 3, 0.272 at 10. It damps the `plm` directly, so it shrinks the anisotropic part of the fODF by construction. This is what drives the large errors in the right hand columns of section 4.1.

3. **Non-negativity does not flatten.** Off to `lambda_nonneg` = 10 moves the ratio 1.116 -> 1.005: it removes the noise-driven inflation rather than pushing the fODF below the truth. That is consistent with what the constraint does mechanically, since it only penalizes amplitude on directions where the fODF is *negative* and leaves the lobes alone.

4. **Past its optimum, non-negativity inflates again.** At `lambda_nonneg` = 100 the ratio is 1.30. A constraint that strong forces the fODF to zero over most of the sphere, and the remaining mass concentrates into narrow spikes that are taller than the true lobes. This is the same failure the accuracy table shows as 0.412 at that row, seen from a different angle.

At the recommended settings the peak height is preserved to 0.2% (ratio 0.998), so the accuracy improvement is not being bought by flattening the fODF. The mean number of detected peaks there is 1.67 per voxel, below 2 because the 40 degree crossing in the condition set sits at the angular resolution limit of an Lmax 6 expansion.

`tau` shows the same pattern: 0.964 at 0.02, 0.998 at 0.1, 1.064 at 0.2, 1.213 at 0.4. The accuracy optimum at 0.1 coincides with the most faithful peak height.

### 4.2 Stage 2, `tau`

| `tau` | 0.02 | 0.05 | **0.1** | 0.2 | 0.4 |
|---|---|---|---|---|---|
| rel err fODF | 0.1204 | 0.1108 | **0.1052** | 0.1285 | 0.2231 |

Interior optimum at the pre-existing default of 0.1. Unchanged.

### 4.3 Stage 3, `Lmax_init`

| `Lmax_init` | 2 | **4** | 6 |
|---|---|---|---|
| rel err fODF | 0.1133 | **0.1052** | 0.1348 |

Interior optimum at the pre-existing default of 4. Unchanged. Initializing at the full order
(6) is worse than initializing at 4, which is the point of the low-order initialization in
Tournier et al.

### 4.4 Stage 4, Tikhonov matrix

| matrix | best `lambda_tikhonov` | rel err fODF |
|---|---|---|
| `identity` | 0.3 | 0.1052 |
| `laplacebeltrami` | 1 | **0.1026** |

Laplace-Beltrami is better by 2.5%, but only after re-optimizing its weight, which it needs:
it is normalized by `max(l(l+1))`, so at equal `lambda_tikhonov` it damps much more weakly
than the identity.

### 4.5 SNR dependence

Best point of the stage 1 grid, per SNR:

| SNR | `lambda_nonneg` | `lambda_tikhonov` | rel err fODF |
|---|---|---|---|
| 20 | 10 | 0 | 0.1358 |
| 30 | 10 | 0.3 | 0.1028 |
| 50 | 10 | 0.3 | 0.0767 |

The non-negativity weight is stable at 10 across the whole range. Only the Tikhonov weight
moves, and at the noisiest setting the sweep prefers none at all: with SNR 20 data the
non-negativity constraint alone already regularizes the problem enough, and adding isotropic
damping on top costs more in bias than it buys in variance.

---

## 5. Recommended vs original settings

### 5.1 End to end comparison

Same 270 voxels, all four scores:

| setting | rel err fODF | RMSE(plm) | negative mass | peak error |
|---|---|---|---|---|
| unregularized | 0.6148 | 0.3244 | 0.1589 | 12.04 deg |
| **original**: `lambda_nonneg=1`, `lambda_tikhonov=0.3`, identity | 0.1991 | 0.1101 | 0.0314 | 9.29 deg |
| **recommended**: `lambda_nonneg=10`, `lambda_tikhonov=0.3`, identity | 0.1052 | 0.0637 | 0.0044 | 8.16 deg |
| **best**: `lambda_nonneg=10`, `lambda_tikhonov=1`, Laplace-Beltrami | **0.1026** | 0.0622 | 0.0042 | 8.16 deg |

Against the previously suggested settings the recommended ones cut the fODF error roughly in
half (0.199 to 0.103, -48%), and the negative mass by a factor of 7 (0.031 to 0.004). For
reference, the negative mass of the *ground truth* fODF is 0.0034: the regularized estimate
is essentially as non-negative as the target, whereas the original settings were leaving
about 9x more spurious negative mass than the target has.

The peak angular error improves more modestly, 9.3 to 8.2 degrees. This is worth being clear
about: **the peak directions are the least sensitive of the four scores.** Most of the gain
is in the shape and the physicality of the fODF rather than in where its lobes point, and a
residual ~8 degrees remains at Lmax 6 regardless of regularization because the 40 degree
crossing in the condition set is near the angular resolution limit of a 6th order expansion.

### 5.2 Parameter by parameter

| parameter | original | recommended | changed? | evidence |
|---|---|---|---|---|
| `lambda_nonneg` | 1 | **10** | **yes** | 0.199 -> 0.105; interior optimum, turnover at 30 and 100 |
| `tau` | 0.1 | 0.1 | no | interior optimum confirmed (4.2) |
| `Lmax_init` | 4 | 4 | no | interior optimum confirmed (4.3) |
| `lambda_tikhonov` | 0 (off) | 0 (off) as a default; 0.3 with identity or 1 with LB when used | no | stays opt-in, see 5.3 |
| `TikhonovMatrix` | `identity` | `identity` | no | LB is 2.5% better but rescales the meaning of the weight, see 5.3 |
| `Ndirs` | 300 | 300 | no | not swept, cost parameter |
| `Niter` | 50 | 50 | no | not swept, converges in 4-6 |

### 5.3 What was changed in `SMI.m`, and what deliberately was not

**Changed:** the default of `lambda_nonneg`, from 1 to 10, in
`SMI.fODF_RegularizationDefaults`.

This only affects a user who switches the constraint on with `flag_nonneg = 1` and does not
choose a weight. That user previously got a setting measurably ~2x worse than achievable.

**Not changed, deliberately:**

- `lambda_tikhonov` stays 0 and `flag_nonneg` stays 0, so **both regularizers remain opt-in**.
  A fit that does not set `options.fODF_regularization` is bit-identical to before the sweep,
  verified: difference exactly 0. Turning a regularizer on by default would silently change
  every existing fODF result in the toolbox.
- `TikhonovMatrix` stays `identity` despite Laplace-Beltrami scoring 2.5% better. Switching it
  would silently rescale what any given `lambda_tikhonov` means, because the two matrices
  differ by the `max(l(l+1))` normalization: a user who had tuned `lambda_tikhonov = 0.3`
  against the identity would quietly get far weaker damping. A 2.5% gain does not justify
  that. It is documented in the README as the better choice when the weight is tuned with it.

### 5.4 Effect on the shipped example

`example_fODF_regularization.m` enables non-negativity without naming a weight, so it picks
up the new default. Its published numbers move accordingly (SNR 30, single condition):

| | rel err fODF, before | after |
|---|---|---|
| non-negativity | 0.1922 | **0.1041** |
| nonneg + tikhonov | 0.1875 | **0.1037** |

with negative mass falling from 0.0331 to 0.0044. The unregularized and Tikhonov-only rows
are unchanged, as they must be. The README figures were updated to match.

---

## 6. Limitations

Worth stating plainly, since these bound how far the numbers travel:

1. **One kernel.** All conditions use `f = 0.6, Da = 2, Depar = 2, Deperp = 0.5`. The optimum
   was not tested for sensitivity to the kernel itself, only to crossing angle and SNR.
2. **One protocol.** 3 shells x 64 directions at b = 1, 2, 3. A protocol with fewer
   directions or lower maximum b is more ill-conditioned and would likely want more
   regularization. Rerun the sweep with the protocol block edited rather than reusing these
   numbers.
3. **Two fibers only.** No single-fiber, three-fiber, or isotropic voxels in the condition
   set. A weight that suits crossings may over-smooth a single coherent bundle.
4. **Gaussian noise.** Real DWI magnitude data is Rician; at SNR 20 that difference is not
   negligible. The sweep does not include the Rician bias correction path.
5. **`smoothing = 0.05` fixes the sharpness of the ground truth.** A sharper target (less
   smoothing) is harder to represent at Lmax 6 and would shift the balance.
6. **Gaussian-noise-optimal is not tractography-optimal.** The primary score is an L2 error
   on the fODF. The peak error column shows the orientation gain is much smaller than the
   shape gain, so if what you care about is streamline direction, expect a smaller
   improvement than the headline number.

---

## 7. Reproducing

```matlab
example_fODF_regularization_sweep;
```

Prints every table in section 4, the per-SNR breakdown and a recommended settings block, and
produces four figures. Edit the protocol, kernel or condition block at the top to re-measure
for a different setup.

Figures are saved to `figures_fODF_sweep/` as PNG and EPS at 300 dpi, with lettered panels,
for direct use in a manuscript. Set `saveFigures = false` to display without writing files.

| figure | content |
|---|---|
| F1 `landscape` | the `lambda_nonneg` x `lambda_tikhonov` grid as four heat maps: relative fODF error, negative mass, peak angular error, peak amplitude ratio, optimum marked |
| F2 `peak_shrinkage` | peak amplitude ratio against each weight separately, and the accuracy vs shrinkage trade-off over the whole grid |
| F3 `secondary_params` | `tau`, `Lmax_init`, and the two Tikhonov matrices |
| F4 `snr` | how the optimum and the error curves move with SNR |

All numbers in this report were produced with GNU Octave 9 against the current `SMI.m`.
