# README for Claude

Handoff for the next agent working on fODF modulation for tractography through
edema. Written after the modulation work stalled, deliberately structured so
that the stall is visible and the dead ends are not re-walked.

**Read section 7 before proposing anything.** Most of the obvious ideas have
been tried and measured, and the measurements are counter-intuitive in several
places.

---

## 1. The goal, stated precisely

Track fibres **through peritumoral edema** using SMI fODFs in MRtrix.

The operational problem: find a per-voxel scalar weight `w`, computable from
`SMI.fit` output, such that after `fODF -> w*fODF`

- **all** white matter keeps amplitude above the tractography cutoff, including
  fibre crossings and including white matter inside edema, and
- grey matter and CSF fall below it.

The edema constraint is the hard one and is non-negotiable: **the user has
rejected tissue-type criteria three times**, because anything keyed to "is this
healthy white matter" deletes exactly the region of interest. Any proposal must
be tested against a voxel with collapsed axonal fraction and high free water but
coherent fibres. If it suppresses that, it is dead on arrival.

---

## 2. Status snapshot

| | |
|---|---|
| Repo | `iprentis777/SMI` |
| Branch | `claude/fodf-modulation-edema-38eh81` (restarted from master after PR #3 merged) |
| Merged | PR [#3](https://github.com/iprentis777/SMI/pull/3) as `e0f6314` |
| In the repo | `SMI.fODF_ModulationWeight`, `SMI.modulate_fODF`, `SMI.fODF_ModulationDefaults`, `grab_pl`, `grab_kernel_pl`, the opt-in `options.fODF_modulation` flag, `example_fODF_modulation.m`, `fODF_modulation_helpers.m`, `REPORT_fODF_modulation.md`, `0004fODFmodulation.patch` |
| Not in the repo | everything from section 6 onward: the crossing-fibre findings, anisotropic power, the noise correction. All scratch. |

**The shipped default `p2product` is now known to be harmful in crossings.**
It is still the default in `SMI.m`. It has not been changed because nothing has
yet been proven better across the board. Do not treat it as endorsed; treat it
as the thing that motivated everything in section 7.

Modulation is opt-in (`flag_modulate` defaults to 0) and `out.plm`, `out.pl`,
`out.kernel` are bit-identical whether it is on or off, so nothing shipped is
currently doing damage to a default fit.

---

## 3. Conventions that must not be got wrong

Pinned analytically and verified numerically (round-trip `max|err| = 1.2e-15`).

`out.plm` is in the **normalized** convention `p_00 = 1`, covering `l = 2..Lmax`
only (the `l=0` term is not stored).

```
SH coefficients :  f_lm = plm .* sqrt((2l+1)/(4*pi))
l = 0 term      :  f_00 = 1/sqrt(4*pi)        (constant in every voxel)
fODF amplitude  :  A(u) = 1/(4*pi) + sum_{l>=2,m} f_lm Y_lm(u)
Lmax from count :  Lmax = sqrt(2*Nlm + 9/4) - 3/2
```

Consequences that keep mattering:

1. **The fODF integrates to 1 in every voxel.** It carries no density
   information at all. A CSF voxel and a coherent WM voxel have equal mass.
2. **The isotropic floor is a fixed `1/(4*pi) = 0.0796`**, which is *above*
   MRtrix's default `iFOD2 -cutoff 0.05`. An unmodulated SMI fODF therefore
   passes the tractography termination test **in every voxel of the brain**,
   CSF included. This is the original motivation for the whole effort.
3. To rescale a peak to target `T`, scale the anisotropic part by
   `s = (T - 1/(4pi))/(peak - 1/(4pi))`, never `T/peak`.
4. Scaling all `l>0` coefficients uniformly preserves peak orientation exactly.
   Clipping SH coefficients **independently** does not and can swing peak
   orientations. Never do the latter.

**Two different p2 exist and they are not the same quantity.** Same for p4.

- `out.kernel(:,:,:,ip2)` — fitted jointly with the kernel by polynomial
  regression on rotational invariants. `ip2 = 6` without T2 fitting, `8` with.
  The map count alone is ambiguous; `SMI.grab_kernel_pl` resolves it from
  `out.shells` row 4 (TE).
- `out.pl(:,:,:,1)` — the l=2 invariant of the **deconvolved** fODF, a raw
  Euclidean norm, unclipped.

With `Lmax = [0 8 8 8]`, `RotInv_Lmax` defaults to 4, so **`out.kernel` has 7
maps, not 6**: `[f Da Depar Deperp fw p2 p4]`.

---

## 4. Environment and validation

Linux container, **GNU Octave 9, no MATLAB**. Working setup:

```
apt-get update && apt-get install -y octave octave-statistics octave-image
pkg load statistics; pkg load image;
```

`SMI.fit` **does run end-to-end in Octave** with three shims (they live in
`scratchpad/mod/stubs/`, recreate if the container was recycled):

| shim | why |
|---|---|
| `round(x,n)` | MATLAB two-arg form, used in `Group_dwi_in_shells_b_beta_TE` |
| `discretize(x,edges)` | not implemented in Octave, used in `StandardModel_MLfit_RotInvs` |
| `datetime()` | not implemented, used only by the log writer |

Other traps:

- **Octave cannot call functions defined at the end of a script**, MATLAB
  requires them there. Incompatible. That is why `fODF_modulation_helpers.m`
  and `fODF_regularization_score.m` are separate files. Subfunctions of a
  *function* file work in both.
- **`SMI.vectorize` takes a different branch if any spatial dimension is a
  singleton** (`ismatrix(S)` is true for a collapsed 3D array). Always build
  simulation volumes with all three dimensions > 1.
- Graphics are unusable (broken text renderer). Stub `figure/histogram/
  subplot/xlabel/ylabel/title/box/grid/drawnow` to test plotting code paths.
- `strel('disk',r,0)` — Octave requires the explicit `N`; `N=0` is an exact
  disk and is the better choice in MATLAB too.

**Parse-checking is not enough.** Two runtime errors reached the user that way,
because the interesting code sits inside `try` blocks. The pattern that works:
extract the post-fit section of the user's script *from the file itself*, stub
the NIfTI I/O, feed a synthetic `out`, and run all flag combinations. See
`scratchpad/script/harness.m`.

---

## 5. The simulation that matters

Forward model is the exact construction `SMI.get_plm_from_S_and_kernel`
inverts, so it round-trips to 1.2e-15:

```
S(u_i)/S0 = sum_{lm} K_l(b_i) * p_lm * Y_lm(u_i) * sqrt((2l+1)*4pi)
```

Eight classes, 80 realisations each, HCP protocol (b = 0,1,2,3 ms/um^2,
18/90/90/90 dirs, `Lmax = [0 8 8 8]`), Rician noise, fitted through the real
`SMI.fit` at SNR 30 and 15.

| class | kernel `[f Da Depar Deperp fw]` | fODF |
|---|---|---|
| 1 WM 1 fibre | `0.60 2.0 2.0 0.50 0.05` | Watson κ=16 |
| 2 WM 2 fibre 60° | same | two Watsons |
| 3 **WM 3 fibre orthogonal** | same | three orthogonal Watsons |
| 4 WM 3 fibre realistic | same | unequal weights, non-orthogonal |
| 5 WM edema 1 fibre | `0.30 2.0 2.2 1.00 0.50` | Watson κ=12 |
| 6 **WM edema 2 fibre** | same | two Watsons |
| 7 GM | `0.15 1.2 1.0 0.80 0.10` | Watson κ=0.8 |
| 8 CSF | `0.02 2.0 3.0 3.00 0.95` | isotropic |

Classes 1-6 must all survive; 7-8 must be suppressed. **Class 3 is the
algebraic worst case, class 6 the clinically hardest.** Scripts:
`scratchpad/mod/{ap_sim,ap_score,sig_sim,sig_score,diag_crossing,ap_truth}.m`.

Primary score: the fraction of each class left **below** MRtrix's default
cutoff of 0.05. Secondary: an adaptive cutoff retaining 95% of all WM, then the
false-positive rate on GM+CSF — scale-free, so it does not reward a weight for
simply shrinking everything.

---

## 6. Ground truth: what is actually there to detect

No fitting, no noise. This is the ceiling for any method.

| fODF | p2 | p4 | p6 | p8 | AP |
|---|---|---|---|---|---|
| 1 fibre | 0.903 | 0.713 | 0.495 | 0.304 | 0.552 |
| 2 fibre 60° | 0.597 | 0.425 | 0.402 | 0.207 | 0.376 |
| **3 fibre orthogonal** | **0.000** | **0.544** | 0.175 | 0.218 | 0.297 |
| 3 fibre realistic | 0.436 | 0.541 | 0.225 | 0.181 | 0.330 |
| GM (κ=0.8) | 0.114 | 0.009 | 0.000 | 0.000 | 0.039 |
| isotropic | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

**`p2` is identically zero for an orthogonal three-way crossing.** Not small —
zero, at any SNR, with any regularizer. The second-moment tensor of three
equal orthogonal fibre populations is isotropic. The l=4 invariant is not
(cubic symmetry has a non-zero l=4 term), so `p4 = 0.544` carries the entire
crossing.

`AP` = total anisotropic power, normalised so a delta fibre gives 1:

```
AP = sqrt( sum_{l>=2} (2l+1) p_l^2 / sum_{l>=2} (2l+1) )
```

---

## 7. Dead ends, with evidence

### 7.1 `p2product` (the shipped default) destroys crossings

Fraction below the MRtrix cutoff, SNR 30, Lmax 8, regularized:

| class | none | **p2product** | AP l=2,4 | APnc l=2,4 |
|---|---|---|---|---|
| WM 1 fibre | 0.0% | 0.0% | 0.0% | 0.0% |
| WM 2 fibre 60° | 0.0% | 0.0% | 0.0% | 0.0% |
| **WM 3 fibre orthogonal** | 0.0% | **100.0%** | 0.0% | 0.0% |
| WM 3 fibre realistic | 0.0% | 0.0% | 0.0% | 0.0% |
| WM edema 1/2 fibre | 0.0% | 0.0% | 0.0% | 0.0% |
| GM | 0.0% | 100.0% | 65.0% | 96.2% |
| CSF | 0.0% | 100.0% | 33.8% | 100.0% |

**Nuance that matters and is easy to miss:** `p2product` only *catastrophically*
fails the perfectly orthogonal equal-weight case. The realistic three-fibre
class survives the cutoff — but its weight is 0.166 against 0.679 for a single
fibre, a **4x dimming**. That is consistent with the user's mrview screenshot,
where the corona radiata is dim but present, not absent. If a future diagnosis
of "missing fibres" does not match a 4x dimming, look for a second cause.

### 7.2 The kernel p4 is structurally blind to crossings

`SMI.Get_uniformly_distributed_SM_prior` draws `p4 = rand * p2 * 0.9`, so the
training prior enforces `p4 <= 0.9*p2`. With `p2 = 0` the polynomial regression
**cannot represent** `p4 = 0.54`; it is off the training manifold. Measured
kernel-based anisotropic power for the orthogonal crossing: **0.000**.

Any kernel-derived anisotropy is therefore dead for crossings. Only the
deconvolved `out.pl` sees them. Changing this would mean retraining the prior.

### 7.3 Raw anisotropic power fails — noise does not cancel

`p_l` is a norm of `(2l+1)` noisy coefficients, so it is positive when the truth
is zero. Measured `pl4` in CSF: **0.196** at SNR 30, **0.349** at SNR 15,
against a true 0.00. AP sums squares, so noise accumulates. At SNR 15 raw AP is
**inverted**: CSF 0.207 > orthogonal crossing 0.170. FP GM+CSF 78.8%,
separation 0.28 — worse than no weight at all.

Higher bands are worse (noisiest, and weighted by `(2l+1)`). Including l=6,8
drops 22-30% of the orthogonal crossing below cutoff. **Restrict to l = 2,4.**

### 7.4 The signal-domain (sigma-free) route does not work

The idea: `p_l = S_l/(K_l S_0)`, and dividing by a collapsing `K_l` is what
amplifies noise. `S_l/S_0` is measured directly, needs no sigma, and should be
doubly suppressed in CSF (both `p_4` and `K_4` collapse) but only singly in a
crossing. `out.RotInvs.S0/S2/S4` are already stored.

**Measured: it does not separate.** Orthogonal crossing 0.0447 vs CSF 0.0365 —
1.2x at SNR 30, 1.07x at SNR 15. FP 13.1% / 46.9%, separation 0.81 / 0.46.

Why the argument failed: `K_4` attenuates the l=4 signal hard *even in healthy
white matter* (`a4 = 0.080` in the crossing vs 0.023 in CSF noise). Dividing by
`K_4` amplifies signal and noise equally. **The information content of the l=4
band is fixed; changing domain does not create any.** This generalises — do not
expect any re-parameterisation of the same measurement to escape the noise
floor.

`Nnegative_dirs` as a noise probe: also marginal (FP 9.4% / 36.2%).

### 7.5 Tissue-fraction weights fail edema, as expected

Adaptive cutoff retaining 95% of WM, SNR 15:

| weight | FP GM+CSF | edema retained |
|---|---|---|
| `f` | 0.0% | **94%** |
| `1-fw` | 31.5% | **85%** |
| any p2-based | 0.0-100% | **100%** |

`1-fw` additionally fails at its own job, leaving 63% of grey matter above
threshold. This is the third confirmation; stop proposing tissue fractions.

### 7.6 Things that turned out not to be problems

- **Tikhonov does not destroy the l>=4 bands.** AP for the crossing: 0.163 at
  `lambda_tikhonov = 0.3` vs 0.164 at 0. Feared twice, wrong twice.
- **Regularization is load-bearing for blow-ups, and modulation cannot replace
  it.** At SNR 15 unregularized, CSF fODF peaks reached **2.2e14**, reproducing
  the 1e13 amplitudes seen in real data; the regularized fit caps the same
  voxels at 0.449. The mechanism is the deconvolution gain `g_2 = 1/||K_2||`,
  measured 3.6 in CSF vs 0.6 in WM. Since the weight is clipped at 1,
  `w * 1e13` is still `1e13`.
- **The default `degenerate = 'clip'` is a bad default.** In a blown-up voxel
  `pl2 >> 1`, so the raw value exceeds the clip and the voxel is assigned weight
  **exactly 1.0 — the maximum in the volume**, while healthy voxels get 0.3-0.7.
  The pathological voxels are the *least* attenuated. `'reject'` (weight 0)
  exists and is arguably the right default. Not yet changed.

---

## 8. The one thing that currently works, and why it is not shipped

**Noise-floor-subtracted anisotropic power over l = 2,4.**

`p_lm` is recovered by dividing the l-th signal harmonic by `K_l`, so with
roughly uniform directions

```
Var(p_lm) ~ sigma^2 / ((2l+1) ||K_l||^2)
E[p_l^2 | noise only] ~ (sigma * g_l)^2,   g_l = 1/||K_l(b)||

w = sqrt( ( 5*max(p2^2 - (sigma g_2)^2, 0)
          + 9*max(p4^2 - (sigma g_4)^2, 0) ) / 14 )
```

Both `sigma` and the kernel are already in `out`. Measured: **0.0% FP on GM+CSF
at both SNRs, all six WM classes fully retained, separation 2.03 / 2.06**,
against 0.00 for `p2product` and 0.48-0.81 for no weight. It is the only
candidate that has ever satisfied every constraint simultaneously.

**The sigma dependence is less fragile than it looks, and the failure is
one-sided:**

| sigma used | SNR 30 FP / sep | SNR 15 FP / sep |
|---|---|---|
| correct | 0.0% / 2.03 | 0.0% / 2.06 |
| underestimated 2x | 0.6% / 1.37 | **10.6% / 0.57** |
| overestimated 1.5x | 0.0% / 4.85 | **0.0% / 7.52** |

Over-estimating sigma *improves* it. So an accurate sigma is not required — an
**upper bound** is, and erring high is nearly free. This matters for HCP
specifically: eddy/topup interpolation spatially correlates noise, which makes a
b=0-repetition estimate read **low**, the dangerous direction. A deliberate
safety factor is well motivated rather than arbitrary.

**Why it is not shipped:** simulation only; sigma sensitivity not re-tested at
Lmax 6; the safety factor is a free parameter nobody has justified; and the
whole thing has never touched real data.

### Independent of any weighting: Lmax 8 is hurting

| SNR 30, regularized | Lmax 8 | Lmax 6 |
|---|---|---|
| CSF peak | 0.352 | 0.271 |
| 3-fibre orthogonal peak | 0.343 | 0.294 |
| crossing vs CSF | **inverted** | correct |
| unweighted FP GM+CSF | 52.5% | **21.9%** |

At Lmax 8 a genuine three-way crossing is **dimmer than CSF noise** before any
weighting — no shape-derived weight can fix that. The user runs
`Lmax = [0 8 8 8]` with 90 directions per shell, and the regularization defaults
were tuned at Lmax 6. Dropping to `[0 6 6 6]` may be the cheapest single
improvement available and is untested on real data.

---

## 9. Where this stalled, stated sharply

Every approach so far computes a **per-voxel scalar from the fitted model**, and
they all run into the same wall:

- **Anything derived from fODF shape** (p2, p4, AP, max p_l) cannot separate
  "isotropic because CSF" from "isotropic because fibres cross" — both genuinely
  are spread out on the sphere. p2 is the extreme case, being exactly zero for a
  symmetric crossing.
- **Anything derived from the kernel** (f, 1-fw, kernel p2/p4) separates
  crossings from CSF fine, but fails edema, because edema genuinely does have low
  `f` and high free water. That is physics, not a measurement artefact.
- **Anything that fights the noise floor** needs a noise estimate. Changing
  domain (section 7.4) does not escape it.

The information needed to separate "crossing white matter" from "CSF" is
**not present in a single voxel's fitted model** with enough margin at Lmax 8.
That is the stall.

---

## 10. Untried directions worth brainstorming

None of these have been tested. They are listed because they break one of the
assumptions above, not because they are known to work.

1. **Spatial context.** *Every approach so far is strictly per-voxel.* Noise is
   spatially incoherent; tracts are not. A voxel whose principal peak agrees
   with its neighbours' is structurally different from one whose peak is noise,
   and this is completely independent of amplitude, anisotropy, and sigma. This
   is the biggest unexplored axis and probably the first thing to try.
   Caution: `SMI.fit` currently couples no voxels, so this must be post-hoc,
   and it interacts with the pre-fit mask erosion already in the driver scripts.

2. **Stop modulating amplitude; constrain geometry instead.** MRtrix has `-mask`,
   `-seed_image`, `-include/-exclude` and ACT (`-act` with a 5TT image). Tumour
   tractography commonly uses a lesion-aware 5TT. This decouples "where may
   tracks go" from "what does fODF amplitude say" and sidesteps the entire
   problem. The reason it was not pursued is that a 5TT segmentation is a
   tissue-type criterion — but applied as a *tracking mask* rather than a
   per-voxel weight it may be acceptable where a weight was not. Worth asking
   the user directly.

3. **Subtract the isotropic floor at export.** `A(u) - 1/(4pi)`, clipped at 0,
   makes CSF ~0 and restores meaning to MRtrix's default cutoff with no weight
   and no sigma. **But note:** it is a monotone shift, so it cannot change the
   *ordering* of voxels and adds no discrimination — it only fixes the
   "everything passes 0.05" problem structurally. Cheap, honest, partial.

4. **Multi-tissue deconvolution.** MSMT-CSD solves exactly this by giving the
   isotropic signal its own compartment rather than forcing it into the fODF.
   SMI already estimates `fw`. Whether an SMI-flavoured 3-tissue deconvolution
   is tractable is unknown, and the edema objection to `fw` would need
   re-examining in that framing.

5. **Peak-based rather than amplitude-based export.** Extract fixels/peaks and
   threshold on peak count, sharpness, or separation angle rather than
   modulating the continuous fODF. `fODF_peak_amplitudes` in the driver script
   already does most of the extraction.

6. **Is the corona radiata really the problem?** The perfectly orthogonal
   equal-weight crossing is a measure-zero degenerate case; the realistic
   three-fibre class survived `p2product`. Before optimising for it, confirm on
   the user's data that the dim regions really are 3-way crossings and that the
   dimming is the ~4x that `p2product` predicts. If it is much more than 4x,
   there is a second cause not yet identified.

---

## 11. Open questions that need real data, not simulation

- Does `p2` actually hold up in the user's edematous voxels? The whole approach
  assumes displaced fibres stay coherent. **Never verified.** Concrete first
  experiment: map `p2`, `p4` and `f` in a known edema ROI and see which degrades.
- **Is SMI's SH basis MRtrix's?** `out.fODF_modulated` is in SMI's own basis
  (`SMI.get_even_SH` with `out.CS_phase`). The ordering and Condon-Shortley
  convention have **never been checked** against MRtrix. Everything downstream
  depends on this and it is a half-hour job.
- Are the stray giant glyphs inside `mask3D`? If the weight map reads 1.0 there
  it is the clip-to-max path of 7.6; if 0, they are outside the mask and
  something else is rendering them.
- Does SMI's internal b=0 sigma estimate suffice, given HCP interpolation?

---

## 12. Working preferences observed

- Wants to be told when a suggestion is wrong, **with evidence**. Has pushed
  back correctly four times: free-water filtering, the percentile-vs-outlier
  argument, p2/p4 as the modulator, and preferring l=2,4 over high orders.
- Prefers measured over assumed. The regularization sweep exists because a
  documented "around 1 is reasonable" turned out to be 2x off; the crossing
  finding exists because a simulation lacked the case that mattered.
- **Nothing goes in the repo until it is proven.** Stated explicitly. Honour it.
- Asks where changes were made — point at file and line numbers.
- Figures matter; a manuscript is the eventual target. A manuscript-grade
  simulation script is wanted but blocked: the current simulation has no
  ground-truth *density*, since every voxel is a single kernel. It needs a
  second tissue type before density-weighted scoring means anything.
- Two driver scripts live only on the user's machine and are sent as files:
  `run_smi_hcp.m` (single subject, now carries `useFODFmod`) and
  `run_smi_batch.m` (5-subject edema batch, supplies sigma from
  `level_noise.nii.gz`, has fODF peak truncation that `run_smi_hcp.m` lacks).
  Peak truncation must run **before** modulation — its thresholds are absolute
  amplitudes.
