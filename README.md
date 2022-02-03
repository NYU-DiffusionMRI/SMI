![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) `Still not ready to go online (THIS IS A DRAFT)`


# Standard Model Imaging (SMI) toolbox
This MATLAB toolbox contains all necessary functions for parameter estimation of the Standard Model of diffusion in white matter

## SMI input data
This implementation of the SM supports as input a 4D array of diffusion-weighted data (3D spatial arrangement of voxels and diffusion measurements along 4th dimention). Measurements can have varying:
- b-values (b)
- b-tensor shapes (β)
- echo times (TE)
 
Thus, each measurement is fully specified by: a b-value, a b-tensor shape, a unit direction (axis of symmetry of **B**), and TE. See the equation below to understand how these parameters make a b-tensor **B**

<img width="1206" alt=" AxSymB_wEq" src="https://user-images.githubusercontent.com/54751227/152437987-d79193d1-1ecc-4707-bdc3-f7cd2dec6ad6.png">

- If no β information is supplied the code assumes β=1 (linear tensor encoding measurements).
- If no TE is supplied, the code assumes that the TE is the same across all measurements.


## Recommended usage
We recommend using the default options but the code provides users with some flexibility.
Recommended inputs:
- Diffusion data (4D array) + protocol information + mask (binary 3D array) + noise map (3D array).
  - If a mask is not provided the fit will be performed in all voxels present in the first 3 dimensions of the 4D array.
  - If a noise map is not provided, it will be estimated using the repetitions of the non-diffusion-weighted images.
  - We recommend preprocessing your raw data with [DESIGNER](https://github.com/NYU-DiffusionMRI/DESIGNER) before doing SM estimation (DESIGNER outputs a robust noise estimation).

- The current SMI implementation is written in Matlab, future work may translate it to other languages.

- For technical details please look at the following publication:
  - Arxiv SM reproducibility link

## Example usage[^note] (THIS IS NOT READY
```
% Load data and protocol
nii=load_untouch_nii(fullfile(pathFiles,'dwi_preproc.nii.gz'));
dwi=abs(double(nii.img));
bval=load(fullfile(pathFiles,'dwi_preproc.bval'));
TE=load(fullfile(pathFiles,'dwi_preproc.TE'));
bshape=load(fullfile(pathFiles,'dwi_preproc.bshape'));
dirs=load(fullfile(pathFiles,'dwi_preproc.bvec'));
% Load mask
nii_mask=load_untouch_nii(fullfile(pathFiles,'mask.nii.gz'));
mask=logical(nii_mask.img);
% Load noise map and p-map
sigma_nii=load_untouch_nii(fullfile(pathFiles,'sigma_dki.nii'));
sigma=abs(sigma_nii.img);

% Keep only fixed TE data (TE=92ms)
keep_fixed_TE=TE==92;
dwi=dwi(:,:,:,keep_fixed_TE);
bval=bval(keep_fixed_TE);
dirs=dirs(:,keep_fixed_TE);
bshape=bshape(keep_fixed_TE);
TE=[];

% Perform Spherical Harmonics fit
[Slm,Sl,~,table_4D_sorted] = STARDOM_debug.Fit2D4D_LLS_RealSphHarm_wSorting_norm_varL(dwi,mask,bval,dirs,bshape,TE,Lmax);

% Perform PR training on Rotational Invariants and fitting SM kernel
KERNEL = STARDOM_debug.StandardModel_PR_fit_RotInvs(RotInvs,mask,sigma,bval,dirs,bshape,TE,lb_training,ub_training,Lmax_train,Ntraining,Nlevels,[0 0.4]);
```
## Advanced usage options
The code provides some additional flexibility:
- Batch processing (multiple datasets with identical protocols)
- Variable number of compartments ('IAS', 'EAS', 'FW', 'DOT')
- Rician bias correction
- Output spherical harmonic decomposition of the ODF for fiber tracking (normalized for MRtrix3)


## Theory: The Standard Model of diffusion in white matter
Multiple approaches to model the physics of water diffusion in white matter rely on similar assumptions. This led to the unifying framework dubbed Standard Model (SM) of diffusion in WM as formulated in ([Novikov et al., 2019](https://doi.org/10.1002/mrm.27101)).

<img width="1665" alt="kernel_wEqConvolution" src="https://user-images.githubusercontent.com/54751227/152442326-32f53e80-42c2-4f79-a1f6-eef7b1196844.png">

Briefly, axons (and possibly glial processes) are represented by impermeable zero-radius cylinders (the so-called “sticks”) arranged in locally coherent fiber fascicles. The diffusion in the extra-axonal space of each fascicle is assumed to be Gaussian and described by an axially symmetric diffusion tensor. The third, optional tissue compartment is the cerebro-spinal fluid (CSF). Such multicomponent fascicles (also called kernel) are distributed in a voxel according to an arbitrary fiber orientation distribution function (ODF). All fascicles in a voxel are assumed to have the same compartment fractions and diffusivities, and differ from each other only by orientation.

The SM encompasses a number of WM models made of anisotropic Gaussian compartments with axons represented by sticks ([Kroenke et al., 2004]( https://doi.org/10.1002/mrm.20260); [Jespersen et al., 2007](https://doi.org/10.1016/j.neuroimage.2006.10.037), [2010](https://doi.org/10.1016/j.neuroimage.2009.08.053); [Fieremans et al., 2011](https://doi.org/10.1016/j.neuroimage.2011.06.006); [Zhang et al., 2012](https://doi.org/10.1016/j.neuroimage.2012.03.072); [Sotiropoulos et al., 2012](https://doi.org/10.1016/j.neuroimage.2012.01.056); [Jensen et al., 2016](https://doi.org/10.1016/j.neuroimage.2015.09.049); [Jelescu et al., 2016a](https://doi.org/10.1002/nbm.3450); [Kaden et al., 2016](https://doi.org/10.1016/j.neuroimage.2016.06.002); [Reisert et al., 2017](https://doi.org/10.1016/j.neuroimage.2016.09.058); [Novikov et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.03.006); [Veraart et al., 2018](https://doi.org/10.1016/j.neuroimage.2017.09.030), to mention a few). From the SM point of view, earlier models impose constraints either on compartment parameters or the functional form of the fiber ODF; such constraints improve robustness but may introduce biases into the estimation of remaining parameters.

# SMI Authors
- Santiago Coelho
- Jelle Veraart
- Els Fieremans
- Dmitry Novikov

Do not hesitate to reach out to Santiago.Coelho@nyulangone.org (or [@santicoelho](https://twitter.com/santicoelho) in Twitter) for feedback, suggestions, questions, or comments[^note].

[^note]:
    Please cite these works if you use the SMI toolbox in your publication:
    - ARXIV COELHO et al 2022
    - Dmitry S. Novikov, Jelle Veraart, Ileana O. Jelescu, Els Fieremans, Rotationally-invariant mapping of scalar and orientational metrics of neuronal microstructure with diffusion MRI, NeuroImage, Volume 174, 2018, Pages 518-538.
    - Marco Reisert, Elias Kellner, Bibek Dhital, Jürgen Hennig, Valerij G. Kiselev, Disentangling micro from mesostructure by diffusion MRI: A Bayesian approach, NeuroImage, Volume 147, 2017, Pages 964-975.

