![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) `Still not ready to go online`


# Standard Model Imaging (SMI) toolbox
This MATLAB toolbox contains all necessary functions for parameter estimation of the Standard Model of diffusion in white matter

## The Standard Model of diffusion in white matter
Multiple approaches to model the physics of water diffusion in white matter rely on similar assumptions. This led to the unifying framework dubbed Standard Model (SM) of diffusion in WM as formulated in (Novikov et al., 2019).

<img width="1904" alt="SM_kernel_voxel_diagram_v2" src="https://user-images.githubusercontent.com/54751227/152257942-d8097c8b-8574-4b2e-8641-a614f9522edc.png">




Briefly, axons (and possibly glial processes) are represented by impermeable zero-radius cylinders (the so-called “sticks”) arranged in locally coherent fiber fascicles. The diffusion in the extra-axonal space of each fascicle is assumed to be Gaussian and described by an axially symmetric diffusion tensor. The third, optional tissue compartment is the cerebro-spinal fluid (CSF). Such multicomponent fascicles (also called kernel) are distributed in a voxel according to an arbitrary fiber orientation distribution function (ODF). All fascicles in a voxel are assumed to have the same compartment fractions and diffusivities, and differ from each other only by orientation.

The SM encompasses a number of WM models made of anisotropic Gaussian compartments with axons represented by sticks (Kroenke et al., 2004; Jespersen et al., 2007, 2010; Fieremans et al., 2011; Zhang et al., 2012; Sotiropoulos et al., 2012; Jensen et al., 2016; Jelescu et al., 2016a; Kaden et al., 2016; Reisert et al., 2017; Novikov et al., 2018; Ve- raart et al., 2018). From the SM point of view, earlier models impose constraints either on compartment parameters or the functional form of the fiber ODF; such constraints improve robustness but may introduce biases into the estimation of remaining parameters.


## SMI input data
This implementation of the SM supports as input a 4D array of diffusion data (3D spatial arrangement of voxels and diffusion measurements along 4th dimention). Measurements can have varying b-values (b), b-tensor shapes (β), anc echo times (TE). If no β information is supplied the code assumes β=1 (linear tensor encoding measurements). If no TE is supplied, the code assumes that the TE is the same across all measurements.

<img width="841" alt="AxSymB" src="https://user-images.githubusercontent.com/54751227/152258293-ac5827be-8bf9-4bf3-8948-7963c84778ff.png">


We recommend using the default options but the code is sufficiently flexible to provide users with some flexibility.

The current SMI implementation is written in Matlab, future work may translate it to other languages.

For technical details please look at the following publication:

- Arxiv SM reproducibility link

# Example usage
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


# Authors
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

