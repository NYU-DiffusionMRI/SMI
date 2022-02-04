![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) `Still not ready to go online (THIS IS A DRAFT)`


# Standard Model Imaging (SMI) toolbox
This MATLAB toolbox contains all necessary functions for parameter estimation of the Standard Model (SM) of diffusion in white matter. Check [our recent paper](https://arxiv.org/) for details on this implementation and on the Standard Model in general.

<br>

## SMI input data
This implementation of the SM supports as input a 4D array of diffusion-weighted data (3D spatial arrangement of voxels + diffusion measurements along 4th dimention). Measurements can have:
- Multiple b-values (b). This input can be a column or row vector with the same length as the 4th dimension of the input data. We recommend microstructural units [milliseconds / (squared micrometers)]. Note that b=1000 [seconds / (squared millimeters)] = 1 [milliseconds / (squared micrometers)].
- Multiple b-vectors or directions (each measurement must have its own direction, [3 x N] array).
- Multiple **B**-tensor shapes (β). This input can be a column or row vector with the same length as the 4th dimension of the input data. Only axially symmetric b-tensors are supported. β is a unitless scalar between -0.5 and 1 that indicates the **B**-tensor shape.
- Multiple echo times (TE, assumed to be in [milliseconds]). This input can be a column or row vector with the same length as the 4th dimension of the input data.

Each measurement is thus fully specified by: a b-value (b), a unit direction (**u**) (axis of symmetry of **B**), a b-tensor shape (β), and TE. See the figure and equation below to understand how these parameters make a b-tensor **B**:
<p align="center">
  <img width="550" alt=" AxSymB_wEq" src="https://user-images.githubusercontent.com/54751227/152437987-d79193d1-1ecc-4707-bdc3-f7cd2dec6ad6.png">
</p>

  - b-values and directions must be supplied for each measurement.
  - If no β is supplied (this input can be empty array) the code assumes β=1 (linear tensor encoding, LTE, measurements).
  - If no TE is supplied (this input can be empty array), the code assumes that the TE is the same across all measurements. In this case compartmental T2 values will not be outputted.

<br>


## SMI outputs
Standard Model parameters (diffusion) + compartmental T2 values (only if multiple TE data was provided). See figure below for some examples:
<img width="1690" alt="SM_maps_github" src="https://user-images.githubusercontent.com/54751227/152456997-5f24f886-03f9-4eb2-a5f8-1dd767134eae.png">
- Fractions (f, fw) and anisotropy (p2) are adimensional.
- Diffusivities are returned in microstructural units [squared micrometers / milliseconds].
- Compartmental T2 values are returned in [milliseconds].
- Note that compartmental T2 maps will only be outputed if variable TE data was used.


<br>

## Recommended usage
We recommend using the default options but the code provides users with some flexibility.
Recommended inputs:
- Diffusion data (4D array) + protocol information + mask (binary 3D array) + noise map (3D array).
  - If a mask is not provided the fit will be performed in all voxels present in the first 3 dimensions of the 4D array.
  - If a noise map is not provided, it will be estimated using the repetitions of the non-diffusion-weighted images.
  - We recommend preprocessing your raw data with [DESIGNER](https://github.com/NYU-DiffusionMRI/DESIGNER) before doing SM estimation (DESIGNER outputs a robust noise estimation).

- The current SMI implementation is written in Matlab, future work may translate it to other languages.

- For technical details please look at our recent publication: [Reproducibility of the Standard Model of diffusion in white matter on clinical MRI systems, (2022), ArXiv](https://arxiv.org/).

<br>

## Example usage[^note]
We provide an example dataset [here](https://drive.google.com/drive/folders/1TQzZGM7PdTf1kplwfWLIIRn8q0kE3Nix?usp=sharing) that was preprocessed with [DESIGNER](https://github.com/NYU-DiffusionMRI/DESIGNER). This contains multiple b-values, tensor shapes, and echo times. See the example file 'example.m' with the SMI fit of a subset of the data with only 1 TE, and with the full dataset. 

```
% Load data and protocol
nii=load_untouch_nii(fullfile(pathFiles,'dwi_preproc.nii.gz'));
dwi=abs(double(nii.img));
bval=load(fullfile(pathFiles,'dwi_preproc.bval'));
TE=load(fullfile(pathFiles,'dwi_preproc.TE'));
bshape=load(fullfile(pathFiles,'dwi_preproc.bshape'));
dirs=load(fullfile(pathFiles,'dwi_preproc.bvec'));
% Load mask

% Perform Spherical Harmonics fit
[Slm,Sl,~,table_4D_sorted] = STARDOM_debug.Fit2D4D_LLS_RealSphHarm_wSorting_norm_varL(dwi,mask,bval,dirs,bshape,TE,Lmax);

% Perform PR training on Rotational Invariants and fitting SM kernel
KERNEL = STARDOM_debug.StandardModel_PR_fit_RotInvs(RotInvs,mask,sigma,bval,dirs,bshape,TE,lb_training,ub_training,Lmax_train,Ntraining,Nlevels,[0 0.4]);
```
## Some advanced usage options
The code provides some additional flexibility:

- Batch processing. Parameter estimation for multiple datasets with identical protocols. Here the machine learning training is done only once, regression coefficients are stored and applied to all.
- Variable number of compartments: 'IAS', 'EAS', 'FW', 'DOT'. Any combination of these is allowed. Only these maps will be outputted.
- User defined parameter distributions for the training data (for the machine learning estimator that does RotInvs -> kernel).
- Rician bias correction (to de-bias the DWI before the spherical harmonics fit).
- **(NOT READY YET)** Output spherical harmonic decomposition of the ODF for fiber tracking (normalized for using it with [MRtrix3](https://mrtrix.readthedocs.io/en/0.3.16/workflows/global_tractography.html)).

<br>

## Some basic theory: The Standard Model of diffusion in white matter
Multiple approaches to model the physics of water diffusion in white matter rely on similar assumptions. This led to the unifying framework dubbed Standard Model (SM) of diffusion in WM as formulated in ([Novikov et al., 2019](https://doi.org/10.1002/mrm.27101)).

<img width="1657" alt="kernel_wEqConvolution_v2" src="https://user-images.githubusercontent.com/54751227/152564788-fc6a0fe0-1002-4354-b75e-3f962303a9ad.png">

Briefly, axons (and possibly glial processes) are represented by impermeable zero-radius cylinders (the so-called “sticks”) arranged in locally coherent fiber fascicles. The diffusion in the extra-axonal space of each fascicle is assumed to be Gaussian and described by an axially symmetric diffusion tensor. The third, optional tissue compartment is the cerebro-spinal fluid (CSF). Such multicomponent fascicles (also called kernel) are distributed in a voxel according to an arbitrary fiber orientation distribution function (ODF). All fascicles in a voxel are assumed to have the same compartment fractions and diffusivities, and differ from each other only by orientation.

The SM encompasses a number of WM models made of anisotropic Gaussian compartments with axons represented by sticks ([Kroenke et al., 2004]( https://doi.org/10.1002/mrm.20260); [Jespersen et al., 2007](https://doi.org/10.1016/j.neuroimage.2006.10.037), [2010](https://doi.org/10.1016/j.neuroimage.2009.08.053); [Fieremans et al., 2011](https://doi.org/10.1016/j.neuroimage.2011.06.006); [Zhang et al., 2012](https://doi.org/10.1016/j.neuroimage.2012.03.072); [Sotiropoulos et al., 2012](https://doi.org/10.1016/j.neuroimage.2012.01.056); [Jensen et al., 2016](https://doi.org/10.1016/j.neuroimage.2015.09.049); [Jelescu et al., 2016a](https://doi.org/10.1002/nbm.3450); [Kaden et al., 2016](https://doi.org/10.1016/j.neuroimage.2016.06.002); [Reisert et al., 2017](https://doi.org/10.1016/j.neuroimage.2016.09.058); [Novikov et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.03.006); [Veraart et al., 2018](https://doi.org/10.1016/j.neuroimage.2017.09.030), to mention a few). From the SM point of view, earlier models impose constraints either on compartment parameters or the functional form of the fiber ODF; such constraints improve robustness but may introduce biases into the estimation of remaining parameters.

<br>

# SMI Authors
- Santiago Coelho
- Jelle Veraart
- Els Fieremans
- Dmitry Novikov

Do not hesitate to reach out to Santiago.Coelho@nyulangone.org (or [@santicoelho](https://twitter.com/santicoelho) in Twitter) for feedback, suggestions, questions, or comments[^note].

[^note]:
    Please cite these works if you use the SMI toolbox in your publication:
    - ARXIV COELHO et al 2022.
    - Dmitry S. Novikov, Jelle Veraart, Ileana O. Jelescu, Els Fieremans, Rotationally-invariant mapping of scalar and orientational metrics of neuronal microstructure with diffusion MRI, NeuroImage, Volume 174, 2018, Pages 518-538.
    - Marco Reisert, Elias Kellner, Bibek Dhital, Jürgen Hennig, Valerij G. Kiselev, Disentangling micro from mesostructure by diffusion MRI: A Bayesian approach, NeuroImage, Volume 147, 2017, Pages 964-975.


## LICENSE

A [US patent](https://patents.google.com/patent/US20160343129A1/en) contains some of the related developments. 

```
%  Authors: Santiago Coelho (santiago.coelho@nyulangone.org), Jelle Veraart, Els Fieremans, Dmitry Novikov
%  Copyright (c) 2022 New York University
%              
%   Permission is hereby granted, free of charge, to any non-commercial entity ('Recipient') obtaining a 
%   copy of this software and associated documentation files (the 'Software'), to the Software solely for
%   non-commercial research, including the rights to use, copy and modify the Software, subject to the 
%   following conditions: 
% 
%     1. The above copyright notice and this permission notice shall be included by Recipient in all copies
%     or substantial portions of the Software. 
% 
%     2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
%     NOT LIMITED TO THE WARRANTIESOF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
%     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE
%     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
% 
%     3. In no event shall NYU be liable for direct, indirect, special, incidental or consequential damages
%     in connection with the Software. Recipient will defend, indemnify and hold NYU harmless from any 
%     claims or liability resulting from the use of the Software by recipient. 
% 
%     4. Neither anything contained herein nor the delivery of the Software to recipient shall be deemed to
%     grant the Recipient any right or licenses under any patents or patent application owned by NYU. 
% 
%     5. The Software may only be used for non-commercial research and may not be used for clinical care. 
% 
%     6. Any publication by Recipient of research involving the Software shall cite the references listed
%     below.
%
% REFERENCES:
% - Coelho, S., Baete, S., Lemberskiy, G., Ades-Aron, B., Barrol, G., Veraart, J., Novikov, D.S., Fieremans, E., 2022. Reproducibility of the Standard Model of diffusion in white matter on clinical MRI systems, ArXiv
% - Novikov, D.S., Veraart, J., Jelescu, I.O., Fieremans, E., 2018. Rotationally-invariant mapping of scalar and orientational metrics of neuronal microstructure with diffusion MRI. NeuroImage 174, 518 – 538
% - Reisert, M., Kellner, E., Dhital, B., Hennig, J., Kiselev, V.G., 2017. Disentangling micro from mesostructure by diffusion MRI: A Bayesian approach. NeuroImage 147, 964 – 975.
```
