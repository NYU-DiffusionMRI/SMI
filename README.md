# Standard Model Imaging (SMI) toolbox
[This MATLAB toolbox](https://github.com/NYU-DiffusionMRI/SMI/blob/master/SMI.m) contains all necessary functions for parameter estimation of the Standard Model (SM) of diffusion in white matter[^note]. Check [our recent paper](https://arxiv.org/pdf/2202.02399.pdf) for details on this implementation and on the Standard Model in general. Below we provide instructions on how to run the toolbox. See the '[example.m](https://github.com/NYU-DiffusionMRI/SMI/blob/master/example.m)' script that performs the parameter estimation in an [example dataset](https://cai2r.net/resources/standard-model-of-diffusion-in-white-matter-the-smi-toolbox/).

<br>

## Overview: The Standard Model of diffusion in white matter

Over the last 15-20 years, multiple approaches aimed to model the physics of water diffusion in white matter have relied on similar assumptions. This led to the unifying framework dubbed Standard Model (SM) of diffusion in WM as formulated in ([Novikov et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.03.006))[^2]. In a nutshell, this model disentangles signal contributions from different structures, i.e. compartments, present in a white matter voxel. 

<img width="1657" alt="kernel_wEqConvolution_v2" src="https://user-images.githubusercontent.com/54751227/152564788-fc6a0fe0-1002-4354-b75e-3f962303a9ad.png">

Briefly, axons (and possibly glial processes) are represented by impermeable zero-radius cylinders (the so-called “sticks”) arranged in locally coherent fiber fascicles. The diffusion in the extra-axonal space of each fascicle is assumed to be Gaussian and described by an axially symmetric diffusion tensor. The third, optional tissue compartment is the cerebro-spinal fluid. Such multicomponent fascicles (also called kernel) are distributed in a voxel according to an arbitrary fiber orientation distribution function (ODF). All fascicles in a voxel are assumed to have the same compartment fractions and diffusivities, and differ from each other only by orientation.

The SM encompasses a number of WM models made of anisotropic Gaussian compartments with axons represented by sticks ([Kroenke et al., 2004]( https://doi.org/10.1002/mrm.20260); [Jespersen et al., 2007](https://doi.org/10.1016/j.neuroimage.2006.10.037), [2010](https://doi.org/10.1016/j.neuroimage.2009.08.053); [Fieremans et al., 2011](https://doi.org/10.1016/j.neuroimage.2011.06.006); [Zhang et al., 2012](https://doi.org/10.1016/j.neuroimage.2012.03.072); [Sotiropoulos et al., 2012](https://doi.org/10.1016/j.neuroimage.2012.01.056); [Jensen et al., 2016](https://doi.org/10.1016/j.neuroimage.2015.09.049); [Jelescu et al., 2016a](https://doi.org/10.1002/nbm.3450); [Kaden et al., 2016](https://doi.org/10.1016/j.neuroimage.2016.06.002); [Reisert et al., 2017](https://doi.org/10.1016/j.neuroimage.2016.09.058); [Novikov et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.03.006); [Veraart et al., 2018](https://doi.org/10.1016/j.neuroimage.2017.09.030), to mention a few). From the SM point of view, earlier models impose constraints either on compartment parameters or the functional form of the fiber ODF; such constraints improve robustness but may introduce biases into the estimation of remaining parameters.

For more details please look at our recent publication: [Reproducibility of the Standard Model of diffusion in white matter on clinical MRI systems, (2022), ArXiv](https://arxiv.org/pdf/2202.02399.pdf), or Section 3 in this review by [Novikov et al. (2019)](https://doi.org/10.1002/nbm.3998).

#### Assumptions in a nutshell
- Sufficient coarse-graining for Gaussian diffusion at the fiber bundle level (thus, no time-dependence)
- All fiber bundles inside a voxel have the same microstructural properties
- Negligible water exchange between compartments in a bundle
- No assumptions on the fibre bundle's diffusivities
- No assumptions on the functional form of the ODF

_Note that this model does not apply to gray matter, where exchange and the presence of soma cannot be ignored._

<br>


## SMI input data


### Minimal requrements for linear tensor encoding (LTE)
- This implementation of the SM supports as input a 4D array [nx x ny x nz x N] of diffusion-weighted data, meaning a 3D spatial arrangement of voxels with N diffusion measurements arranged along 4th dimention.
- Protocol files (FSL format) are needed: a b-value (.bval) file, a [1 x N] vector, and a directions file (.bvec), a [3 x N] array.
- A noise map, this must be a 3D array of size [nx x ny x nz]. This can be previously estimated using MPPCA, see [this Matlab implementation](https://github.com/NYU-DiffusionMRI/mppca_denoise). 

_On the b-value units:_ We recommend microstructural units [milliseconds / (squared micrometers)].
Note that typical scanner units are b=1000 [seconds / (squared millimeters)] which equal 1 [milliseconds / (squared micrometers)], thus to convert from scanner units to microstructural units one should divide by 1000. 

### Data with non-LTE encodings and variable echo times (TE)
Input data can also have:
- Multiple **B**-tensor shapes (specified by β). This input must be a [1 x N] vector. Only axially symmetric b-tensors are supported. β is a unitless scalar between -0.5 and 1 that indicates the **B**-tensor shape (see figure below).
- Multiple echo times (TE, assumed to be in [milliseconds]). This input must be a [1 x N] vector.

In this general scenario, each measurement is thus fully specified by: a b-value (b), a unit direction (**u**) (axis of symmetry of **B**), a b-tensor shape (β), and TE. See the figure and equation below to understand how these parameters make a b-tensor **B**:
<p align="center">
  <img width="500" alt=" AxSymB_wEq" src="https://user-images.githubusercontent.com/54751227/152437987-d79193d1-1ecc-4707-bdc3-f7cd2dec6ad6.png">
</p>

  - If no β is supplied (this input can be absent or an empty array) the code assumes β=1 (linear tensor encoding, LTE, measurements).
  - If no TE is supplied (this input can be absent or an empty array), the code assumes that the TE is the same across all measurements. In this case compartmental T2 values will not be outputted and water fractions will be T2-weighted.

<br>


## SMI outputs
Standard Model parameters (diffusion) + compartmental T2 values (only if multiple TE data was provided). See figure below for some examples:
<img width="1690" alt="SM_maps_github" src="https://user-images.githubusercontent.com/54751227/152456997-5f24f886-03f9-4eb2-a5f8-1dd767134eae.png">
- Fractions (f, fw) and anisotropy (p2) are adimensional.
- Diffusivities are returned in microstructural units [squared micrometers / milliseconds] (independently of input units).
- Compartmental T2 values are returned in [milliseconds].
- Note that compartmental T2 maps will only be outputed if variable TE data was used.

<br>

## Recommended usage
We recommend using the default options but the code provides users with some flexibility.
Recommended inputs:
- Diffusion data (4D array) + protocol information + mask (binary 3D array) + noise map (3D array).
  - If a mask is not provided the fit will be performed in all voxels present in the first 3 dimensions of the 4D array.
  - If a noise map is not provided (not recommended), it will be estimated using the repetitions of the non-diffusion-weighted images.
  - We recommend preprocessing your raw data with [DESIGNER](https://github.com/NYU-DiffusionMRI/DESIGNER) before doing SM estimation. DESIGNER outputs a robust noise estimation, otherwise [MPPCA](https://github.com/NYU-DiffusionMRI/mppca_denoise) can be used to robustly estimate the noise level.

- The current SMI implementation is written in Matlab, future work may translate it to other languages.

<br>

## Example usage[^note]
We provide three example datasets [here](https://cai2r.net/resources/standard-model-of-diffusion-in-white-matter-the-smi-toolbox/) that were preprocessed with [DESIGNER](https://github.com/NYU-DiffusionMRI/DESIGNER). These contain
- Multiple b-values
- Multiple b-values and tensor shapes
- Multiple b-values, tensor shapes, and echo times.

See the example file 'example.m' with the SMI fit. Note that data is provided in nifti format, we used [this nifti toolbox](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) but feel free to use any and modify the lines where the data is loaded.

```
% =====================================================================
% Minimal usage:

% Add SMI.m to the path, e.g.:
addpath('/Documents/SantiagoCoelho/Git/SMI')

% Load dwi, protocol, and mask

% Specify protocol information
options.b    = bval;
options.dirs = dirs;

% Specify mask and noise map
options.mask  = logical(mask);
options.sigma = sigma;

% Run SM fitting (dwi is a 4D array)
[out] = SMI.fit(dwi,options);
% =====================================================================
% More advanced usage:

% Add SMI.m to the path, e.g.:
addpath('/Documents/SantiagoCoelho/Git/SMI')

% Load dwi, protocol, and mask

% Specify protocol information
options.b    = bval;
options.beta = beta;
options.dirs = dirs;
options.TE   = TE;

% Specify mask and noise map
options.mask  = logical(mask);
options.sigma = sigma;

% Specify options for the fit
options.compartments = {'IAS','EAS','FW'}; % The order does not matter
options.NoiseBias    = 'Rician'; % use 'None' for zero-mean noise

% Run SM fitting (dwi is a 4D array)
[out] = SMI.fit(dwi,options);
% =====================================================================
```
<br>

## Parameter estimation
First, the rotational invariants of the diffusion signal are estimated. This is done using a least squares estimator. If data has a non-negligible rician bias, we suggest to correct it during the fitting (see advanced options below).

Then, the SM parameters are estimated from the rotational invariants of the signal. Unlike conventional parameter estimation approaches which rely on an analytical forward model, e.g. maximum likelihood, here we use data-driven machine learning (ML) regression. This is done by applying a sufficiently flexible nonlinear regression to _training data_ generated with the forward model of interest, considering a wide distribution of model parameters, the noise level, and the protocol that was used. Then, such regression is applied to the data of interest.

For typical SNR values found in clinical dMRI experiments, the optimal regression, i.e. minimizing MSE of the training data, can be achieved already by a cubic polynomial, as it  captures all relevant degrees of freedom in the data represented by the set of rotational invariants.



## Some advanced usage options
The code provides some additional flexibility:
- Rician bias correction: to de-bias the DWI before the spherical harmonics fit.
- Variable number of compartments: 'IAS', 'EAS', 'FW'. _This will be extended but at the moment the only two options are {'IAS', 'EAS'} or {'IAS', 'EAS', 'FW'}._
- User defined parameter distributions for the training data (for the machine learning estimator that performs RotInvs -> kernel).
- **(NOT READY YET)** Batch processing. Parameter estimation for multiple datasets with identical protocols. Here the machine learning training is done only once, regression coefficients are stored and applied to all.
- **(NOT READY YET)** Output spherical harmonic decomposition of the ODF for fiber tracking (normalized for using it with [MRtrix3](https://mrtrix.readthedocs.io/en/0.3.16/workflows/global_tractography.html)).

<br>

## SMI Authors
- [Santiago Coelho](https://santiagocoelho.github.io/)
- [Jelle Veraart](https://med.nyu.edu/faculty/jelle-veraart)
- [Els Fieremans](https://www.diffusion-mri.com/who-we-are/els-fieremans/)
- [Dmitry Novikov](https://www.diffusion-mri.com/who-we-are/dmitry-novikov/)

Do not hesitate to reach out to Santiago.Coelho@nyulangone.org (or [@santicoelho](https://twitter.com/santicoelho) in Twitter) for feedback, suggestions, questions, or comments[^note].

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
% - Coelho, S., Baete, S., Lemberskiy, G., Ades-Aron, B., Barrol, G., Veraart, J., Novikov, D.S., Fieremans, E., 2022. Reproducibility of the Standard Model of diffusion in white matter on clinical MRI systems, 2022, ArXiv.
% - Novikov, D.S., Veraart, J., Jelescu, I.O., Fieremans, E., 2018. Rotationally-invariant mapping of scalar and orientational metrics of neuronal microstructure with diffusion MRI. NeuroImage 174, 518 – 538.
% - Reisert, M., Kellner, E., Dhital, B., Hennig, J., Kiselev, V.G., 2017. Disentangling micro from mesostructure by diffusion MRI: A Bayesian approach. NeuroImage 147, 964 – 975.
```


[^note]:
    Please cite these works if you use the SMI toolbox in your publication:
    - Coelho, S., Baete, S., Lemberskiy, G., Ades-Aron, B., Barrol, G., Veraart, J., Novikov, D.S., Fieremans, E., 2022. Reproducibility of the Standard Model of diffusion in white matter on clinical MRI systems, 2022, ArXiv.
    - Dmitry S. Novikov, Jelle Veraart, Ileana O. Jelescu, Els Fieremans, Rotationally-invariant mapping of scalar and orientational metrics of neuronal microstructure with diffusion MRI, NeuroImage, Volume 174, 2018, Pages 518-538.
    - Marco Reisert, Elias Kellner, Bibek Dhital, Jürgen Hennig, Valerij G. Kiselev, Disentangling micro from mesostructure by diffusion MRI: A Bayesian approach, NeuroImage, Volume 147, 2017, Pages 964-975.
    
[^2]: The SM name originates from a tongue-in-cheek association with the Standard Model in particle physics, as both encompass a fair bit of previous modeling effort from various groups. Here “standard” refers to common assumptions among modeling approaches and does not imply “exact”.
