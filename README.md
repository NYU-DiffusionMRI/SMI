# Standard Model Imaging (SMI) toolbox
## This MATLAB toolbox contains all necessary functions for parameter estimation of the Standard Model of diffusion in white matter

### The Standard Model of diffusion in white matter
Multiple approaches to model the physics of water diffusion in white matter rely on similar assumptions. This led to the unifying framework dubbed Standard Model (SM) of diffusion in WM as formulated in (Novikov et al., 2019). Briefly, axons (and possibly glial processes) are represented by impermeable zero-radius cylinders (the so-called “sticks”) arranged in locally coherent fiber fascicles. The diffusion in the extra-axonal space of each fascicle is assumed to be Gaussian and described by an axially symmetric diffusion tensor. The third, optional tissue compartment is the cerebro-spinal fluid (CSF). Such multicomponent fascicles are distributed in a voxel according to an arbitrary fiber orientation distribution function (ODF). All fascicles in a voxel are assumed to have the same compartment fractions and diffusivities, and differ from each other only by orientation.

The SM encompasses a number of WM models made of anisotropic Gaussian compartments with axons represented by sticks (Kroenke et al., 2004; Jespersen et al., 2007, 2010; Fieremans et al., 2011; Zhang et al., 2012; Sotiropoulos et al., 2012; Jensen et al., 2016; Jelescu et al., 2016a; Kaden et al., 2016; Reisert et al., 2017; Novikov et al., 2018; Ve- raart et al., 2018). From the SM point of view, earlier models impose constraints either on compartment parameters or the functional form of the fiber ODF; such constraints improve robustness but may introduce biases into the estimation of remaining parameters.


### SMI input data
This implementation of the SM supports diffusion data with varying b-values, b-tensor shapes (\beta), anc echo times (TE).



We recommend using the default options but the code is sufficiently flexible to provide users with some flexibility.

The current SMI implementation is written in Matlab, future work may translate it to other languages.

For technical details please look at the following publications:

- Arxiv SM reproducibility
- RotInv NeuroImage 2018

Do not hesitate to reach out to Santiago.Coelho@nyulangone.org for feedback, suggestions, questions, or comments.
