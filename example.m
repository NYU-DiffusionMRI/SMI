%% Example code for Standard Model (SM) parameter estimation with SMI
%
% EXAMPLE 1 - Fitting SM on fixed TE data (DKI like dataset with only two
%             shells)
%
% EXAMPLE 2 - Fitting SM on fixed TE data (multiple b-values and b-tensors)
%
%
% EXAMPLE 3 - Fitting SM on variable TE data (multiple b-values, b-tensors,
%             and TEs)
%
%
%% EXAMPLE 1
clc,clear,close all

% Path were the test dataset is located
pathFiles='/Users/coelhs01/Documents/SantiagoCoelho/Git/SMI_3datasets/dataset_1';

% Load data and protocol
nii = load_untouch_nii(fullfile(pathFiles,'dwi.nii'));
dwi = abs(double(nii.img));
bval = load(fullfile(pathFiles,'dwi.bval'));
dirs = load(fullfile(pathFiles,'dwi.bvec'));
% Load mask
nii_mask = load_untouch_nii(fullfile(pathFiles,'mask.nii'));
% Load noise map
sigma_nii = load_untouch_nii(fullfile(pathFiles,'sigma_dki.nii'));

% Specify protocol information (no need to specify beta and TE if all
% measurements are LTE and have the same TE)
options.b    = bval;
options.beta = [];
options.dirs = dirs;
options.TE   = [];
options.MergeDistance = 0.1; % If []: default is 0.05 [ms/um^2], this 
% is the threshold for considering different b-values as the same shell
 
% Specify mask and noise map
options.mask  = logical(nii_mask.img);
options.sigma = abs(sigma_nii.img);


% Specify options for the fit
options.compartments = {'IAS','EAS'}; % The order does not matter
options.NoiseBias    = 'Rician'; % this data has Rician noise bias
options.MLTraining.bounds = [0.05, 1, 1, 0.1, 0, 50, 50, 0.05; 0.95, 3, 3, 1.2, 0.5, 150, 120, 0.99];
% The order is: [f, Da, Depar, Deperp, fw, T2a, T2e, p2] (If data has
% fixed TE then the T2a and T2e priors are simply ignored)

% Run SM fitting
tic
[out] = SMI.fit(dwi,options);
t=toc;
fprintf('Time SM fit %f s\n',t)

% Load fa to make an approximate WM mask for plotting results
nii_fa = load_untouch_nii(fullfile(pathFiles,'fa.nii'));
WM_mask=single(nii_fa.img>0.2);
WM_mask(WM_mask(:)==0)=NaN;
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$p_2$','$p_4$'};
clims=[0 1;0 3;0 3;0 1.5;0 1;0 1]; slice=40; Nrows=2;
% Plot results
figure('Position',[506 225 1347 905]), SMI.plotSlices(out.kernel(:,:,:,[1 2 3 4 6 7]).*WM_mask, slice,clims,paramNames,Nrows,[],1,1)

%% EXAMPLE 2
clc,clear,close all

% Path were the test dataset is located
pathFiles='/Users/coelhs01/Documents/SantiagoCoelho/Git/SMI_3datasets/dataset_2';

% Load data and protocol
nii = load_untouch_nii(fullfile(pathFiles,'dwi.nii'));
dwi = abs(double(nii.img));
bval = load(fullfile(pathFiles,'dwi.bval'));
beta = load(fullfile(pathFiles,'dwi.beta'));
dirs = load(fullfile(pathFiles,'dwi.bvec'));
% Load mask
nii_mask = load_untouch_nii(fullfile(pathFiles,'mask.nii'));
% Load noise map
sigma_nii = load_untouch_nii(fullfile(pathFiles,'sigma_dki.nii'));

% Specify protocol information (no need to specify beta and TE if all
% measurements are LTE and have the same TE)
options.b    = bval;
options.beta = beta;
options.dirs = dirs;
options.TE   = [];
options.MergeDistance = 0.1; % If []: default is 0.05 [ms/um^2], this 
% is the threshold for considering different b-values as the same shell

% Specify mask and noise map
options.mask  = logical(nii_mask.img);
options.sigma = abs(sigma_nii.img);


% Specify options for the fit
options.compartments = {'IAS','EAS','FW'}; % The order does not matter
options.D_FW = 3; % Free water diffusivity at body temperature
options.NoiseBias    = 'None'; % the example data has ~ zero-mean noise
options.MLTraining.bounds = [0.05, 1, 1, 0.1, 0, 50, 50, 0.05; 0.95, 3, 3, 1.2, 0.5, 150, 120, 0.99];
% The order is: [f, Da, Depar, Deperp, fw, T2a, T2e, p2] (If data has
% fixed TE then the T2a and T2e priors are simply ignored)

% Run SM fitting
tic
[out] = SMI.fit(dwi,options);
t=toc;
fprintf('Time SM fit %f s\n',t)

% Load fa to make an approximate WM mask for plotting results
nii_fa = load_untouch_nii(fullfile(pathFiles,'fa.nii'));
WM_mask=single(nii_fa.img>0.2);
WM_mask(WM_mask(:)==0)=NaN;
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$','$p_2$','$p_4$'};
clims=[0 1;0 3;0 3;0 1.5;0 1;0 1;0 1]; slice=40; Nrows=2;
% Plot results
figure('Position',[506 225 1347 905]), SMI.plotSlices(out.kernel.*WM_mask, slice,clims,paramNames,Nrows,[],1,1)

%% EXAMPLE 3
clc,clear,close all

% Path were the test dataset is located
pathFiles='/Users/coelhs01/Documents/SantiagoCoelho/Git/SMI_3datasets/dataset_3';

% Load data and protocol
nii = load_untouch_nii(fullfile(pathFiles,'dwi.nii'));
dwi = abs(double(nii.img));
bval = load(fullfile(pathFiles,'dwi.bval'));
TE = load(fullfile(pathFiles,'dwi.TE'));
beta = load(fullfile(pathFiles,'dwi.beta'));
dirs = load(fullfile(pathFiles,'dwi.bvec'));
% Load mask
nii_mask = load_untouch_nii(fullfile(pathFiles,'mask.nii'));
% Load noise map
sigma_nii = load_untouch_nii(fullfile(pathFiles,'sigma_dki.nii'));

% Specify protocol information
options.b    = bval;
options.beta = beta;
options.dirs = dirs;
options.TE   = TE;
options.MergeDistance = 0.1; % If []: default is 0.05 [ms/um^2], this 

% Specify mask and noise map
options.mask  = logical(nii_mask.img);
options.sigma = abs(sigma_nii.img);

% Specify options for the fit
options.compartments = {'IAS','EAS','FW'}; % The order does not matter
options.D_FW = 3; % Free water diffusivity at body temperature
options.NoiseBias    = 'None'; % the example data has ~ zero-mean noise
options.MLTraining.bounds = [0.05   1      1      0.1      0       50    50    0.05;0.95   3      3      1.2      1    150   120    0.99];
% The order is: [f, Da, Depar, Deperp, fw, T2a, T2e, p2] (If data has
% fixed TE then the T2a and T2e priors are simply ignored)
options.RotInv_Lmax=2;

% Run SM fitting
tic
[out] = SMI.fit(dwi,options);
t=toc;
fprintf('Time SM fit %f s\n',t)

% Load fa to make an approximate WM mask for plotting results
nii_fa = load_untouch_nii(fullfile(pathFiles,'fa.nii'));
WM_mask=single(nii_fa.img>0.2);
WM_mask(WM_mask(:)==0)=NaN;
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$',...
            '$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$','$T_2^\mathrm{a}$ [ms]','$T_2^\mathrm{e}$ [ms]','$p_2$'};
clims=[0 1;0 3;0 3;0 1.5;0 1;0 150;0 150;0 1]; slice=40; Nrows=2;
% Plot results
figure('Position',[269 281 1785 821]), SMI.plotSlices(out.kernel.*WM_mask, slice,clims,paramNames,Nrows,[],1,1)

%% EXAMPLE 4 - (same data as EXAMPLE 2)
% This adds a maximum likelihood fit on Slm(b), initialized with machine learning fits of kernel and fODF
clc,clear,close all

% Path were the test dataset is located
pathFiles='/Users/coelhs01/Documents/SantiagoCoelho/Git/SMI_3datasets/dataset_2';

% Load data and protocol
nii = load_untouch_nii(fullfile(pathFiles,'dwi.nii'));
dwi = abs(double(nii.img));
bval = load(fullfile(pathFiles,'dwi.bval'));
beta = load(fullfile(pathFiles,'dwi.beta'));
dirs = load(fullfile(pathFiles,'dwi.bvec'));
% Load mask
nii_mask = load_untouch_nii(fullfile(pathFiles,'mask.nii'));
% Load noise map
sigma_nii = load_untouch_nii(fullfile(pathFiles,'sigma_dki.nii'));

% Specify protocol information (no need to specify beta and TE if all
% measurements are LTE and have the same TE)
options.b    = bval;
options.beta = beta;
options.dirs = dirs;
options.TE   = [];
options.MergeDistance = 0.1; % If []: default is 0.05 [ms/um^2], this 
% is the threshold for considering different b-values as the same shell

% Specify WM mask to run fitting and noise map
nii_fa = load_untouch_nii(fullfile(pathFiles,'fa.nii'));
WM_mask=single(nii_fa.img>0.2);
options.mask  = logical(WM_mask);
options.sigma = abs(sigma_nii.img);

% Keep only 1 slice to make fit run faster
Slice_keep = 42; Nslices = size(WM_mask,3);
options.mask(:,:,[1:(Slice_keep-1) (Slice_keep+1):Nslices]) = 0;

% Specify options for the fit
options.compartments = {'IAS','EAS','FW'}; % The order does not matter
options.D_FW = 3; % Free water diffusivity at body temperature
options.NoiseBias    = 'None'; % the example data has ~ zero-mean noise
options.MLTraining.bounds = [0.05, 1, 1, 0.1, 0, 50, 50, 0.05; 0.95, 3, 3, 1.2, 0.5, 150, 120, 0.99];
% The order is: [f, Da, Depar, Deperp, fw, T2a, T2e, p2] (If data has fixed TE then the T2a and T2e priors are simply ignored)

% Fit fODF
options.flag_fit_fODF = 1;

% options.Nlevels = 1;
options.RotInv_Lmax = 4;
options.MLTraining.Ntraining = 5e4;

% Run Maximum Likelihood fit starting from the polynomial regression fit (takes significantly more time)
options.run_maximumlikelihood = 1;

% Run SM fitting
tic
[out] = SMI.fit(dwi,options);
t=toc;
fprintf('Time SM fit %f s\n',t)

% Plot results from Polynomial Regression
WM_mask(WM_mask(:)==0)=NaN;
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$','$p_2$','$p_4$'};
clims=[0 1;0 3;0 3;0 1.5;0 1;0 1;0 1]; Slice_keep=42; Nrows=2;
figure('Position',[506 225 1347 905]), SMI.plotSlices(out.kernel(:,:,:,1:6).*WM_mask, Slice_keep,clims,paramNames,Nrows,[],1,1)
sgt= sgtitle('Machine Learning fit'); sgt.FontSize = 40; sgt.Interpreter='latex';

% Plot results from Maximum Likelihood
WM_mask(WM_mask(:)==0)=NaN;
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{e}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$','$p_2$','$p_4$'};
clims=[0 1;0 3;0 3;0 1.5;0 1;0 1;0 1]; Slice_keep=42; Nrows=2;
figure('Position',[506 225 1347 905]), SMI.plotSlices(cat(4,out.KERNEL_maxlik(:,:,:,2:6),out.p2_maxlik).*WM_mask, Slice_keep,clims,paramNames,Nrows,[],1,1)
sgt= sgtitle('Maximum Likelihood fit'); sgt.FontSize = 40; sgt.Interpreter='latex';















