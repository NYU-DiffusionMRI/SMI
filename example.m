%% Example code for Standard Model (SM) parameter estimation
%
% EXAMPLE 1 - Fitting SM on fixed TE data
% EXAMPLE 2 - Fitting SM on variable TE data
%
% Download test dataset from 'https://drive.google.com/drive/folders/1TQzZGM7PdTf1kplwfWLIIRn8q0kE3Nix?usp=sharing'
%

%% EXAMPLE 1 - Fitting SM on fixed TE data
clc,clear,close all

% Path were the test dataset is located
pathFiles='/Users/coelhs01/Documents/SantiagoCoelho/Git/exampleDataSMI';

% Load data and protocol
nii = load_untouch_nii(fullfile(pathFiles,'dwi_preproc.nii.gz'));
dwi = abs(double(nii.img));
bval = load(fullfile(pathFiles,'dwi_preproc.bval'));
TE = load(fullfile(pathFiles,'dwi_preproc.TE'));
beta = load(fullfile(pathFiles,'dwi_preproc.bshape'));
dirs = load(fullfile(pathFiles,'dwi_preproc.bvec'));
% Load mask
nii_mask = load_untouch_nii(fullfile(pathFiles,'mask.nii.gz'));
% Load noise map
sigma_nii = load_untouch_nii(fullfile(pathFiles,'sigma_dki.nii'));

% Keep only the subset of the data that has TE=92ms
keep_fixed_TE = TE==92;
dwi = dwi(:,:,:,keep_fixed_TE);
bval = bval(keep_fixed_TE);
dirs = dirs(:,keep_fixed_TE);
beta = beta(keep_fixed_TE);
TE = []; % We do not need TE since all data has equal TE

% Specify protocol information
options.b    = bval;
options.beta = beta;
options.dirs = dirs;
options.TE   = TE;

% Specify mask and noise map
options.mask  = logical(nii_mask.img);
options.sigma = abs(sigma_nii.img);


% Specify options for the fit
options.compartments = {'IAS','EAS','FW'}; % The order does not matter
options.NoiseBias    = 'None'; % the example data has ~ zero-mean noise
options.MLTraining.bounds = [0.05   1      1      0.1      0       50    50    0.05;0.95   3      3      1.2      1    150   120    0.99];

% Run SM fitting
tic
[out] = SMI.fit(dwi,options);
t=toc;
fprintf('Time SM fit %f s\n',t)

% Load fa to make an approximate WM mask for plotting results
nii_fa = load_untouch_nii(fullfile(pathFiles,'fa.nii.gz'));
WM_mask=single(nii_fa.img>0.2);
WM_mask(WM_mask(:)==0)=NaN;
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{2}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{2}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$','$p_2$'};
clims=[0 1;0 3;0 3;0 1.5;0 1;0 1]; slice=40; Nrows=2;
% Plot results
figure('Position',[506 225 1347 905]), SMI.plotSlices(out.kernel.*WM_mask, slice,clims,paramNames,Nrows,[],1,1)

%% EXAMPLE 2 - Fitting SM on variable TE data
clc,clear,close all

% Path were the test dataset is located
pathFiles='/Users/coelhs01/Documents/SantiagoCoelho/Git/exampleDataSMI';

% Load data and protocol
nii = load_untouch_nii(fullfile(pathFiles,'dwi_preproc.nii.gz'));
dwi = abs(double(nii.img));
bval = load(fullfile(pathFiles,'dwi_preproc.bval'));
TE = load(fullfile(pathFiles,'dwi_preproc.TE'));
beta = load(fullfile(pathFiles,'dwi_preproc.bshape'));
dirs = load(fullfile(pathFiles,'dwi_preproc.bvec'));
% Load mask
nii_mask = load_untouch_nii(fullfile(pathFiles,'mask.nii.gz'));
% Load noise map
sigma_nii = load_untouch_nii(fullfile(pathFiles,'sigma_dki.nii'));

% Specify protocol information
options.b    = bval;
options.beta = beta;
options.dirs = dirs;
options.TE   = TE;

% Specify mask and noise map
options.mask  = logical(nii_mask.img);
options.sigma = abs(sigma_nii.img);

% Specify options for the fit
options.compartments = {'IAS','EAS','FW'}; % The order does not matter
options.NoiseBias    = 'None'; % the example data has ~ zero-mean noise
options.MLTraining.bounds = [0.05   1      1      0.1      0       50    50    0.05;0.95   3      3      1.2      1    150   120    0.99];

% Run SM fitting
tic
[out] = SMI.fit(dwi,options);
t=toc;
fprintf('Time SM fit %f s\n',t)

% Load fa to make an approximate WM mask for plotting results
nii_fa = load_untouch_nii(fullfile(pathFiles,'fa.nii.gz'));
WM_mask=single(nii_fa.img>0.2);
WM_mask(WM_mask(:)==0)=NaN;
paramNames={'$f$','$D_\mathrm{a}\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{2}^\|\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$D_\mathrm{2}^\perp\,[\mathrm{\mu m}^2/\mathrm{ms}]$','$f_\mathrm{w}$','$T_2^\mathrm{a}$ [ms]','$T_2^\mathrm{e}$ [ms]','$p_2$'};
clims=[0 1;0 3;0 3;0 1.5;0 1;0 150;0 150;0 1]; slice=40; Nrows=2;
% Plot results
figure('Position',[269 281 1785 821]), SMI.plotSlices(out.kernel.*WM_mask, slice,clims,paramNames,Nrows,[],1,1)










