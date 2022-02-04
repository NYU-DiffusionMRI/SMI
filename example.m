%% Example code for Standard Model (SM) parameter estimation
%
% EXAMPLE 1 - Fitting SM on fixed TE data
% EXAMPLE 2 - Fitting SM on variable TE data
%
%
% Load test dataset from 'https://drive.google.com/drive/folders/1TQzZGM7PdTf1kplwfWLIIRn8q0kE3Nix?usp=sharing'
%


%% EXAMPLE 1 - Fitting SM on fixed TE data
clc,clear,close all
pathFiles='/Users/coelhs01/Documents/SantiagoCoelho/Git/exampleDataSMI';
% =========================================================================
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
% =========================================================================
% Keep only fixed TE data (TE=92ms)
keep_fixed_TE=TE==92;
dwi=dwi(:,:,:,keep_fixed_TE);
bval=bval(keep_fixed_TE);
dirs=dirs(:,keep_fixed_TE);
bshape=bshape(keep_fixed_TE);
TE=[];
% =========================================================================

% Perform Spherical Harmonics fit
Lmax = SMI.GetDefaultLmax(bval,bshape,TE);
tic
[Slm,Sl,~,table_4D_sorted] = SMI.Fit2D4D_LLS_RealSphHarm_wSorting_norm_varL(dwi,mask,bval,dirs,bshape,TE,Lmax);
t=toc;
fprintf('Time SH fit %f s\n',t)
S0nn=squeeze(Sl(:,:,:,1,:));
S2nn=squeeze(Sl(:,:,:,2,:));
S4nn=squeeze(Sl(:,:,:,3,:));
S0_lowest=S0nn(:,:,:,1);
IMGUI(S0nn./S0_lowest,[0 0.5]), set(gcf,'Position',[40.5714 26.3333 294.1429 59.8000])
IMGUI(S2nn./S0_lowest,[0 0.1]), set(gcf,'Position',[40.5714 26.3333 294.1429 59.8000])
IMGUI(S4nn./S0_lowest,[0 0.02]), set(gcf,'Position',[40.5714 26.3333 294.1429 59.8000])
% IMGUI(cat(4,20*S2nn(:,:,:,6)./S0_lowest,50*S4nn(:,:,:,6)./S0_lowest),[0 1]), set(gcf,'Position',[40.5714 26.3333 294.1429 59.8000])
% =========================================================================
% Perform PR training on Rotational Invariants and fitting SM kernel
S0nn=squeeze(Sl(:,:,:,1,:));
S2nn=squeeze(Sl(:,:,:,2,:));
RotInvs=cat(4,S0nn,S2nn);
lb_training=[0.05   1      1      0.1      0       50    50    0.05];
ub_training=[0.95   3      3      1.2      0.3    150   120    0.99];
Lmax_train=6;
Ntraining=2.5e4;
Nlevels=10;
tic
KERNEL = SMI.StandardModel_PR_fit_RotInvs(RotInvs,mask,sigma,bval,dirs,bshape,TE,lb_training,ub_training,Lmax_train,Ntraining,Nlevels,[0 0.4]);
t=toc;
fprintf('Time PR training fit %f s\n',t)
KERNEL_resc=KERNEL.*permute([1 1/3 1/3 1 1 1],[1 3 4 2]);
IMGUI(KERNEL_resc,[0 1]), set(gcf, 'Position', [30.5714 7.4667 258.2857 71.8000])


%% EXAMPLE 2 - Fitting SM on variable TE data
clc,clear,close all
pathFiles='/Users/coelhs01/Documents/SantiagoCoelho/Git/exampleDataSMI';
% =========================================================================
% Load data and protocol
nii=load_untouch_nii(fullfile(pathFiles,'dwi_preproc.nii.gz'));
bval=load(fullfile(pathFiles,'dwi_preproc.bval'));
TE=load(fullfile(pathFiles,'dwi_preproc.TE'));
bshape=load(fullfile(pathFiles,'dwi_preproc.bshape'));
dirs=load(fullfile(pathFiles,'dwi_preproc.bvec'));
dwi=abs(double(nii.img));
% save_3D4D_nii_with_template_nii_header(nii,fullfile(pathFiles,'dwi_preproc.nii.gz'),dwi(:,:,40:43,:))

Ndwi=size(dwi,4);
nii.img=0;
% Load mask
nii_mask=load_untouch_nii(fullfile(pathFiles,'mask.nii.gz'));
mask=logical(nii_mask.img);
se=strel('cube',2);
emask=imerode(mask,se);
% Load noise map and p-map
sigma_nii=load_untouch_nii(fullfile(pathFiles,'sigma_dki.nii'));
sigma=abs(sigma_nii.img);
sigma_mf = medfilt3(sigma,[3 3 3]);
nii_pmap=load_untouch_nii(fullfile(pathFiles,'Npars.nii'));
p=nii_pmap.img;
ps=7^3;
sigma_factor=sqrt(p*(1/Ndwi+1/ps));
% =========================================================================
% Perform Spherical Harmonics fit
slice=33;
Lmax=[zeros(1,4) 2 2 2 0 4 4 4 4 4 6];
tic
[Slm,Sl,~,table_4D_sorted] = SMI.Fit2D4D_LLS_RealSphHarm_wSorting_norm_varL(dwi,mask,bval,dirs,bshape,TE,Lmax);
t=toc;
fprintf('Time SH fit %f s\n',t)
S0nn=squeeze(Sl(:,:,:,1,:));
S2nn=squeeze(Sl(:,:,:,2,:));
S0_lowest=S0nn(:,:,:,1);
IMGUI(S0nn./S0_lowest,[0 0.5]), set(gcf,'Position',[40.5714 26.3333 294.1429 59.8000])
IMGUI(S2nn./S0_lowest,[0 0.1]), set(gcf,'Position',[40.5714 26.3333 294.1429 59.8000])
% =========================================================================
% Perform PR training on Rotational Invariants and fitting SM kernel
S0nn=squeeze(Sl(:,:,:,1,:));
S2nn=squeeze(Sl(:,:,:,2,:));
RotInvs=cat(4,S0nn,S2nn);
lb_training=[0.05   1      1      0.1      0       50    50    0.05];
ub_training=[0.95   3      3      1.2      1      150   120    0.99];
Lmax_train=6;
Ntraining=2.5e4;
Nlevels=10;
tic
KERNEL = SMI.StandardModel_PR_fit_RotInvs(RotInvs,mask,sigma_mf,bval,dirs,bshape,TE,lb_training,ub_training,Lmax_train,Ntraining,Nlevels,[0 0.4]);
t=toc;
fprintf('Time PR training fit %f s\n',t)
KERNEL_resc=KERNEL.*permute([1 1/3 1/3 1 1 1/150 1/150 1],[1 3 4 2]);
IMGUI(KERNEL_resc,[0 1]), set(gcf, 'Position', [30.5714 7.4667 258.2857 71.8000])















