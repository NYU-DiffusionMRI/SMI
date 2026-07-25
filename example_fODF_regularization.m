% =========================================================================
% Example: regularized fODF deconvolution in the SMI toolbox
%
% The fODF is obtained by deconvolving the SM kernel from the DWI, which is
% an ill-conditioned problem: the kernel rotational invariants Kl decay
% quickly with l, so the high order plm are dominated by noise and the
% resulting fODF typically has large negative lobes.
%
% Two regularizers can be added to the deconvolution:
%
%   - Non-negativity, as in constrained spherical deconvolution
%     (Tournier et al., NeuroImage 2007). The fODF amplitude is iteratively
%     penalized on the directions where it becomes negative.
%   - Tikhonov damping of the plm, which suppresses the coefficients that
%     the kernel attenuates the most.
%
% This script builds a synthetic two-fiber voxel, generates its signal with
% the SM forward model, and compares the deconvolution with and without
% regularization. It requires no data and runs in a few seconds.
%
% By: Santiago Coelho
% =========================================================================
clear; close all
% addpath('/Documents/SantiagoCoelho/Git/SMI')

CS_phase = 1; D_FW = 3; Lmax = 6;
L_all = repelem(0:2:Lmax,2*(0:2:Lmax)+1);

% ---- Protocol: 3 shells with 64 directions each + 4 b0s ----------------
Ndirs_shell = 64;
dirs_shell = SMI.GetUniformHemisphereDirs(Ndirs_shell);
bvals = [0 1 2 3]; % [ms/um^2]
b = []; dirs = [];
for ii = 1:length(bvals)
    if bvals(ii)==0
        b = [b zeros(1,4)]; dirs = [dirs; repmat([0 0 1],4,1)];
    else
        b = [b bvals(ii)*ones(1,Ndirs_shell)]; dirs = [dirs; dirs_shell];
    end
end
beta = ones(size(b)); TE = zeros(size(b));
Ndwi = length(b); Lmax_shells = Lmax*ones(1,length(bvals));

% ---- Ground truth: two crossing fibers at 50 degrees --------------------
% The fODF is a pair of smoothed delta functions, so it is non-negative and
% its plm (in the normalized convention used by SMI, where p_00=1) are known
fibers = [1 0 0; cosd(50) sind(50) 0]; weights = [0.6 0.4]; smoothing = 0.05;
plm_SH = (weights*SMI.get_even_SH(fibers,Lmax,CS_phase)).*exp(-smoothing*L_all.*(L_all+1));
plm_all = plm_SH./sqrt((2*L_all+1)/(4*pi));
plm_true = plm_all(2:end);

% ---- Kernel and noise-free signal --------------------------------------
f = 0.6; Da = 2; Depar = 2; Deperp = 0.5; fw = 0; T2a = 1; T2e = 1;
S = SMI.SM_wFW_b_beta_TE_RealSphHarm(f,Da,Depar,Deperp,1-f-fw,T2a,T2e,plm_true,b,dirs,beta,TE,CS_phase,D_FW)';
kernel = reshape([f Da Depar Deperp fw T2a T2e],[1 1 1 7]);
mask = true(1,1,1);

% ---- Directions where the fODF is evaluated for the comparison ---------
dirs_odf = SMI.GetUniformHemisphereDirs(500);
Y_odf = SMI.get_even_SH(dirs_odf,Lmax,CS_phase).*sqrt((2*L_all+1)/(4*pi));
fODF = @(plm) Y_odf(:,1) + Y_odf(:,2:end)*plm(:); % plm are normalized, p_00 contributes 1/(4*pi)
negative_mass = @(plm) -sum(min(fODF(plm),0))/sum(abs(fODF(plm)));
fODF_true = fODF(plm_true);

% ---- Noise-free check: the deconvolution is exact ----------------------
plm_hat = squeeze(SMI.get_plm_from_S_and_kernel(reshape(S,[1 1 1 Ndwi]),Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW))';
fprintf('Noise free, unregularized: max|plm-plm_true| = %.2e\n\n',max(abs(plm_hat-plm_true)));

% ---- Noisy comparison ---------------------------------------------------
SNR = 30; Nrep = 40;
rng(1)
noise = randn(Nrep,Ndwi)/SNR;

% Options of each deconvolution that is compared
reg_tikhonov.lambda_tikhonov = 1;
reg_nonneg.flag_nonneg = 1;
reg_both.flag_nonneg = 1; reg_both.lambda_tikhonov = 0.3; reg_both.TikhonovMatrix = 'laplacebeltrami';
configs = {[],reg_tikhonov,reg_nonneg,reg_both};
names = {'unregularized  ','tikhonov       ','non-negativity ','nonneg+tikhonov'};

fprintf('SNR = %d, %d noise realizations, Lmax = %d\n',SNR,Nrep,Lmax)
fprintf('(negative mass is the fraction of the fODF absolute mass that is negative,\n')
fprintf(' it is %.4f for the ground truth)\n\n',negative_mass(plm_true));
fprintf('%s   RMSE(plm)   rel err fODF   negative mass\n',repmat(' ',1,15))
plm_examples = zeros(length(configs),length(plm_true));
for cc = 1:length(configs)
    err = 0; err_odf = 0; negmass = 0;
    for rr = 1:Nrep
        dwi = reshape(S+noise(rr,:),[1 1 1 Ndwi]);
        if isempty(configs{cc})
            p = SMI.get_plm_from_S_and_kernel(dwi,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW);
        else
            p = SMI.get_plm_from_S_and_kernel(dwi,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW,configs{cc});
        end
        p = squeeze(p)';
        err = err + norm(p-plm_true)^2;
        err_odf = err_odf + norm(fODF(p)-fODF_true)^2/norm(fODF_true)^2;
        negmass = negmass + negative_mass(p);
        if rr==1, plm_examples(cc,:) = p; end
    end
    fprintf('%s   %9.4f   %11.4f   %13.4f\n',names{cc},sqrt(err/Nrep),sqrt(err_odf/Nrep),negmass/Nrep);
end

% ---- Plot one realization along the plane containing both fibers -------
phi = linspace(0,2*pi,360)';
dirs_plane = [cos(phi) sin(phi) zeros(size(phi))];
Y_plane = SMI.get_even_SH(dirs_plane,Lmax,CS_phase).*sqrt((2*L_all+1)/(4*pi));
amp = @(plm) Y_plane(:,1) + Y_plane(:,2:end)*plm(:);
figure('Color','w'), hold on
plot(phi*180/pi,amp(plm_true),'k','LineWidth',2)
for cc = 1:length(configs)
    plot(phi*180/pi,amp(plm_examples(cc,:)),'LineWidth',1.2)
end
plot(phi*180/pi,zeros(size(phi)),'k--')
xlabel('angle in the fiber plane [degrees]'), ylabel('fODF')
legend(['ground truth',names],'Location','best'), box on
title('fODF along the plane containing both fibers (single noise realization)')

% =========================================================================
% Usage inside SMI.fit:
%
% options.flag_fit_fODF = 1;
% options.fODF_regularization.flag_nonneg = 1;
% options.fODF_regularization.lambda_tikhonov = 0.3;
% [out] = SMI.fit(dwi,options);
%
% out.fODF_regularization keeps the options that were used together with the
% number of iterations, the number of constrained directions, and a
% convergence flag for each voxel.
% =========================================================================
