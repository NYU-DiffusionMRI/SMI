% =========================================================================
% Example: choosing the fODF deconvolution regularization parameters
%
% example_fODF_regularization.m compares a few fixed settings. This script
% sweeps them, so that the values in the documentation are measured rather
% than guessed, and so that they can be re-measured for a different
% protocol by editing the protocol block below.
%
% The fODF is estimated from synthetic two-fiber voxels whose ground truth
% plm are known exactly, over a range of crossing angles and SNRs, and the
% estimate is scored against that ground truth. The sweep is staged:
%
%   1. joint grid over lambda_nonneg x lambda_tikhonov
%   2. tau, at the best point of that grid
%   3. Lmax_init, at the best point of that grid
%   4. identity vs Laplace-Beltrami Tikhonov matrix
%
% Scores (all averaged over crossing angles, SNRs and noise realizations):
%
%   rel err fODF   relative L2 error of the fODF over the sphere. This is
%                  the primary score: it is what "reconstructs the fODF
%                  similar to the ground truth" means, and it is
%                  dimensionless.
%   RMSE(plm)      error of the coefficients themselves
%   negative mass  fraction of the absolute mass of the fODF that is
%                  negative, i.e. how unphysical the estimate is
%   peak error     mean angle between each true fiber direction and the
%                  closest peak of the estimated fODF. This is the
%                  quantity tractography actually consumes.
%
% Scoring lives in fODF_regularization_score.m (a separate file rather
% than a local function, so that the script also runs in GNU Octave).
%
% All the deconvolutions for one parameter setting are batched into a
% single call to SMI.get_plm_from_S_and_kernel by stacking the noise
% realizations along the first image dimension, which is much faster than
% calling it once per realization and gives bit-identical results.
%
% Requires no data. Runtime is a few minutes.
%
% By: Santiago Coelho
% =========================================================================
clear; close all

CS_phase = 1; D_FW = 3; Lmax = 6;
L_all = repelem(0:2:Lmax,2*(0:2:Lmax)+1);
L_rest = L_all(2:end);

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

% ---- Kernel ------------------------------------------------------------
f = 0.6; Da = 2; Depar = 2; Deperp = 0.5; fw = 0; T2a = 1; T2e = 1;
kernel_1vox = reshape([f Da Depar Deperp fw T2a T2e],[1 1 1 7]);

% ---- Conditions the parameters have to work across ---------------------
% A single crossing angle at a single SNR would give an optimum tuned to
% that one case. Averaging over a range gives settings that are not.
angles    = [40 60 90];   % crossing angle [degrees]
SNRs      = [20 30 50];
weights   = [0.6 0.4];    % volume fractions of the two fibers
smoothing = 0.05;         % makes the ground truth fODF non-negative
Nrep      = 30;           % noise realizations per condition

rng_seed = 1;
if exist('rng','builtin') || exist('rng','file'), rng(rng_seed), else, randn('seed',rng_seed); end

% ---- Build every voxel of the experiment once --------------------------
Ncond = numel(angles)*numel(SNRs);
Nvox  = Ncond*Nrep;
DWI   = zeros(Nvox,Ndwi);
PLM_T = zeros(Nvox,numel(L_rest));
COND  = zeros(Nvox,2);          % [angle index, SNR index]
FIB   = cell(numel(angles),1);
vv = 0;
for ia = 1:numel(angles)
    fibers = [1 0 0; cosd(angles(ia)) sind(angles(ia)) 0];
    FIB{ia} = fibers;
    plm_SH  = (weights*SMI.get_even_SH(fibers,Lmax,CS_phase)).*exp(-smoothing*L_all.*(L_all+1));
    plm_all = plm_SH./sqrt((2*L_all+1)/(4*pi));
    plm_a   = plm_all(2:end);
    S_a = SMI.SM_wFW_b_beta_TE_RealSphHarm(f,Da,Depar,Deperp,1-f-fw,T2a,T2e, ...
              plm_a,b,dirs,beta,TE,CS_phase,D_FW)';
    for is = 1:numel(SNRs)
        for rr = 1:Nrep
            vv = vv+1;
            DWI(vv,:)   = S_a + randn(1,Ndwi)/SNRs(is);
            PLM_T(vv,:) = plm_a;
            COND(vv,:)  = [ia is];
        end
    end
end
dwi4  = reshape(DWI,[Nvox 1 1 Ndwi]);
kernel = repmat(kernel_1vox,[Nvox 1 1 1]);
mask   = true(Nvox,1,1);
fprintf('Sweep over %d conditions x %d noise realizations = %d voxels per setting\n', ...
    Ncond, Nrep, Nvox);
fprintf('  crossing angles: %s deg | SNR: %s | Lmax = %d\n\n', ...
    mat2str(angles), mat2str(SNRs), Lmax);

% ---- Directions used to score the fODF ---------------------------------
dirs_odf = SMI.GetUniformHemisphereDirs(500);
Y_odf    = SMI.get_even_SH(dirs_odf,Lmax,CS_phase).*sqrt((2*L_all+1)/(4*pi));
AMP_T    = Y_odf(:,1) + Y_odf(:,2:end)*PLM_T.';        % [Ndirs x Nvox]
nrm_T    = sqrt(sum(AMP_T.^2,1));

% Neighbour list on the scoring directions, for peak detection
nNeighbors = 8;
G = abs(dirs_odf*dirs_odf.'); G(1:size(G,1)+1:end) = -Inf;
[~,ordG] = sort(G,2,'descend');
NB = ordG(:,1:nNeighbors);

% =========================================================================
%  STAGE 1: joint grid, lambda_nonneg x lambda_tikhonov
% =========================================================================
lam_nn  = [0 0.3 1 3 10 30 100]; % 0 = non-negativity constraint disabled
lam_tik = [0 0.1 0.3 1 3 10];  % 0 = Tikhonov disabled

E_odf = nan(numel(lam_nn),numel(lam_tik));
E_plm = nan(size(E_odf)); E_neg = nan(size(E_odf)); E_pk = nan(size(E_odf));
fprintf('STAGE 1: lambda_nonneg x lambda_tikhonov  (primary score: rel err fODF)\n');
fprintf('%12s','lam_nn\\lam_tik');
fprintf('%9.3g',lam_tik); fprintf('\n');
for in = 1:numel(lam_nn)
    fprintf('%12.3g',lam_nn(in));
    for it = 1:numel(lam_tik)
        reg = struct();
        reg.flag_nonneg     = double(lam_nn(in) > 0);
        if reg.flag_nonneg, reg.lambda_nonneg = lam_nn(in); end
        reg.lambda_tikhonov = lam_tik(it);
        [E_odf(in,it),E_plm(in,it),E_neg(in,it),E_pk(in,it)] = ...
            fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                          PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB);
        fprintf('%9.4f',E_odf(in,it));
    end
    fprintf('\n');
end
[bestE,ix] = min(E_odf(:));
[in_b,it_b] = ind2sub(size(E_odf),ix);
fprintf('\n  best: lambda_nonneg = %g, lambda_tikhonov = %g  -> rel err fODF %.4f\n', ...
    lam_nn(in_b), lam_tik(it_b), bestE);
fprintf('  at that point: RMSE(plm) %.4f | negative mass %.4f | peak error %.2f deg\n\n', ...
    E_plm(in_b,it_b), E_neg(in_b,it_b), E_pk(in_b,it_b));

best.flag_nonneg     = double(lam_nn(in_b) > 0);
best.lambda_nonneg   = max(lam_nn(in_b),1);
best.lambda_tikhonov = lam_tik(it_b);

% =========================================================================
%  STAGE 2: tau
% =========================================================================
taus = [0.02 0.05 0.1 0.2 0.4];
E_tau = nan(size(taus));
fprintf('STAGE 2: tau (at the stage 1 optimum)\n');
for k = 1:numel(taus)
    reg = best; reg.tau = taus(k);
    E_tau(k) = fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                             PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB);
    fprintf('  tau = %5.2f -> rel err fODF %.4f\n', taus(k), E_tau(k));
end
[~,k_b] = min(E_tau); best.tau = taus(k_b);
fprintf('  best tau = %g\n\n', best.tau);

% =========================================================================
%  STAGE 3: Lmax_init
% =========================================================================
Lin = 2:2:Lmax;
E_lin = nan(size(Lin));
fprintf('STAGE 3: Lmax_init (at the stage 1-2 optimum)\n');
for k = 1:numel(Lin)
    reg = best; reg.Lmax_init = Lin(k);
    E_lin(k) = fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                             PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB);
    fprintf('  Lmax_init = %d -> rel err fODF %.4f\n', Lin(k), E_lin(k));
end
[~,k_b] = min(E_lin); best.Lmax_init = Lin(k_b);
fprintf('  best Lmax_init = %d\n\n', best.Lmax_init);

% =========================================================================
%  STAGE 4: Tikhonov matrix
% =========================================================================
mats = {'identity','laplacebeltrami'};
fprintf('STAGE 4: Tikhonov matrix (at the optimum, re-optimizing lambda_tikhonov)\n');
E_mat = nan(numel(mats),numel(lam_tik));
for m = 1:numel(mats)
    for it = 1:numel(lam_tik)
        reg = best; reg.TikhonovMatrix = mats{m}; reg.lambda_tikhonov = lam_tik(it);
        E_mat(m,it) = fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                                    PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB);
    end
    [em,ib] = min(E_mat(m,:));
    fprintf('  %-16s best lambda_tikhonov = %-5g -> rel err fODF %.4f\n', mats{m}, lam_tik(ib), em);
end
[~,m_b] = min(min(E_mat,[],2));
[~,it_b2] = min(E_mat(m_b,:));
best.TikhonovMatrix = mats{m_b}; best.lambda_tikhonov = lam_tik(it_b2);

% =========================================================================
%  Per SNR optimum: the best amount of regularization depends on the noise
% =========================================================================
fprintf('\nPer-SNR optimum over the stage 1 grid (rel err fODF)\n');
for is = 1:numel(SNRs)
    sel = COND(:,2)==is;
    Es = nan(numel(lam_nn),numel(lam_tik));
    for in = 1:numel(lam_nn)
        for it = 1:numel(lam_tik)
            reg = struct();
            reg.flag_nonneg = double(lam_nn(in) > 0);
            if reg.flag_nonneg, reg.lambda_nonneg = lam_nn(in); end
            reg.lambda_tikhonov = lam_tik(it);
            Es(in,it) = fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask,b,beta,TE,dirs,CS_phase,D_FW, ...
                                      PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB,sel);
        end
    end
    [e,ix] = min(Es(:)); [ii,jj] = ind2sub(size(Es),ix);
    fprintf('  SNR %2d: lambda_nonneg = %-4g lambda_tikhonov = %-4g -> %.4f\n', ...
        SNRs(is), lam_nn(ii), lam_tik(jj), e);
end

% =========================================================================
%  Recommendation
% =========================================================================
fprintf('\n==================== RECOMMENDED ====================\n');
fprintf('options.fODF_regularization.flag_nonneg      = %d;\n', best.flag_nonneg);
fprintf('options.fODF_regularization.lambda_nonneg    = %g;\n', best.lambda_nonneg);
fprintf('options.fODF_regularization.tau              = %g;\n', best.tau);
fprintf('options.fODF_regularization.Lmax_init        = %d;\n', best.Lmax_init);
fprintf('options.fODF_regularization.lambda_tikhonov  = %g;\n', best.lambda_tikhonov);
fprintf('options.fODF_regularization.TikhonovMatrix   = ''%s'';\n', best.TikhonovMatrix);
fprintf('=====================================================\n');

% ---- Plots -------------------------------------------------------------
figure('Color','w','Name','fODF regularization sweep','Position',[80 80 1250 420])

subplot(1,3,1)
imagesc(E_odf), colorbar, hold on
plot(it_b,in_b,'w*','MarkerSize',12,'LineWidth',1.5)
set(gca,'XTick',1:numel(lam_tik),'XTickLabel',lam_tik, ...
        'YTick',1:numel(lam_nn), 'YTickLabel',lam_nn)
xlabel('\lambda_{tikhonov}'), ylabel('\lambda_{nonneg}  (0 = off)')
title('rel err fODF (lower is better)')

subplot(1,3,2), hold on
plot(taus,E_tau,'o-','LineWidth',1.4)
plot(best.tau,min(E_tau),'r*','MarkerSize',12)
xlabel('\tau'), ylabel('rel err fODF'), box on, grid on
title('tau, at the stage 1 optimum')

subplot(1,3,3), hold on
plot(lam_tik,E_mat(1,:),'o-','LineWidth',1.4)
plot(lam_tik,E_mat(2,:),'s-','LineWidth',1.4)
set(gca,'XScale','log')
xlabel('\lambda_{tikhonov}'), ylabel('rel err fODF'), box on, grid on
legend(mats,'Location','best'), title('Tikhonov matrix')
