% example_fODF_modulation.m
%
% Modulating the SMI fODF by orientational coherence (p2), so that its
% amplitude separates coherent tissue from grey matter and CSF the way an
% MRtrix FOD does -- without using a tissue type criterion, which would delete
% peritumoral edema, the region of interest.
%
% Full methodology and results: REPORT_fODF_modulation.md
%
% THE PROBLEM
%
% out.plm is stored in the normalized convention p_00 = 1, so the fODF
% integrates to 1 in EVERY voxel and its isotropic floor is a fixed
% 1/(4*pi) = 0.0796. MRtrix's default iFOD2 termination threshold is
% -cutoff 0.05. The floor is above the threshold, so an unmodulated SMI fODF
% passes MRtrix's termination test everywhere in the brain, CSF included, and
% the tracker loses its ability to tell tissue from non-tissue by amplitude.
%
% THE EXPERIMENT
%
% Seven voxel classes, 200 noise realizations each, run through SMI.fit twice
% (unregularized and regularized deconvolution) at two SNRs. The third class,
% white matter inside edema, is the one that must survive: its axonal fraction
% and free water fraction are edema-like, but its fibres stay coherent. A
% weight that suppresses GM and CSF while keeping that class is doing what is
% wanted; a weight that suppresses it too is a tissue type criterion in
% disguise and must be rejected.
%
% Helper functions live in fODF_modulation_helpers.m rather than at the end of
% this script, so that this example runs under Octave as well as MATLAB.
%
% Runtime is a few minutes per arm. Set QUICK = 1 for a smaller version.

clear; close all
QUICK = 0;

H = fODF_modulation_helpers();

CS_phase = 1; D_FW = 3; Lmax_fod = 4;
Nreal = 200; if QUICK, Nreal = 20; end
MRTRIX_CUTOFF = 0.05;

% ---- protocol: 6 shells, Lmax = [0 2 2 2 4 4] ----------------------------
bvals = [0 1 2 3 4.5 6];
Ndirs = [8 30 30 60 60 60];
Lmax_shell = [0 2 2 2 4 4];
b = []; dirs = [];
for k = 1:numel(bvals)
    dirs = [dirs; H.dirs(Ndirs(k))];
    b    = [b; bvals(k)*ones(Ndirs(k),1)];
end
beta = ones(size(b)); TE = zeros(size(b)); Ndwi = numel(b);

K_WM  = [0.60 2.0 2.0 0.50 0.05 1 1];
K_ED  = [0.30 2.0 2.2 1.00 0.50 1 1];   % edema: low f, high fw, coherent fibres
K_GM  = [0.15 1.2 1.0 0.80 0.10 1 1];
K_CSF = [0.02 2.0 3.0 3.00 0.95 1 1];

%          name                kernelA kapA nfibA  kernelB kapB  fraction of A
tissue = { {'WM single',        K_WM,  16,  1,  [],    0,   1.0}, ...
           {'WM crossing 60',   K_WM,  16,  2,  [],    0,   1.0}, ...
           {'WM in edema',      K_ED,  12,  1,  [],    0,   1.0}, ...
           {'GM',               K_GM,   0.8,1,  [],    0,   1.0}, ...
           {'CSF',              K_CSF,  0,  1,  [],    0,   1.0}, ...
           {'WM/CSF interface', K_WM,  16,  1,  K_CSF, 0,   0.5}, ...
           {'WM/GM interface',  K_WM,  16,  1,  K_GM,  0.8, 0.5} };
Nclass = numel(tissue);
names  = cellfun(@(t) t{1}, tissue, 'UniformOutput', false);

dq   = H.dirs(3000);                       % quadrature grid
de   = H.dirs(2000); de = de(de(:,3)>=0,:); % evaluation grid (hemisphere)
Nlm  = Lmax_fod*(Lmax_fod+3)/2;
Nvox = Nreal*Nclass;
% A genuinely 3D volume: SMI.vectorize takes a different code path if any
% spatial dimension collapses to a singleton.
vol  = [Nreal/2 Nclass 2];

SNRs = [30 15];
store = struct();

for isnr = 1:numel(SNRs)
  SNR = SNRs(isnr);
  try, rng(7); catch, rand('seed',7); randn('seed',7); end   % rng is MATLAB only
  S_flat   = zeros(Ndwi,Nvox);
  class_id = zeros(1,Nvox);
  i = 0;
  for c = 1:Nclass
    kA = tissue{c}{2}; kapA = tissue{c}{3}; nfibA = tissue{c}{4};
    kB = tissue{c}{5}; kapB = tissue{c}{6}; vA    = tissue{c}{7};
    for r = 1:Nreal
        i = i+1;
        n1 = randn(1,3); n1 = n1/norm(n1);
        fod = H.watson(dq,n1,kapA);
        if nfibA == 2
            t = randn(1,3); t = t-(t*n1')*n1; t = t/norm(t);
            fod = fod + H.watson(dq, cos(pi/3)*n1+sin(pi/3)*t, kapA);
        end
        pA = H.mixture_plm(fod,dq,Lmax_fod,CS_phase);
        S  = vA*H.signal(pA,kA,b,beta,TE,dirs,Lmax_fod,CS_phase,D_FW);
        if ~isempty(kB)
            m1 = randn(1,3); m1 = m1/norm(m1);
            pB = H.watson_plm(dq,m1,kapB,Lmax_fod,CS_phase);
            S  = S + (1-vA)*H.signal(pB,kB,b,beta,TE,dirs,Lmax_fod,CS_phase,D_FW);
        end
        % Rician noise; S0 = 1 by construction
        S_flat(:,i) = sqrt((S+randn(Ndwi,1)/SNR).^2 + (randn(Ndwi,1)/SNR).^2);
        class_id(i) = c;
    end
  end

  opts = struct('b',b,'dirs',dirs,'beta',beta,'TE',TE, ...
                'sigma',ones(vol)/SNR,'mask',true(vol),'CS_phase',CS_phase, ...
                'D_FW',D_FW,'Lmax',Lmax_shell,'flag_fit_fODF',1);
  opts.compartments = {'IAS','EAS','FW'};

  arms = {{'unregularized', struct('flag_nonneg',0,'lambda_tikhonov',0)}, ...
          {'regularized',   struct('flag_nonneg',1,'lambda_nonneg',10,'lambda_tikhonov',0.3)}};
  for arm = 1:numel(arms)
      o = opts;
      o.fODF_regularization = arms{arm}{2};
      o.filename_log = sprintf('log_modulation_SNR%d_%s.txt',SNR,arms{arm}{1});
      fprintf('fitting SNR %d, %s deconvolution ...\n', SNR, arms{arm}{1});
      store(isnr,arm).name     = arms{arm}{1};
      store(isnr,arm).SNR      = SNR;
      store(isnr,arm).out      = SMI.fit(reshape(S_flat',[vol Ndwi]), o);
      store(isnr,arm).class_id = class_id;
  end
end

%% ------------------------------------------------------------------------
% Scoring. The question tractography actually asks is how many GM and CSF
% voxels stay above the termination threshold.
for isnr = 1:numel(SNRs)
 for arm = 1:2
    out      = store(isnr,arm).out;
    class_id = store(isnr,arm).class_id;

    W = {'none',        ones(vol); ...
         'kernel p2',   SMI.fODF_ModulationWeight(out,struct('source','kernel_p2')); ...
         'kernel p2^2', SMI.fODF_ModulationWeight(out,struct('source','kernel_p2','exponent',2)); ...
         'pl2',         SMI.fODF_ModulationWeight(out,struct('source','pl2')); ...
         'pl4',         SMI.fODF_ModulationWeight(out,struct('source','pl4')); ...
         'p2product',   SMI.fODF_ModulationWeight(out,struct('source','p2product'))};

    plm = reshape(out.plm,[Nvox size(out.plm,4)])';
    pk  = zeros(1,Nvox);
    for i = 1:Nvox, pk(i) = H.peak(plm(:,i),de,Lmax_fod,CS_phase); end

    fprintf('\n=== SNR %d, %s: %% of voxels above the MRtrix cutoff %.2f ===\n', ...
            SNRs(isnr), store(isnr,arm).name, MRTRIX_CUTOFF);
    fprintf('%-19s',''); fprintf('%13s',W{:,1}); fprintf('\n');
    for c = 1:Nclass
        m = class_id==c;
        fprintf('%-19s',names{c});
        for k = 1:size(W,1)
            wk = W{k,2}(:)';
            fprintf('%12.1f%%',100*mean(wk(m).*pk(m) > MRTRIX_CUTOFF));
        end
        fprintf('\n');
    end
 end
end

%% ------------------------------------------------------------------------
% Applying the modulation. The result is the SH coefficient array INCLUDING
% the l=0 term, which out.plm does not store: without it the density weighting
% would be lost.
%
% Either post hoc on an already fitted out ...
out = store(end,2).out;
[sh,w,info] = SMI.modulate_fODF(out, struct('source','p2product','mode','density'));
fprintf(['\nmodulated fODF: %d SH coefficients, weight in [%.3f %.3f], ' ...
         '%d degenerate voxels, mode ''%s''\n'], ...
        size(sh,4), min(w(:)), max(w(:)), info.Ndegenerate, info.mode);

% ... or as a flag on the fit itself, which is the way to run one instance
% with the modulation and one without. out.plm, out.pl and out.kernel are
% identical either way; only out.fODF_modulated is added.
o = opts;
o.fODF_regularization = arms{2}{2};
o.filename_log        = 'log_modulation_flag_on.txt';
o.fODF_modulation     = struct('flag_modulate',1);   % default 0, off
out_on = SMI.fit(reshape(S_flat',[vol Ndwi]), o);
fprintf('flag on: out.fODF_modulated is %s, out.plm unchanged: %d\n', ...
        mat2str(size(out_on.fODF_modulated)), ...
        isequal(out_on.plm, out.plm));
fprintf('NOTE: in density mode p_00 is no longer 1. Record which convention\n');
fprintf('      a saved NIfTI is in, and check SMI''s SH basis against MRtrix''s\n');
fprintf('      ordering and Condon-Shortley convention before running tckgen.\n');
