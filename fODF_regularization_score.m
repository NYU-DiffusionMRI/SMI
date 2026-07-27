function [e_odf,e_plm,e_neg,e_pkamp,e_pkratio,e_peak,e_npk] = fODF_regularization_score(reg,dwi4,Lmax_shells,kernel,mask, ...
                                     b,beta,TE,dirs,CS_phase,D_FW,PLM_T,Y_odf,AMP_T,nrm_T,dirs_odf,NB,COND,FIB,sel)
% [e_odf,e_plm,e_neg,e_pkamp,e_pkratio,e_peak,e_npk] = fODF_regularization_score(reg,dwi4,...)
%
% The cheap scores come first on purpose: the peak ANGULAR error needs a
% per-voxel loop over local maxima, so it is only computed when the caller
% actually asks for a 6th output.
%
% Scores one fODF deconvolution regularization setting against a known
% ground truth. Used by example_fODF_regularization_sweep.m.
%
% Every voxel of the experiment is deconvolved in a single call, with the
% noise realizations stacked along the first image dimension. That is much
% faster than one call per realization and gives identical results.
%
% reg          options passed to SMI.get_plm_from_S_and_kernel
% dwi4         [Nvox 1 1 Ndwi] synthetic signals
% PLM_T        [Nvox x Nlm-1] ground truth plm of each voxel
% Y_odf        [Ndirs x Nlm] SH basis, already scaled by sqrt((2l+1)/(4pi))
% AMP_T,nrm_T  ground truth fODF amplitudes and their L2 norms
% dirs_odf,NB  scoring directions and their nearest neighbour indices
% COND         [Nvox x 2] condition index of each voxel [angle, SNR]
% FIB          cell array, FIB{ia} is the [Nfib x 3] true fiber directions
% sel          optional logical mask selecting a subset of the voxels
%
% e_odf     relative L2 error of the fODF over the sphere (primary score)
% e_plm     RMSE of the plm coefficients
% e_neg     fraction of the absolute fODF mass that is negative
% e_peak    mean angle [deg] from each true fiber to the closest fODF peak
% e_pkamp   mean peak (maximum over directions) fODF amplitude
% e_pkratio mean of peak(estimate)/peak(ground truth), per voxel. This is
%           the shrinkage: 1 means the regularizer preserves the height of
%           the fODF, below 1 means it flattens it. Regularization buys
%           accuracy partly by damping the coefficients, so this is the
%           cost side of that trade and is worth reporting next to e_odf.
% e_npk     mean number of detected peaks per voxel (2 is correct for the
%           two-fiber ground truth, above the angular resolution limit)
%
P = SMI.get_plm_from_S_and_kernel(dwi4,Lmax_shells,kernel,mask, ...
                                  b,beta,TE,dirs,CS_phase,D_FW,reg);
% reshape rather than squeeze: with a single voxel squeeze would collapse
% [1 1 1 Ncoeff] to a column and silently transpose the problem
P = reshape(P, [], size(P,4));
if nargin < 20 || isempty(sel), sel = true(size(P,1),1); end
P    = P(sel,:);
PT   = PLM_T(sel,:);
AT   = AMP_T(:,sel);
nT   = nrm_T(sel);
cond = COND(sel,:);

AMP = Y_odf(:,1) + Y_odf(:,2:end)*P.';                  % [Ndirs x Nvox]

e_odf = mean(sqrt(sum((AMP-AT).^2,1))./nT);
e_plm = sqrt(mean(sum((P-PT).^2,2)));
e_neg = mean(-sum(min(AMP,0),1)./sum(abs(AMP),1));

% Peak height and its shrinkage relative to the ground truth. Cheap, so
% always computed: this is the quantity that shows what the regularizer
% costs in fODF amplitude while it buys accuracy.
pk_est   = max(AMP,[],1);
pk_true  = max(AT ,[],1);
e_pkamp  = mean(pk_est);
e_pkratio = mean(pk_est./pk_true);

if nargout < 6, e_peak = NaN; e_npk = NaN; return, end

% Peak angular error. Peaks are local maxima above a quarter of the voxel
% maximum, which keeps the genuine lobes and drops the SH truncation
% ripples (a band limited fODF has a local maximum on every ripple).
A    = AMP.';                                           % [Nvox x Ndirs]
isPk = true(size(A));
for j = 1:size(NB,2)
    isPk = isPk & (A > A(:,NB(:,j)));
end
isPk = isPk & (A > max(0.25*max(A,[],2),0));

ang_err = nan(size(A,1),1);
for v = 1:size(A,1)
    pk = dirs_odf(isPk(v,:),:);
    if isempty(pk), continue, end
    fibv = FIB{cond(v,1)};
    d = acosd(min(abs(fibv*pk.'),1));                   % even SH are antipodal
    ang_err(v) = mean(min(d,[],2));
end
e_peak = mean(ang_err(~isnan(ang_err)));
e_npk  = mean(sum(isPk,2));
end
