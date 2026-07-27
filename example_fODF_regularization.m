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

% =========================================================================
% VISUALIZATION (single noise realization, the one stored in plm_examples)
%
% All the plots below show the fODF amplitude
%
%     fODF(u) = 1/(4*pi) + sum_{l>=2,m} plm_lm*Ylm(u)*sqrt((2l+1)/(4*pi))
%
% i.e. the same quantity the non-negativity constraint acts on. Negative
% amplitudes are physically meaningless, so they are always drawn in a way
% that makes them explicit rather than clipped away.
% =========================================================================
labels = cellfun(@strtrim,names,'UniformOutput',false);
cols = [0.85 0.33 0.10;   % unregularized
        0.00 0.45 0.74;   % tikhonov
        0.47 0.67 0.19;   % non-negativity
        0.49 0.18 0.56];  % nonneg+tikhonov

% ---- 2D: profile in the plane containing both fibers -------------------
phi = linspace(0,2*pi,721)';
dirs_plane = [cos(phi) sin(phi) zeros(size(phi))];
Y_plane = SMI.get_even_SH(dirs_plane,Lmax,CS_phase).*sqrt((2*L_all+1)/(4*pi));
amp = @(plm) Y_plane(:,1) + Y_plane(:,2:end)*plm(:);
ang = phi*180/pi;

figure('Color','w','Name','fODF in the fiber plane','Position',[80 80 1150 460])

% Signed amplitude: the clearest way to see the spurious negative lobes
subplot(1,2,1), hold on
plot(ang,amp(plm_true),'k','LineWidth',2)
for cc = 1:length(configs)
    plot(ang,amp(plm_examples(cc,:)),'LineWidth',1.3,'Color',cols(cc,:))
end
plot(ang,zeros(size(phi)),'k--')
ylims = get(gca,'YLim');
for kk = 1:size(fibers,1) % true fiber orientations (and their antipodes)
    a0 = atan2d(fibers(kk,2),fibers(kk,1));
    plot([a0 a0],ylims,':','Color',[.6 .6 .6])
    plot([a0 a0]+180,ylims,':','Color',[.6 .6 .6])
end
set(gca,'XTick',0:60:360), xlim([0 360]), ylim(ylims), box on
xlabel('angle in the fiber plane [degrees]'), ylabel('fODF amplitude')
legend(['ground truth',labels],'Location','best')
title('signed amplitude (dotted = true fiber directions)')

% Polar view: the fODF shape, with the negative part folded out separately
subplot(1,2,2), hold on
a_gt = amp(plm_true);
plot(max(a_gt,0).*cos(phi),max(a_gt,0).*sin(phi),'k','LineWidth',2)
for cc = 1:length(configs)
    a = amp(plm_examples(cc,:));
    plot(max(a,0).*cos(phi),max(a,0).*sin(phi),'LineWidth',1.3,'Color',cols(cc,:))
end
for cc = 1:length(configs) % |negative part|, dashed in the same colour
    rn = max(-amp(plm_examples(cc,:)),0);
    if any(rn>1e-12)
        plot(rn.*cos(phi),rn.*sin(phi),'--','LineWidth',0.9,'Color',cols(cc,:))
    end
end
axis equal, grid on, box on
xlabel('x'), ylabel('y')
title('polar view (solid = positive part, dashed = |negative part|)')

% ---- 3D: fODF glyphs ----------------------------------------------------
% Radius along u is the fODF amplitude, so the glyph is the familiar peanut
% shaped surface whose lobes point along the fiber directions
Ntheta = 61; Nphi = 121;
[TH,PH] = meshgrid(linspace(0,pi,Ntheta),linspace(0,2*pi,Nphi));
Xs = sin(TH).*cos(PH); Ys = sin(TH).*sin(PH); Zs = cos(TH);
Y_mesh = SMI.get_even_SH([Xs(:) Ys(:) Zs(:)],Lmax,CS_phase).*sqrt((2*L_all+1)/(4*pi));
amp_mesh = @(plm) reshape(Y_mesh(:,1)+Y_mesh(:,2:end)*plm(:),size(Xs));

amp_all = [{amp_mesh(plm_true)}, cell(1,length(configs))];
for cc = 1:length(configs), amp_all{cc+1} = amp_mesh(plm_examples(cc,:)); end
% Common radial scale so that the glyphs are directly comparable
rmax = max(cellfun(@(a) max(a(:)),amp_all));
titles = ['ground truth',labels];

figure('Color','w','Name','fODF glyphs','Position',[80 80 1200 760])
for kk = 1:length(amp_all)
    subplot(2,3,kk)
    SMI_plot_fODF_glyph(amp_all{kk},Xs,Ys,Zs,rmax,titles{kk},fibers)
end
subplot(2,3,6), axis off
text(0,0.5,{'Glyph radius = fODF amplitude,','all panels share the same scale.','', ...
            'Surface colour = orientation','(red = x, green = y, blue = z).','', ...
            'Translucent red lobes are the','NEGATIVE part of the fODF, drawn','at radius |fODF|.','', ...
            'Black lines are the two true','fiber directions.'},'FontSize',10,'VerticalAlignment','middle')

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

function SMI_plot_fODF_glyph(a,Xs,Ys,Zs,rmax,ttl,fibers)
% SMI_plot_fODF_glyph(a,Xs,Ys,Zs,rmax,ttl,fibers)
%
% Draws one fODF glyph on the current axes. The positive part of the
% amplitude is rendered as a surface with radius a(u) coloured by
% orientation, and the negative part (which is what the non-negativity
% constraint suppresses) is overlaid as a translucent red surface at radius
% |a(u)|, so that spurious lobes are visible instead of being clipped.
%
% a       [Ntheta x Nphi] fODF amplitude on the spherical mesh
% Xs,Ys,Zs  the unit sphere mesh the amplitude was evaluated on
% rmax    radial scale shared by all the glyphs that are compared
% fibers  [Nfibers x 3] true fiber directions, drawn as reference axes
%
rp = max(a,0);
surf(rp.*Xs,rp.*Ys,rp.*Zs,cat(3,abs(Xs),abs(Ys),abs(Zs)),...
     'EdgeColor','none','FaceColor','interp'), hold on
rn = max(-a,0);
if any(rn(:)>1e-12)
    surf(rn.*Xs,rn.*Ys,rn.*Zs,'EdgeColor','none','FaceColor',[0.85 0.1 0.1],'FaceAlpha',0.35)
end
for kk = 1:size(fibers,1)
    u = 1.05*rmax*fibers(kk,:);
    plot3([-u(1) u(1)],[-u(2) u(2)],[-u(3) u(3)],'k-','LineWidth',1)
end
axis(1.1*rmax*[-1 1 -1 1 -1 1]), daspect([1 1 1]), axis off
view(-25,30), camlight headlight, lighting gouraud, material dull
title(ttl,'FontWeight','normal')
end
