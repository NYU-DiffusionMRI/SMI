%% Removing b3000 shell and fitting DKI to SH signal
clc,clear,close all
root='/Volumes/cbi05data';
NII0=load_untouch_nii([root,'/data1/Hamster/Ben/atlas/SH/mean_b0.nii']);
NII1=load_untouch_nii([root,'/data1/Hamster/Ben/atlas/SH/mean_elastix_sh1.nii']);
NII2=load_untouch_nii([root,'/data1/Hamster/Ben/atlas/SH/mean_elastix_sh2.nii']);
nii_mask=load_untouch_nii([root,'/data1/Hamster/Ben/atlas/brainmask.nii']);
se = strel('cube',3);
mask = imerode(nii_mask.img,se);
% Normalize s_{lm} for S_l computation
b0=(NII0.img);
b1=(NII1.img)./(b0*sqrt(4*pi));
b2=(NII2.img)./(b0*sqrt(4*pi));
% Compute rotational invariants from s_lm (no need to generate directions and fit) 
S0=cat(4,b1(:,:,:,1),b2(:,:,:,1));
S2=sqrt(cat(4,sum(b1(:,:,:,2:6).^2,4),sum(b2(:,:,:,2:6).^2,4))/5);
S4=sqrt(cat(4,sum(b1(:,:,:,7:15).^2,4),sum(b2(:,:,:,7:15).^2,4))/9);
S6=sqrt(cat(4,sum(b1(:,:,:,16:28).^2,4),sum(b2(:,:,:,16:28).^2,4))/13);

sh1 = vec(NII1.img,logical(nii_mask.img));
sh2 = vec(NII2.img,logical(nii_mask.img));
sz=size(nii_mask.img);
Nb0s=4;
% Ndirs1=20;
% Ndirs2=60;
Ndirs1=40;
Ndirs2=100;
% Load precomputed directions for the shells
dirs1=load([root,'/data1/Hamster/Santi/MyMatlabTools/Directions_from_Koay/',num2str(Ndirs1),'_points_antipodally_symmetric.txt']);
dirs256=load([root,'/data1/Hamster/Santi/MyMatlabTools/Directions_from_Koay/',num2str(Ndirs2),'_points_antipodally_symmetric.txt']);
b1 = SHproj(sh1,dirs1,6);
b2 = SHproj(sh2,dirs256,6);
b1 = unvec(b1,logical(nii_mask.img));
b2 = unvec(b2,logical(nii_mask.img));

b0 = repmat(NII0.img,[1,1,1,Nb0s]);
b=[zeros(1,Nb0s) ones(1,Ndirs1) 2*ones(1,Ndirs2)];
Ndwi=length(b);
bvec=[repmat([0;0;0],1,Nb0s),dirs1',dirs256'];
bshape=0*b+1; % Only LTE data
B = ConstructAxiallySymmetricBset(b,bshape,bvec);
% These are the generated dwi (without noise, note that the signals are not normalized)
signal_SH = cat(4,b0,b1,b2);

se=strel('sphere',2);
emask=imerode(mask,se);

tic
grad=[bvec;b]';
[b0, dt] = dki_fit(signal_SH, grad, logical(emask), [0 0 0]);
[fa, md, rd, ad, fe, mk,  rk, ak] = dki_parameters(dt, logical(emask));
[~, ~, ~, ~, ~, mw,  rw, aw] = dwi_parameters(dt, logical(emask));
toc
IMGUI(cat(4,md, rd, ad,mw,rw,aw),[0 3])

Dlm=wrapper_STF_decomposition_backandforth(dt(:,:,:,1:6),'cart2STF',logical(mask),'JV');
Wlm=wrapper_STF_decomposition_backandforth(dt(:,:,:,7:21),'cart2STF',logical(mask),'JV');


% % % double-checking my functions
% % Dlm_back=wrapper_STF_decomposition_backandforth(dt(:,:,:,1:6),'cart2STF',logical(mask),'JV');
% % Wlm_back=wrapper_STF_decomposition_backandforth(dt(:,:,:,7:21),'cart2STF',logical(mask),'JV');
% % Dij_back=wrapper_STF_decomposition_backandforth(Dlm_back,'STF2cart',logical(mask),'SC');
% % Wijkl_back=wrapper_STF_decomposition_backandforth(Wlm_back,'STF2cart',logical(mask),'SC');
% % IMGUI(cat(4,dt(:,:,:,[1 4 6 2 3 5])-Dij_back,dt(:,:,:,[1 11 15 4 6 2 3 5 13 7 8 12 9 10 14]+6)-Wijkl_back),[-1 1]*1e-3)
% % tic
% % [~,Dlm,Wlm,~,~,~,~,~] = Fit_DKI_LTE(double(signal_SH),grad(:,4)',grad(:,1:3)', logical(mask), 0, 0);
% % toc
% % IMGUI(cat(4,Dlm-Dlm_back,Wlm-Wlm_back),[-1 1]*1e-3)

mask2=emask;mask2(:,:,[1:70 81:145])=0;
tic
[dt_rJV, b0_rJV, mkpred_rJV, list_rJV] = RobustDKIFitting(double(signal_SH), grad, logical(mask2));
dt_rJV(~isfinite(dt_rJV(:)))=0;
[~, ~, ~, ~, ~, mw_rJV,  rw_rJV, aw_rJV] = dwi_parameters(dt_rJV, logical(mask));
toc
IMGUI(cat(4, mw,mw_rJV, rw, rw_rJV,aw, aw_rJV),[0 3])

data.signal_SH=signal_SH;
data.gt.b0=b0;
data.gt.dt=dt;
data.gt.Dlm=Dlm;
data.gt.Wlm=Wlm;
data.grad=grad;
data.mask=mask;
data.mask2=mask2;
% save('/Volumes/labspace/Santiago/MyRobustDKI/BensPhantom_DKI_v3.mat','-struct','data')
save('/Volumes/labspace/Santiago/MyRobustDKI/BensPhantom_DKI_v4.mat','-struct','data')

%% Generating DKI phantom from fits and tuning it to avoid black voxels
clc,clear,close all
load('/Volumes/labspace/Santiago/MyRobustDKI/BensPhantom_DKI_v3.mat')
b0=gt.b0;
dt=gt.dt;
Dlm=gt.Dlm;
Wlm=gt.Wlm;
clearvars gt
dt(~isfinite(dt(:)))=0;
signal_dki = dki_signal([],dt,grad,logical(mask),'JV');
[~, ~, ~, ~, ~, mk,  rk, ak] = dki_parameters(dt, logical(mask));
[fa, md, rd, ad, ~, mw,  rw, aw] = dwi_parameters(dt, logical(mask));
% tic
% [b02,Dlm2,Wlm2,Dl2,Wl2,DTI_scalars,DKI_scalars,ExtraScalars] = Fit_DKI_LTE(double(signal_SH),grad(:,4)',grad(:,1:3)', logical(mask2), 1, 0);
% toc
% IMGUI(cat(4,md,DTI_scalars.md,rd,ExtraScalars.rd,ad,ExtraScalars.ad,mw,DKI_scalars.mw,rw,ExtraScalars.rw,aw,ExtraScalars.aw),[0 3])
% IMGUI(cat(4,mk,ExtraScalars.mk,rk,ExtraScalars.rk,ak,ExtraScalars.ak),[0 3])

D0 = Dlm(:,:,:,1)/sqrt(4*pi);
D2=sqrt(sum(Dlm(:,:,:,2:6).^2,4))/sqrt(5*4*pi);
W0 = Wlm(:,:,:,1)/sqrt(4*pi);
W2=sqrt(sum(Wlm(:,:,:,2:6).^2,4))/sqrt(5*4*pi);
W4=sqrt(sum(Wlm(:,:,:,7:15).^2,4))/sqrt(9*4*pi);
aw_axsym=(W0+5*W2+9*W4).*D0.^2./(D0+5*D2).^2;
rw_axsym=(W0-5/2*W2+27/8*W4).*D0.^2./(D0-5/2*D2).^2;
IMGUI(cat(4,mw,rw,rw_axsym,aw,aw_axsym),[0 3])

nametags={'$W_0$','$W_2$','$W_4$','MW','RW','AW'};
clims=[0 2.5;0 0.6;0 0.1;0 3;0 3;0 3];
figure('Position',[42 128 2364 1062]), colormap gray
WrapperPlotManySlices(cat(4,W0,W2,W4,mw,  rw, aw), 67,clims,nametags,2,[]),
sgt = sgtitle('raw $W_\ell$ and W maps from WLLS fit on Spherical Harmonics signals', 'Interpreter', 'latex'); sgt.FontSize = 30;

[iD,iK] = get_Iso_K_4Darray(signal_SH,grad(:,4),logical(mask2));


%% Remove black voxels and run training
clc,close all
se=strel('sphere',3); 
emask=imerode(mask,se);
% mask_test=emask;
mask3=emask;mask3(:,:,[1:60 86:145])=0;
mask_test=mask3;

% Grab ugly voxels
Dlm_2D=vectorize_JV(Dlm,logical(mask_test));
Wlm_hat_2D=vectorize_JV(Wlm.*md.^2,logical(mask_test));
dirs256 = get256dirs_JV();
Y_LM_matrix2 = get_even_SH(dirs256,4,0);
Y00=Y_LM_matrix2(:,1);
Y2m_2=Y_LM_matrix2(:,2:6);
Y4m_2=Y_LM_matrix2(:,7:15);
adc = [Y00 Y2m_2]*Dlm_2D;
akc = [Y00 Y2m_2 Y4m_2]*Wlm_hat_2D./(adc.^2);
neg_Kn=any(akc<0);
neg_Dn=any(adc<0);
neg_iK=vectorize_JV(iK,logical(mask_test))<0;
W0_bounded=vectorize_JV(W0,logical(mask_test))<3;
idx= ~neg_Kn & ~neg_Dn & ~neg_iK & W0_bounded;

% Doing the PR on the svd of the dwi
signal_dki_vect=vectorize_JV(signal_dki,logical(mask_test));
Wl_train=vectorize_JV(cat(4,W0,W2,W4),logical(mask_test));
trainingSize=5e4;
id_rand_training=randperm(size(signal_dki_vect,2),trainingSize);
tic
[coeff, score] = pca(signal_dki_vect(:,id_rand_training)');
time_pca=toc;
% s=svd(signal_dki_vect(:,id_rand)');
% semilogy(s,'-o')
fprintf('pca on training computed in %f s\n', time_pca);

npars=21;
reducedDimension = coeff(:,1:npars);
data_training = signal_dki_vect(:,id_rand_training)' * reducedDimension;
% =========================================================================
% Applying training transformation to test dataset (optimal)
data_testing = signal_dki_vect' * reducedDimension;
% =========================================================================
% Polynomial Regression FITTING using svd
Degree=2;
X_train = Compute_extended_moments(data_training,Degree);
X_fit_pca = Compute_extended_moments(data_testing,Degree);
tic
Pinv_X=pinv(X_train);
time_pinv=toc;
fprintf('pseudoinversion computed in %f s\n', time_pinv);
coeffs_W0=Pinv_X*Wl_train(1,id_rand_training)';
coeffs_W2=Pinv_X*Wl_train(2,id_rand_training)';
coeffs_W4=Pinv_X*Wl_train(3,id_rand_training)';
fit_W0_trainingdata=X_train*coeffs_W0;
fit_W2_trainingdata=X_train*coeffs_W2;
fit_W4_trainingdata=X_train*coeffs_W4;
clearvars pinvX

figure('Position',[197 462 2035 584])
subplot(131), plot(Wl_train(1,id_rand_training)',fit_W0_trainingdata,'.',Wl_train(1,id_rand_training)',Wl_train(1,id_rand_training)','r-.'), set(gca,'FontSize',20), title('$W_0$','interpreter','latex'), xlabel('dki fit','interpreter','latex'), ylabel('RotInv fit','interpreter','latex'), axis([0 4 0 4])
subplot(132), plot(Wl_train(2,id_rand_training)',fit_W2_trainingdata,'.',Wl_train(2,id_rand_training)',Wl_train(2,id_rand_training)','r-.'), set(gca,'FontSize',20), title('$W_2$','interpreter','latex'), xlabel('dki fit','interpreter','latex'), ylabel('RotInv fit','interpreter','latex'), axis([0 1 0 1])
subplot(133), plot(Wl_train(3,id_rand_training)',fit_W4_trainingdata,'.',Wl_train(3,id_rand_training)',Wl_train(3,id_rand_training)','r-.'), set(gca,'FontSize',20), title('$W_4$','interpreter','latex'), xlabel('dki fit','interpreter','latex'), ylabel('RotInv fit','interpreter','latex'), axis([0 0.15 0 0.15])

% signal_dki_sm=signal_dki;
% for ii=1:size(signal_dki_sm,4)
%     signal_dki_sm(:,:,:,ii) = smooth3(signal_dki(:,:,:,ii), 'gaussian', [3 3 3],0.8);
% end
% IMGUI(signal_dki,[0 1])
% IMGUI(signal_dki_sm,[0 1])
% signal_dki_vect=vectorize_JV(signal_dki_sm,logical(mask));
signal_dki_vect=vectorize_JV(signal_dki,logical(mask));
data_testing = signal_dki_vect' * reducedDimension;
X_fit_pca = Compute_extended_moments(data_testing,Degree);
fit_W0=vectorize_JV((X_fit_pca*coeffs_W0)',logical(mask)); fit_W0(fit_W0(:)<0)=0;
fit_W2=vectorize_JV((X_fit_pca*coeffs_W2)',logical(mask)); fit_W2(fit_W2(:)<0)=0;
fit_W4=vectorize_JV((X_fit_pca*coeffs_W4)',logical(mask)); fit_W4(fit_W4(:)<0)=0;

% % Using W0,W2,W4 from PR
% aw_axsym_PR=(fit_W0+5*fit_W2+9*fit_W4).*D0.^2./(D0+5*D2).^2;
% rw_axsym_PR=(fit_W0-5/2*fit_W2+27/8*fit_W4).*D0.^2./(D0-5/2*D2).^2;

% % Using W0,W2 from PR and W4 from WLLS
% aw_axsym_PR=(fit_W0+5*fit_W2+9*W4).*D0.^2./(D0+5*D2).^2;
% rw_axsym_PR=(fit_W0-5/2*fit_W2+27/8*W4).*D0.^2./(D0-5/2*D2).^2;

% Using W0 from PR and W2,W4 from WLLS
aw_axsym_PR=(fit_W0+5*W2+9*W4).*D0.^2./(D0+5*D2).^2;
rw_axsym_PR=(fit_W0-5/2*W2+27/8*W4).*D0.^2./(D0-5/2*D2).^2;


IMGUI(cat(4,W0/2,fit_W0/2,W2/0.6,fit_W2/0.6,W4/0.1,fit_W4/0.1,rw_axsym/3,rw_axsym_PR/3,aw_axsym/3,aw_axsym_PR/3),[0 1])





rw_axsym_PR_tuned=rw_axsym_PR;
rw_axsym_PR_tuned(rw_axsym_PR_tuned(:)<0.2)=3;
% IMGUI(cat(4,rw_axsym_PR,rw_axsym_PR_tuned),[0 3])

nametags={'MW','RW','AW','MW','RW','AW'};
clims=[0 3;0 3;0 3;0 3;0 3;0 3];
figure('Position',[42 128 2364 1062]), colormap gray
WrapperPlotManySlices(cat(4,mw,  rw, aw,fit_W0,rw_axsym_PR,aw_axsym_PR), 67,clims,nametags,2,[]),
sgt = sgtitle('raw $W_\ell$ and W maps from WLLS fit on Spherical Harmonics signals', 'Interpreter', 'latex'); sgt.FontSize = 30;

rw_axsym_PR_nan=rw_axsym_PR;
rw_axsym_PR_nan(rw_axsym_PR_nan(:)<1)=nan;
IMGUI(cat(4,rw_axsym_PR,rw_axsym_PR_nan),[0 3])
tic
[rw_axsym_PR_filt] = nonlocalFilter(rw_axsym_PR_nan, [5 5 5], 70, 'median');
toc
load('/Volumes/labspace/Santiago/MyRobustDKI/BensPhantom_DKI_ROIbb.mat');
rw_axsym_PR_tuned=rw_axsym_PR;
rw_axsym_PR_tuned(ROI(:)==1)=rw_axsym_PR_filt(ROI(:)==1);
IMGUI(cat(4,rw_axsym_PR,rw_axsym_PR_tuned),[0 3])
% tic
% [rw_axsym_PR_tuned_v2] = nonlocalFilter(rw_axsym_PR_tuned, [5 5 5], 70, 'mean');
% toc
% IMGUI(cat(4,rw_axsym_PR,rw_axsym_PR_tuned,rw_axsym_PR_tuned_v2),[0 3])
rw_axsym_PR_tuned(rw_axsym_PR_tuned(:)<0)=eps;


% Fix W0
a1=(D0+5*D2).^2./D0.^2.*aw_axsym_PR-fit_W0;
b1=(D0-5/2*D2).^2./D0.^2.*rw_axsym_PR_tuned-fit_W0;
W0_back=fit_W0;
W4_back=4/63*(a1+2*b1);
W4_back(W4_back(:)<0)=0;
W2_back=1/5*(a1-9*W4_back);
W2_back(W2_back(:)<0)=0;

IMGUI(cat(4,W0/2,W0_back/2,W2/0.6,W2_back/0.6,W4/0.1,W4_back/0.1),[0 1])

aw_back=(W0_back+5*W2_back+9*W4_back).*D0.^2./(D0+5*D2).^2;
rw_back=(W0_back-5/2*W2_back+27/8*W4_back).*D0.^2./(D0-5/2*D2).^2;
[sum(W0_back(:)<0) sum(rw_back(:)<0) sum(aw_back(:)<0)]


DKImaps=cat(4,W0,rw_axsym,aw_axsym,W0_back,rw_back,aw_back);nametags={'MW','RW','AW','MW','RW','AW','MW','RW','AW'};
clims=[0 3;0 3;0 3;0 3;0 3;0 3;0 3;0 3;0 3];
figure('Position',[182 310 2056 853]), colormap gray
WrapperPlotManySlices(DKImaps, 67,clims,nametags,2,[]),
sgt = sgtitle('DKI phantom: original and tuned', 'Interpreter', 'latex'); sgt.FontSize = 30;

DKImaps=cat(4,W0,W2,W4,W0_back,W2_back,W4_back);nametags={'$W_0$','$W_2$','$W_4$','$W_0$','$W_2$','$W_4$'};
clims=[0 2.5;0 0.6;0 0.1;0 2.5;0 0.6;0 0.1];
figure('Position',[182 310 2056 853]), colormap gray
WrapperPlotManySlices(DKImaps, 67,clims,nametags,2,[]),
sgt = sgtitle('DKI phantom: original and tuned', 'Interpreter', 'latex'); sgt.FontSize = 30;


%% Save readjusted W tensor
clc,close all
IMGUI(cat(4,W0_back./W0,W2_back./W2,W4_back./W4),[0 2])

scaling=cat(4,W0_back./W0,repmat(W2_back./W2,[1 1 1 5]),repmat(W4_back./W4,[1 1 1 9]));
Wlm_rescaled=Wlm.*scaling;


tic
maps = get_maps_from_4D_DW_tensors(cat(4,Dlm,Wlm_rescaled),logical(mask2),'STF',[],1);
toc
% IMGUI(cat(4,maps.fa,maps.md/2,maps.rd/2,maps.rd_axsym/2,maps.ad/2,maps.ad_axsym/2,maps.mw/2.5,maps.mk/2.5,maps.rw_axsym/3,maps.rw/3,maps.rk/3,maps.aw_axsym/3,maps.aw/3,maps.ak/3),[0 1])
% IMGUI(cat(4,maps.mw/2.5,maps.mk/2.5,maps.rw_axsym/3,maps.rw/3,maps.rk/3,maps.aw_axsym/3,maps.aw/3,maps.ak/3),[0 1])
% IMGUI(cat(4,maps.rw-maps.rk,maps.aw-maps.ak),[-0.01 0.01])


Dij_final=wrapper_STF_decomposition_backandforth(Dlm,'STF2cart',logical(mask),'SC');
Wijkl_final=wrapper_STF_decomposition_backandforth(Wlm_rescaled,'STF2cart',logical(mask),'SC');
dwi_synthetic = dki_signal([],cat(4,Dij_final,Wijkl_final),grad,logical(mask),'SC');
[min(dwi_synthetic(:)) max(dwi_synthetic(:))]

aaa=vectorize_JV(dwi_synthetic,logical(emask));
[min(aaa(:)) max(aaa(:))]
aaaa=sort(aaa,'descend');


% IMGUI(Wlm_rescaled)
% dwi_synthetic = dki_signal(gt.b0,gt.dt,grad,mask);
% dwi_synthetic = dki_signal(gt.b0,cat(4,dki.Dij,dki.Wijkl),grad,mask);
% IMGUI(dwi_synthetic-signal_dki,[-1 1]*1e-8)

tic
maps1 = get_maps_from_4D_DW_tensors(cat(4,dki.Dlm,dki.Wlm),logical(mask),'STF',[],0);
toc

gt.fa=maps.fa;
gt.md=maps.md;
gt.rd_axsym=maps.rd_axsym;
gt.ad_axsym=maps.ad_axsym;
gt.mw=maps.mw;
gt.rw_axsym=maps.rw_axsym;
gt.aw_axsym=maps.aw_axsym;

% dki.dwi=single(dwi_synthetic);
B=ConstructAxiallySymmetricBset(grad(:,4)',1+0*grad(:,4)',grad(:,1:3)');
dki.B=B;
dki.grad=grad;
dki.mask=logical(mask);
dki.b0=b0;
dki.Dlm=Dlm;
dki.Wlm=Wlm_rescaled;
dki.Dij=Dij_final;
dki.Wijkl=Wijkl_final;
dki.gt=gt;
save('/Volumes/labspace/Santiago/MyRobustDKI/BensPhantom_DKI_v2_tuned.mat','-struct','dki')

tic
maps1 = get_maps_from_4D_DW_tensors(cat(4,dki.Dlm,dki.Wlm),logical(mask),'STF',[],0);
toc
tic
maps2 = get_maps_from_4D_DW_tensors(cat(4,dki.Dij,dki.Wijkl),logical(mask),'cart',[],0);
toc

IMGUI(cat(4,maps1.fa-maps2.fa,maps1.md-maps2.md,maps1.mw-maps2.mw,maps1.rw_axsym-maps2.rw_axsym,maps1.aw_axsym-maps2.aw_axsym),[-1 1]*1e-6)




%%
tic
% maps = get_maps_from_4D_DW_tensors(dt,logical(mask2),'cart','JV',1);
maps = get_maps_from_4D_DW_tensors(cat(4,Dlm,Wlm),logical(mask2),'STF',[],1);
toc
% IMGUI(cat(4,maps.fa,maps.md/2,maps.rd/2,maps.rd_axsym/2,maps.ad/2,maps.ad_axsym/2,maps.mw/2.5,maps.mk/2.5,maps.rw_axsym/3,maps.rw/3,maps.rk/3,maps.aw_axsym/3,maps.aw/3,maps.ak/3),[0 1])
% IMGUI(cat(4,maps.mw/2.5,maps.mk/2.5,maps.rw_axsym/3,maps.rw/3,maps.rk/3,maps.aw_axsym/3,maps.aw/3,maps.ak/3),[0 1])
% IMGUI(cat(4,maps.rw-maps.rk,maps.aw-maps.ak),[-0.01 0.01])


%%
clc,close all
colormap gray
slices=62:79;
figure('Position',[88 198 2192 971]), colormap gray
for ii=1:length(slices)
%     imagesc(rw_axsym_PR(:,:,slices(ii)),[0 3]), title(['slice = ',num2str(slices(ii))]), set(gca,'FontSize',20)
%     set(gcf,'Position',[529 98 1290 1120])    
%     BW1=roipoly;  
%     BW2=roipoly;  
%     BW(:,:,ii)=BW1|BW2;
    subplot(121),     imagesc(rw_axsym_PR(:,:,slices(ii)),[0 3])
    subplot(122),     imagesc(rw_axsym_PR(:,:,slices(ii)).*(~BW(:,:,ii)),[0 3])
    pause(0.2)
end
ROI=0*mask;
ROI(:,:,slices)=BW;
IMGUI(ROI&mask,[0 1])
% save('/Volumes/labspace/Santiago/MyRobustDKI/BensPhantom_DKI_ROIbb.mat','ROI')

%% smooth gt
clc,close all
rw_smooth=rw_axsym_PR;
rw_smooth(rw_smooth(:)<1)=nan;
% IMGUI(cat(4,rw_axsym_PR,rw_smooth),[0 3])
t = 6; f = 3; h1 = 1; h2 = .1; selfsim = 0;
% t = 3; f = 1; h1 = 1/6; h2 = .1/6; selfsim = 0;
slices=77;%62:79;
for ii=1:length(slices)
    rw_smooth(:,:,slices(ii)) = simple_nlm(double(rw_smooth(:,:,slices(ii))),t,f,h1,h2,selfsim);
end
rw_smooth(ROI(:)==0)=rw_axsym_PR(ROI(:)==0);
% IMGUI(cat(4,rw_axsym_PR,rw_smooth),[0 3])

figure('Position',[88 198 2192 971]), colormap gray
subplot(121),     imagesc(rw_axsym_PR(:,:,77),[0 3])
subplot(122),     imagesc(rw_smooth(:,:,77),[0 3])


rw_ROI=rw_axsym_PR(ROI(:)==1);
histogram(rw_ROI,linspace(-1,4,20))

img=rw_axsym_PR(:,:,77);
img(img(:)<0)=0;
img(img(:)>3)=3;
img=uint8(img*85);
imagesc(img)
imwrite(img,'/Volumes/labspace/Santiago/MyRobustDKI/test_slice_rw.png')

a=imread('/Volumes/labspace/Santiago/MyRobustDKI/test_slice_rw.png');



%%
% % % % %% smooth gt
% % % % clc,close all
% % % % Wl_smooth=Wl_back;
% % % % t = 3; f = 1; h1 = 1/6; h2 = .1/6; selfsim = 0;
% % % % for ii=61:89%1:size(W0_back,3)
% % % % %     Wl_smooth(:,:,ii,1) = simple_nlm(double(W0_back(:,:,ii)),t,f,h1,h2,selfsim);
% % % % %     Wl_smooth(:,:,ii,2) = simple_nlm(double(W2_back(:,:,ii)),t,f,h1,h2,selfsim);
% % % %     Wl_smooth(:,:,ii,3) = simple_nlm(double(W4_back(:,:,ii)),t,f,h1,h2,selfsim);
% % % % end
% % % % 
% % % % IMGUI(cat(4,Wl_back(:,:,:,1)/2.5,Wl_back(:,:,:,2)/0.6,Wl_back(:,:,:,3)/0.1,Wl_smooth(:,:,:,3)/0.1),[0 1])
