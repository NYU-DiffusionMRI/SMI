function H = fODF_modulation_helpers()
% H = fODF_modulation_helpers()
%
% Returns a struct of function handles used by example_fODF_modulation.m:
%
%   H.dirs(N)                              N nearly uniform directions on the sphere
%   H.watson_plm(dirs_q,n,kappa,Lmax,CS)   normalized plm of a Watson fODF about n
%   H.mixture_plm(fodf,dirs_q,Lmax,CS)     normalized plm of an arbitrary sampled fODF
%   H.watson(dirs_q,n,kappa)               unnormalized Watson amplitudes
%   H.signal(plm,kernel,b,beta,TE,dirs,Lmax,CS,D_FW)   forward DWI signal
%   H.peak(plm,dirs_e,Lmax,CS)             fODF peak amplitude and direction
%
% These live in their own file, rather than as local functions at the end of
% the example script, because MATLAB requires local functions at the END of a
% script while Octave cannot call them there at all. Subfunctions of a function
% file work in both, which keeps the example runnable and testable outside
% MATLAB. Same reasoning as fODF_regularization_score.m.
H = struct();
H.dirs        = @uniform_sphere_dirs;
H.watson      = @watson_amp;
H.watson_plm  = @watson_plm;
H.mixture_plm = @fodf_to_plm;
H.signal      = @forward_signal;
H.peak        = @fodf_peak;
end

% =====================================================================
function dirs = uniform_sphere_dirs(N)
% Golden angle spiral over the full sphere, usable as an equal area
% quadrature grid.
ga = pi*(3-sqrt(5));
k  = (0:N-1)';
z  = 1 - 2*(k+0.5)/N;
r  = sqrt(max(0,1-z.^2));
th = ga*k;
dirs = [r.*cos(th), r.*sin(th), z];
end

% =====================================================================
function f = watson_amp(dirs_q,n,kappa)
% Watson distribution about axis n, unnormalized. kappa <= 0 gives isotropic.
if kappa <= 0
    f = ones(size(dirs_q,1),1);
else
    f = exp(kappa*((dirs_q*n(:)/norm(n)).^2));
end
end

% =====================================================================
function plm = watson_plm(dirs_q,n,kappa,Lmax,CS_phase)
plm = fodf_to_plm(watson_amp(dirs_q,n,kappa),dirs_q,Lmax,CS_phase);
end

% =====================================================================
function plm = fodf_to_plm(fodf,dirs_q,Lmax,CS_phase)
% Project a sampled fODF onto normalized plm, in SMI's p_00 = 1 convention:
%   f_lm = plm .* sqrt((2l+1)/(4*pi)),  f_00 = 1/sqrt(4*pi)
% The returned vector covers l = 2..Lmax, matching out.plm.
Y  = SMI.get_even_SH(dirs_q,Lmax,CS_phase);
wq = 4*pi/size(dirs_q,1);                       % equal area quadrature weight
fodf = fodf(:)/(sum(fodf(:))*wq);               % so that \int fODF dOmega = 1
flm  = (Y'*fodf)*wq;
L_all   = repelem(0:2:Lmax, 2*(0:2:Lmax)+1)';
plm_all = flm./sqrt((2*L_all+1)/(4*pi));
plm     = plm_all(2:end)/plm_all(1);
end

% =====================================================================
function S = forward_signal(plm,kernel_vec,b,beta,TE,dirs,Lmax,CS_phase,D_FW)
% Forward model, built from exactly the construction that
% SMI.get_plm_from_S_and_kernel inverts:
%   S(u_i)/S0 = sum_{lm} K_l(b_i) * p_lm * Y_lm(u_i) * sqrt((2l+1)*4pi)
% kernel_vec is [f Da Depar Deperp fw T2a T2e].
L_all = repelem(0:2:Lmax, 2*(0:2:Lmax)+1);
N_l   = sqrt((2*L_all+1)*(4*pi));
Y     = SMI.get_even_SH(dirs,Lmax,CS_phase);
Kell  = zeros(numel(b),Lmax/2+1);
for il = 0:2:Lmax
    Kell(:,il/2+1) = SMI.RotInv_Kell_wFW_b_beta_TE_numerical(il,b,beta,TE,kernel_vec(:)',D_FW);
end
Kmat = Kell(:, repelem(1:(Lmax/2+1), 2*(0:2:Lmax)+1));
S    = (Kmat.*(Y.*N_l))*[1; plm(:)];
end

% =====================================================================
function [pk,pkdir] = fodf_peak(plm,dirs_e,Lmax,CS_phase)
% Peak amplitude of the fODF. plm holds l = 2..Lmax; the l=0 term contributes
% the fixed 1/(4*pi) isotropic floor.
Y     = SMI.get_even_SH(dirs_e,Lmax,CS_phase);
L_all = repelem(0:2:Lmax, 2*(0:2:Lmax)+1)';
[pk,i] = max(Y*([1; plm(:)].*sqrt((2*L_all+1)/(4*pi))));
pkdir  = dirs_e(i,:);
end
