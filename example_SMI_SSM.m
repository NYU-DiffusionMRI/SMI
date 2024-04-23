%% Example code for plotting sensitivity-specificity matrix (SSM) for simulations
% -----------------------------------------------------------------------------------
% please cite:  Liao, Y., Coelho, S., Chen, J., Ades-Aron, B., Pang, M., Stepanov, V., 
%               Osorio, R., Shepherd, T., Lui, Y.W., Novikov, D.S. and Fieremans, E., 2024. 
%               Mapping tissue microstructure of brain white matter in vivo in health and 
%               disease using diffusion MRI. Imaging Neuroscience.
%------------------------------------------------------------------------------------

clear;clc;close all

% placeholder for ground truth and estimates of SM parameters
N_samples = 1e3;
N_params = 6; % [p2, f, Da, Depa, Depe, ffw]
gt = rand(N_samples, N_params); 
x_smi = rand(N_samples, N_params);

% normalized by the mean
gt_norm = gt ./ mean(gt);
x_smi_norm = x_smi ./ mean(x_smi);

% calculate SSM (remove linear regression coefficients for constants)
tmp = [ones(size(gt_norm,1),1), gt_norm] \ x_smi_norm;
ssm = tmp(2:end, :);

% plot SSM
make_it_tight = true; subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.2 0.15], [0.2 0.2]);
h = figure; axis off; subplot(1,1,1);
set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
imagesc(ssm); axis square

% set color palette
c_8 = [0.192156862745098,0.211764705882353,0.584313725490196;
       0.270588235294118,0.458823529411765,0.705882352941177;
       0.670588235294118,0.850980392156863,0.913725490196078;
       1,1,0.900000000000000;
       1,1,0.900000000000000;
       0.992156862745098,0.682352941176471,0.380392156862745;
       0.843137254901961,0.188235294117647,0.152941176470588;
       0.647058823529412,0,0.149019607843137];
caxis([-0.96 0.96]);colormap(c_8);
xticks({}); yticks({});
h_colorbar = colorbar('XTickLabel',{'-0.85','-0.6','-0.35','0','0.35','0.6','0.85'},'XTick',[-0.85 -0.6 -0.35 0 0.35 0.6 0.85]);
h_colorbar.LineWidth = 1; h_colorbar.FontSize = 15; h_colorbar.Position = [0.68 0.2 0.01 0.648];

for i = 1:size(ssm,2)
    for j = 1:size(ssm,1)
        if ssm(j,i) > 0, xpos = i - 0.28; else, xpos = i - 0.32; end
        if i==j, weight = 'bold'; else, weight = 'normal'; end
            text(xpos,j-0.02,num2str(ssm(j,i),'%.2f'),'FontSize',26,'FontWeight',weight);
    end
end

for i = 0.5:size(ssm,2)+0.5
    diagonal_bold(i-0.5,i-0.5)
    for j = 0.5:size(ssm,1)+0.5
        line([j,j],[0.5,size(ssm,2)+0.5],'LineWidth',2,'Color','k')
        line([0.5,size(ssm,1)+0.5],[i,i],'LineWidth',2,'Color','k')
    end
end

para_name = {'$\hat{p}_2$','$\hat{f}$','$\hat{D}_{a}$','$\hat{D}_{e}^{\parallel}$',...
    '$\hat{D}_{e}^{\bot}$','$\hat{f}_{\mathrm{w}}$'};
xpos_para = [0.13, 0.1, 0.2, 0.2, 0.25, 0.17];
ypos_para = [0.1, 0.1, 0.05, 0.05, 0.07, 0.1];
for i = 1:size(ssm,2)
    text(i-xpos_para(i),0.3-ypos_para(i),para_name{i},'Interpreter','latex','FontSize',36)
end
para_name = {'${p}_2$','$f$','${D}_{a}$','${D}_{e}^{\parallel}$','${D}_{e}^{\bot}$','$f_{\mathrm{w}}$'};
xpos_para = [0.07, 0.04, 0.15, 0.18, 0.22, 0.13];
ypos_para = [0.015, 0.02, 0.01, 0.01, 0.02, 0.015];
for i = 1:size(ssm,1)
    text(0.1-xpos_para(i),i-ypos_para(i),para_name{i},'Interpreter','latex','FontSize',36)
end

%% supporting functions

function diagonal_bold(i,j)
    line([i-0.5 i+0.5],[j-0.5 j-0.5],'LineWidth',5,'Color','k')
    line([i-0.5 i+0.5],[j+0.5 j+0.5],'LineWidth',5,'Color','k')
    line([i-0.5 i-0.5],[j-0.5 j+0.5],'LineWidth',5,'Color','k')
    line([i+0.5 i+0.5],[j-0.5 j+0.5],'LineWidth',5,'Color','k')
end

function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
%
% Functional purpose: A wrapper function for Matlab function subplot. Adds the ability to define the gap between
% neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the gap between
% subplots can reach 40% of figure area, which is pretty lavish.  
%
% Input arguments (defaults exist):
%   gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%            relatively large axis. 
%   marg_h  margins in height in normalized units (0...1)
%            or [lower uppper] for different lower and upper margins 
%   marg_w  margins in width in normalized units (0...1)
%            or [left right] for different left and right margins 
%
% Output arguments: same as subplot- none, or axes handle according to function call.
%
% Issues & Comments: Note that if additional elements are used in order to be passed to subplot, gap parameter must
%       be defined. For default gap value use empty element- [].      
%
% Usage example: h=subtightplot((2,3,1:2,[0.5,0.2])

if (nargin<4) || isempty(gap),    gap=0.01;  end
if (nargin<5) || isempty(marg_h),  marg_h=0.05;  end
if (nargin<5) || isempty(marg_w),  marg_w=marg_h;  end
if isscalar(gap),   gap(2)=gap;  end
if isscalar(marg_h),  marg_h(2)=marg_h;  end
if isscalar(marg_w),  marg_w(2)=marg_w;  end
gap_vert   = gap(1);
gap_horz   = gap(2);
marg_lower = marg_h(1);
marg_upper = marg_h(2);
marg_left  = marg_w(1);
marg_right = marg_w(2);

%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);  

% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   

% single subplot dimensions:
height=(1-(marg_lower+marg_upper)-(m-1)*gap_vert)/m;
width =(1-(marg_left+marg_right)-(n-1)*gap_horz)/n;

% merged subplot dimensions:
merged_height=subplot_rows*( height+gap_vert )- gap_vert;
merged_width= subplot_cols*( width +gap_horz )- gap_horz;

% merged subplot position:
merged_bottom=(m-max(subplot_row))*(height+gap_vert) +marg_lower;
merged_left=(min(subplot_col)-1)*(width+gap_horz) +marg_left;
pos_vec=[merged_left merged_bottom merged_width merged_height];

% h_subplot=subplot(m,n,p,varargin{:},'Position',pos_vec);
% Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
h=subplot('Position',pos_vec,varargin{:});

if (nargout < 1),  clear h;  end

end
