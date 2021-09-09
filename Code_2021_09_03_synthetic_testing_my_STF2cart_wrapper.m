%% ORDER 2 (correct in rank2)
clc,clear,close all
y=getY(2);
Y_ij_00 =y{1,1,1};
Y_ij_2m2=y{1,2,1};
Y_ij_2m1=y{1,2,2};
Y_ij_20 =y{1,2,3};
Y_ij_21 =y{1,2,4};
Y_ij_22 =y{1,2,5};

d=[1 rand(1,5)*10];
D=d(1)*Y_ij_00+d(2)*Y_ij_2m2+d(3)*Y_ij_2m1+d(4)*Y_ij_20+d(5)*Y_ij_21+d(6)*Y_ij_22;

C0=sqrt(1/(4*pi));
C2=sqrt(5/(4*pi));
d00 =1/3*(1/C0^2)*Y_ij_00(:)'*D(:);
d2m2=2/3*(1/C2^2)*Y_ij_2m2(:)'*D(:);
d2m1=2/3*(1/C2^2)*Y_ij_2m1(:)'*D(:);
d20 =2/3*(1/C2^2)*Y_ij_20(:)'*D(:);
d21 =2/3*(1/C2^2)*Y_ij_21(:)'*D(:);
d22 =2/3*(1/C2^2)*Y_ij_22(:)'*D(:);

all=[d;d00 d2m2 d2m1 d20 d21 d22];
all(1,:)-all(2,:)

d_stf=[d00 d2m2 d2m1 d20 d21 d22]';
d_cart = wrapper_STF_decomposition_backandforth(d_stf,'STF2cart',[]);
D([1 5 9 2 3 6])-d_cart'
d_stf_back = wrapper_STF_decomposition_backandforth(d_cart,'cart2STF',[]);
(d_stf-d_stf_back)'

% [trace(D) 3*d00/sqrt(4*pi)]


%% ORDER 4 (correct in rank4)
clc,clear,close all
y=getY(4);
Y_ij_00 =y{2,1,1};
Y_ij_2m2=y{2,2,1};
Y_ij_2m1=y{2,2,2};
Y_ij_20 =y{2,2,3};
Y_ij_21 =y{2,2,4};
Y_ij_22 =y{2,2,5};
Y_ij_4m4=y{2,3,1};
Y_ij_4m3=y{2,3,2};
Y_ij_4m2=y{2,3,3};
Y_ij_4m1=y{2,3,4};
Y_ij_40 =y{2,3,5};
Y_ij_41 =y{2,3,6};
Y_ij_42 =y{2,3,7};
Y_ij_43 =y{2,3,8};
Y_ij_44 =y{2,3,9};

d=[rand(1,15)*10];
D=d(1)*Y_ij_00+d(2)*Y_ij_2m2+d(3)*Y_ij_2m1+d(4)*Y_ij_20+d(5)*Y_ij_21+d(6)*Y_ij_22+d(7)*Y_ij_4m4+d(8)*Y_ij_4m3+d(9)*Y_ij_4m2+d(10)*Y_ij_4m1+d(11)*Y_ij_40+d(12)*Y_ij_41+d(13)*Y_ij_42+d(14)*Y_ij_43+d(15)*Y_ij_44;

C0=sqrt(1/(4*pi));
C2=sqrt(5/(4*pi));
C4=sqrt(9/(4*pi));
d00 =1/5*(1/C0^2)*Y_ij_00(:)'*D(:);
d2m2=4/7*(1/C2^2)*Y_ij_2m2(:)'*D(:);
d2m1=4/7*(1/C2^2)*Y_ij_2m1(:)'*D(:);
d20 =4/7*(1/C2^2)*Y_ij_20(:)'*D(:);
d21 =4/7*(1/C2^2)*Y_ij_21(:)'*D(:);
d22 =4/7*(1/C2^2)*Y_ij_22(:)'*D(:);
d4m4=8/35*(1/C4^2)*Y_ij_4m4(:)'*D(:);
d4m3=8/35*(1/C4^2)*Y_ij_4m3(:)'*D(:);
d4m2=8/35*(1/C4^2)*Y_ij_4m2(:)'*D(:);
d4m1=8/35*(1/C4^2)*Y_ij_4m1(:)'*D(:);
d40 =8/35*(1/C4^2)*Y_ij_40(:)'*D(:);
d41 =8/35*(1/C4^2)*Y_ij_41(:)'*D(:);
d42 =8/35*(1/C4^2)*Y_ij_42(:)'*D(:);
d43 =8/35*(1/C4^2)*Y_ij_43(:)'*D(:);
d44 =8/35*(1/C4^2)*Y_ij_44(:)'*D(:);

all=[d;d00 d2m2 d2m1 d20 d21 d22 d4m4 d4m3 d4m2 d4m1 d40 d41 d42 d43 d44];
all(1,:)-all(2,:)

d_stf=[d00 d2m2 d2m1 d20 d21 d22 d4m4 d4m3 d4m2 d4m1 d40 d41 d42 d43 d44]';
d_cart = wrapper_STF_decomposition_backandforth(d_stf,'STF2cart',[]);

idx=[ 1,1,1,1; 2,2,2,2; 3,3,3,3; 1,1,2,2; 1,1,3,3; 1,1,1,2; 1,1,1,3; 1,1,2,3; 2,2,3,3; 2,2,1,2; 2,2,1,3; 2,2,2,3; 3,3,1,2; 3,3,1,3; 3,3,2,3];
ind = sub2ind([3 3 3 3],idx(:,1),idx(:,2),idx(:,3),idx(:,4));
(D(ind)-d_cart)'

d_stf_back = wrapper_STF_decomposition_backandforth(d_cart,'cart2STF',[]);
(d_stf-d_stf_back)'

%% CART -> STF -> CART
clc,clear,close all
a2=symm(randn(3,3)*10);
a2_cart=a2([1 5 9 2 3 6])';
a2_stf  = wrapper_STF_decomposition_backandforth(a2_cart,'cart2STF',[]);
a2_cart_back = wrapper_STF_decomposition_backandforth(a2_stf,'STF2cart',[]);
(a2_cart-a2_cart_back)'

idx=[ 1,1,1,1; 2,2,2,2; 3,3,3,3; 1,1,2,2; 1,1,3,3; 1,1,1,2; 1,1,1,3; 1,1,2,3; 2,2,3,3; 2,2,1,2; 2,2,1,3; 2,2,2,3; 3,3,1,2; 3,3,1,3; 3,3,2,3];
ind = sub2ind([3 3 3 3],idx(:,1),idx(:,2),idx(:,3),idx(:,4));
a4=symm(randn(3,3,3,3)*10);
a4_cart=a4(ind)k;
a4_stf  = wrapper_STF_decomposition_backandforth(a4_cart,'cart2STF',[]);
a4_cart_back = wrapper_STF_decomposition_backandforth(a4_stf,'STF2cart',[]);
(a4_cart-a4_cart_back)'



