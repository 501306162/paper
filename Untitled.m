% data set:
% f_m & f_f : images
% lamda     : regularization    0.1
% w         : LCC kernel bandwith   2.5
% Mpyr      : # of image pyramid levels  log2(m/8)
% Mref      : # of grid refinements
% Miter     : max.# of ADMM iterations
% M_lbgfs   : max.# of LBFGS iterations   5
% kesai_tol :optimization argument tolerance(10^-2)
% 
% ADMM internal parameters : ro=5;mu=10;tao=1.3;
% =========================================================

f_m=randi(10,10,10);
f_f=randi(10,10,10);
[m,n]=size(f_m);
lamda=0.1;
w=2.5;
Mpyr=ceil(log2(m/8));
Mref=5;
Miter=5;
M_lbgfs=5;
