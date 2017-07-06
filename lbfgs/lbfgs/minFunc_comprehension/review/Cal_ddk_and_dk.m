function [dk,ddk]=Cal_ddk_and_dk(kr,kc,m,n)

% b-spline grid spacing in x and y direction
Spacing=[4 4];

% Construct the k_grid
O_m=length(0:Spacing(1):m);
O_n=length(0:Spacing(2):n);

O_trans=zeros(O_m,O_n,2);
O_trans(:,:,1)=kr;
O_trans(:,:,2)=kc;
O_trans=double(O_trans);
 
% Transform the image with the B-spline grid
dk=bspline_transform(O_trans,Spacing,[m,n]);
ddk = Cal_ddk( m,n,Spacing );
end