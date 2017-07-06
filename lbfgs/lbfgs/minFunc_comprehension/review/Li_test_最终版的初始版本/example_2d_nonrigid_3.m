
% Add all function paths
addpaths

% Read two greyscale images of Lena
I1=double(deformation_image); 

% b-spline grid spacing in x and y direction
Spacing=[4 4];
 
O_m=length(0:Spacing(1):m);
O_n=length(0:Spacing(2):n);

O_trans=ones(O_m,O_n,2);
O_trans(:,:,1)=kr;
O_trans(:,:,2)=kc;
O_trans=double(O_trans);
 
% Transform the image with the B-spline grid
[I_tran,dkx,dky]=bspline_transform(O_trans,I1,Spacing);

% Reshape O_trans from a vector to a matrix
k_grid=reshape(k_grid,sizes);

% Transform the input image with the found optimal grid.
Icor=bspline_transform(k_grid,I1,Spacing); 

% Make a (transformed) grid image
Igrid=make_grid_image(Spacing,size(I1));
Igrid=bspline_transform(k_grid,Igrid,Spacing); 

% Show the registration results
figure,
subplot(2,2,1), imshow(I1); title('input image 1');
subplot(2,2,2), imshow(I2); title('input image 2');
subplot(2,2,3), imshow(Icor); title('transformed image 1');
subplot(2,2,4), imshow(Igrid); title('grid');
