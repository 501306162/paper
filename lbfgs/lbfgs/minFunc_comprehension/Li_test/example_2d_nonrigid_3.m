
% Add all function paths
addpaths

% Read two greyscale images of Lena
I1=im2double(f_mj); 
I2=im2double(f_fj);

% b-spline grid spacing in x and y direction
Spacing=[4 4];

% Make the Initial b-spline registration grid
[k_grid]=make_init_grid(Spacing,size(I1));

% Convert all values tot type double
I1=double(I1); I2=double(I2); k_grid=double(k_grid); 

% Smooth both images for faster registration
I1s=imfilter(I1,fspecial('gaussian',[20 20],5));
I2s=imfilter(I2,fspecial('gaussian',[20 20],5));

% Optimizer parameters
optim=optimset('Display','iter','MaxIter',40);

% Reshape O_trans from a matrix to a vector.
sizes=size(k_grid); k_grid=k_grid(:);

% Start the b-spline nonrigid registration optimizer
k_grid = lsqnonlin(@(x)bspline_registration_image(x,sizes,Spacing,I1s,I2s,type),k_grid,[],[],optim);

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
