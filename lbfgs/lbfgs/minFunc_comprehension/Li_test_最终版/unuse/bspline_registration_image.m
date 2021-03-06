function [I_tran,Trans]=bspline_registration_image(O_grid,sizes,Spacing,I1)
% Function Registration_image. This function will create
% an small registration error image after b-spline non-rigid registration
% of two images / volumes.
%
% Ierror=bspline_registration_image(grid,sizes,I1,I2,type);
%
% inputs,
%   grid: b-spline control grid created by make_init_grid.m, (and reshaped
%         to one long vector.)
%   sizes: sizes need tot reshape the grid to matrix format.
%   Spacing: Spacing of the b-spline grid knot vector
%   I1 and I2: The input images, of which I1 is transformed.
%   type: Type of registration error (see registration_error.m)
%
% outputs,
%   Ierror: The error image, resized for usage by lsqnonlin optimizer.
%
% example,
%   I1=im2double(imread('lenag1.png')); 
%   I2=im2double(imread('lenag2.png'));
%   O_trans=make_init_grid([8 8],size(I1)); sizes=size(O_trans);
%   O_trans = lsqnonlin(@(x)bspline_registration_image(x,sizes,I1,I2,'d'),O_trans(:),[],[],optimset('Display','iter','MaxIter',10));
%   Icor=bspline_transform(reshape(O_trans,sizes),I1); 
%   figure, imshow(I1), figure, imshow(I2), figure, imshow(Icor)
%
%
% This function is written by D.Kroon University of Twente (July 2008)


O_grid=reshape(O_grid,sizes);

% Transform the image with the B-spline grid
[I_tran,Trans]=bspline_transform(O_grid,I1,Spacing);

end

