
function [Tx,Ty]=demos(imobserve,imref)
% Read two images
I2=imref;
I1=imobserve;

% Set static and moving image
S=I2; M=I1;

% Alpha (noise) constant
alpha=1;
d=1.5;

% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],10);

% The transformation fields
Tx=zeros(size(M)); Ty=zeros(size(M));

[Sy,Sx] = gradient(S);
for itt=1:200
	    % Difference image between moving and static image
        Idiff=M-S;
        
        gx=d.*Sx;
        gy=d.*Sy;
        [My,Mx] = gradient(M);
        kx=d.*Mx;
        ky=d.*My;             
        Ux = -Idiff.*  ((Sx./((gx.^2+gy.^2)+alpha^2*Idiff.^2))+(Mx./((kx.^2+ky.^2)+alpha^2*Idiff.^2)));
        Uy = -Idiff.*  ((Sy./((gx.^2+gy.^2)+alpha^2*Idiff.^2))+(My./((kx.^2+ky.^2)+alpha^2*Idiff.^2)));
 
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;

        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);

        % Add the new transformation field to the total transformation field.
        Tx=Tx+Uxs;
        Ty=Ty+Uys;
        M=movepixels(I1,Tx,Ty); 
end
end