%  Basic demon registration code. (To easy understand the algorithm)
try
    functionname='basic_demon_example.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions_nonrigid'])

catch me
    disp(me.message);
end

% Clean
clc; clear all; close all;

%% Read two images  
I1=im2double(imread('images/brain1.png'));  
I2=im2double(imread('images/brain2.png'));

% Set static and moving image
M=I1;S=I2; 

% Alpha (noise) constant
alpha=0.5;

%权值
d=1.5;

% The transformation fields
Tx=zeros(size(M)); Ty=zeros(size(M));
Tx1=zeros(size(M)); Ty1=zeros(size(M));
Tx2=zeros(size(M)); Ty2=zeros(size(M));
cx=zeros(size(M)); cy=zeros(size(M));
w=zeros(size(M));

% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],10);

[Sy,Sx] = gradient(S);

for itt=1:200
     % Difference image between moving and static image
        Idiff=M-S;
        
%         %1.0 Default demon force, (Thirion 1998)    
%         Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
%         Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         Extended demon force. With forces from the gradients from both 
% %         moving as static image. (Cachier 1999, He Wang 2005)
%         [My,Mx] = gradient(M);
%         Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
%         Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %% 2.2 Extended demon force.  (2014)   本文方法：引入平衡系数
            gx=d.*Sx;
            gy=d.*Sy;
            [My,Mx] = gradient(M);
            kx=d.*Mx;
            ky=d.*My;             
            Ux = -Idiff.*  ((Sx./((gx.^2+gy.^2)+alpha^2*Idiff.^2))+(Mx./((kx.^2+ky.^2)+alpha^2*Idiff.^2)));
            Uy = -Idiff.*  ((Sy./((gx.^2+gy.^2)+alpha^2*Idiff.^2))+(My./((kx.^2+ky.^2)+alpha^2*Idiff.^2)));
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % When divided by zero
        Ux(isnan(Ux))=0; 
        Uy(isnan(Uy))=0;

        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);

        % Add the new transformation field to the total transformation field.
        Tx=Tx+Uxs;
        Ty=Ty+Uys;
        M=movepixels(I1,Tx,Ty);
        mse1(itt)=myMSE(I2,M);
end
% figure,showgrid(Tx,Ty,4); title('Transformation');

figure,
subplot(1,3,1), imshow(I1,[]); title('moving image');
subplot(1,3,2), imshow(I2,[]); title('static image');
subplot(1,3,3), imshow(M,[]); title('Registered moving image');

% Show the difference 
figure,
subplot(1,2,1), imshow(I2-I1,[]); title('the differnence of before');
subplot(1,2,2), imshow(I2-M,[]); title('the differnence of after');

%% 增加评估MSD并输出
mse=myMSE(S,M);
tss=now();
names=strcat(num2str(hour(tss)),'_',num2str(minute(tss)));
namex=strcat('Yp_',names,'.mat');
fprintf(strcat(namex,'\n'))
save(namex,'mse1','mse','M');


	