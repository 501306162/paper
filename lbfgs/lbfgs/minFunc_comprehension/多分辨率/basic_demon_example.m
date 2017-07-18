%  Basic demon registration code. (To easy understand the algorithm)

% Clean
clc; clear all; close all;

% Compile the mex files
compile_c_files

functionname='basic_demon_example.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir '/functions_nonrigid'])
addpath([functiondir '/functions_affine'])
addpath([functiondir '/functions'])
addpath([functiondir '/lianxi'])
% Read two images
% 
I1=im2double(imread('images/square.bmp'));  %imread('images/lenag1.png')
I2=im2double(imread('images/circle.bmp')); %imread('images/lenag2.png')
% % % 
% I1=im2double(imread('images/lenag1.png'));  %imread('images/lenag1.png')
% I2=im2double(imread('images/lenag2.png')); %imread('images/lenag2.png')
% 
% I1=im2double(imread('images/checkboard4.png'));  %imread('images/lenag1.png')
% I2=im2double(imread('images/checkboard2.png')); %imread('images/lenag2.png')

% Set static and moving image
S=I2; M=I1;

% Alpha (noise) constant
alpha=0.5:0.5:4;
% beta=0.1:0.5:3;zeros(size(M))
% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],10);

Txx=0; Tyy=0;l=1;

for L=3:-1:1
% for aa=1%:length(alpha)
%     for bb=1:length(beta)
%     M=I1;

    M0(:,:,l)=down_sample(I1,L);
    S0(:,:,l)=down_sample(I2,L);
    M=M0(:,:,l);
    S=S0(:,:,l);
% The transformation fields
% Tx=zeros(size(M)); Ty=zeros(size(M));
    t=zeros(size(M(:,:,l)));

    [Sy,Sx] = gradient(S(:,:,l));

    aa=1;
    c=['b','g','r','m', 'y','k','c',':'];
    diary on
% % figure,hold on
Tx=Txx(:,:,l);
Ty=Tyy(:,:,l);

 
   
for itt=1:200
	    % Difference image between moving and static image
        Idiff=M-S;
%         xx(kk,itt)=length(find(sign(Idiff.*Sx)<0));
        [My,Mx] = gradient(M);
        for kk=1:(size(My,1)-2)
            for jj=1:(size(My,2)-2)
               t(kk,jj)=TAR(Mx(kk,jj),My(kk,jj), Mx(kk+1,jj+1),My(kk+1,jj+1)  ,Mx(kk+2,jj+2),My(kk+2,jj+2));
            end
        end
       [ty,tx] = gradient(t);
%         k=Ant_K(M,I1,I2);
% %         Default demon force, (Thirion 1998)
%         Ux = -(Idiff.*Mx)./((Mx.^2+My.^2)+(alpha.*Idiff).^2); %.*((beta(bb))^2
%         Uy = -(Idiff.*My)./((Mx.^2+My.^2)+(alpha.*Idiff).^2);

        Ux = -Idiff.*  ((tx./((tx.^2+ty.^2)+(alpha(aa).^2).*Idiff.^2)) +(Mx./((Mx.^2+My.^2)+((alpha(aa).^2).*Idiff.^2)) ) );
        Uy = -Idiff.*  ((ty./((tx.^2+ty.^2)+(alpha(aa).^2).*Idiff.^2)) +(My./((Mx.^2+My.^2)+((alpha(aa).^2).*Idiff.^2)) ) );
%  
        % Extended demon force. With forces from the gradients from both
% % % %         moving as static image. (Cachier 1999, He Wang 2005)
% %         [My,Mx] = gradient(M);
%         Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+(alpha(aa).*Idiff).^2))+(Mx./((Mx.^2+My.^2)+(alpha(aa).*Idiff).^2)));
%         Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+(alpha(aa).*Idiff).^2))+(My./((Mx.^2+My.^2)+(alpha(aa).*Idiff).^2)));
 
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;  %判断是否是非数值参数

        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);

        % Add the new transformation field to the total transformation field.
        Tx=Tx+Uxs;
        Ty=Ty+Uys;
        M=movepixels(M0(:,:,l),Tx,Ty); %I1 

%         mii(itt)=CalculateMI(M,I2);
%         psni(itt)=psnr(M,I2);
%         msi(itt)=mse(M,I2);
%         ssii(itt)=ssim(M,I2);

end
figure,
subplot(1,3,1), imshow(I1,[]); title('image 1');
subplot(1,3,2), imshow(I2,[]); title('image 2');
subplot(1,3,3), imshow(M,[]); title('Registered image 1');
Txx(:,:,l)=Tx(:,:);
Tyy(:,:,l)=Ty(:,:);

Txx(:,:,l+1)=up_sample(Tx,L);
Tyy(:,:,l+1)=up_sample(Ty,L);

l=l+1;
% figure(1),hold on
% plot(1:itt,mii,c(aa));
% figure(2),hold on
% plot(1:itt,psni,c(aa));
% figure(3),hold on
% plot(1:itt,msi,c(aa));
% figure(4),hold on
% plot(1:itt,ssii,c(aa));
% fprintf('==============================================\n');
% ms(aa)=mse(M,I2);%,bb
% fprintf('the mse is %8.5f\n',ms(aa));
% mi(aa)=CalculateMI(M,I2);
% fprintf('the mi is %8.5f\n',mi(aa));
% psn(aa)=psnr(M,I2);
% fprintf('the psnr is %8.5f\n',psn(aa));
% ssi(aa)=ssim(M,I2);
% fprintf('the ssim is %8.5f\n',ssi(aa));
end 

M=imresize(M,'bicubic');
figure, 
subplot(1,3,1), imshow(I1,[]); title('image 1');
subplot(1,3,2), imshow(I2,[]); title('image 2');
subplot(1,3,3), imshow(M,[]); title('Registered image 1');

      
% end
% 
% c=rand(1,3);
% q=1:50:200;
% plot(q,mi(q),'color',c);
%    MI(kk)=sum(mi)/200; 
% 
% end
% hold off
% figure,
% plot(k,MI);title(' k & MI ');

% [alpha,beta]=meshgrid(alpha,beta);
%  mesh(alpha,beta,ms)

diary off