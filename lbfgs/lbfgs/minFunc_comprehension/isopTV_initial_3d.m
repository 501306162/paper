%% Preperation
cd D:\GitHub\paper\lbfgs\lbfgs\minFunc_comprehension
addpath(genpath(pwd))
% mexAll

clc
clear all
close all

%% Solve problem in 3d
rho_0=5.0;
alpha=1.0;
lamda=0.1;
f_m=im2double(rgb2gray(imread('fang.png')));
f_f=im2double(rgb2gray(imread('yuan.png')));

figure,
subplot(1,2,1), imshow(f_m); title('原浮动图');
subplot(1,2,2), imshow(f_f); title('原固定图');

 % mapped to [0 1] interval
f_mw=(f_m-min(min(f_m)))./(max(max(f_m))-min(min(f_m)));
f_fw=(f_f-min(min(f_f)))./(max(max(f_f))-min(min(f_f)));

diary diary.txt
[k] = isopTV_3d(f_mw,f_fw,lamda, rho_0, alpha);
diary off
%% Reporting                                                                                                 
    [m,n,q]=size(f_m);
    Spacing=[2 2 2];

    k_m=length(0:Spacing(1):m-1);  
    k_n=length(0:Spacing(2):n-1);
    k_q=length(0:Spacing(3):q-1);
    
    % Construct the k_grid
    % 网格位移矩阵
    k_grid=reshape(k,k_m,k_n,k_q,3); % 第一片为Tr,第二片为Tc...
    k_grid=double(k_grid);

    % 由网格位移获得的像素点位移矩阵
    dk=bspline_transform(k_grid,Spacing,[m,n,q]);

    fprintf('---------------------------------------\n');   
    % 移动以后图像 ...
    f_m=movepixels_3d(f_m,dk(:,:,:,1),dk(:,:,:,2),dk(:,:,:,3));
    figure,
    subplot(1,2,1), imshow(f_m); title('浮动图');
    subplot(1,2,2), imshow(f_f); title('固定图');
