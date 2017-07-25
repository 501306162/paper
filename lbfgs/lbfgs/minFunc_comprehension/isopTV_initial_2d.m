%% Preperation
cd D:\GitHub\paper\lbfgs\lbfgs\minFunc_comprehension
addpath(genpath(pwd))
% mexAll

clc
clear all
close all

%% Solve problem in 2d
rho_0=5.0;
alpha=1.0;
lamda=0.1;
f_m=im2double(rgb2gray(imread('fang.png')));
f_f=im2double(rgb2gray(imread('yuan.png')));

figure,
subplot(1,2,1), imshow(f_m); title('ԭ����ͼ');
subplot(1,2,2), imshow(f_f); title('ԭ�̶�ͼ');

 % mapped to [0 1] interval
f_mw=(f_m-min(min(f_m)))./(max(max(f_m))-min(min(f_m)));
f_fw=(f_f-min(min(f_f)))./(max(max(f_f))-min(min(f_f)));

diary diary.txt
[k] = isopTV_2d(f_mw,f_fw,lamda, rho_0, alpha);
diary off
%% Reporting                                                                                                 
    [m,n]=size(f_m);
    Spacing=[2 2];

    k_m=length(0:Spacing(1):m-1);  
    k_n=length(0:Spacing(2):n-1);

    % Construct the k_grid
    % ����λ�ƾ���
    k_grid=reshape(k,k_m,k_n,2); % ��һƬΪTr,�ڶ�ƬΪTc
    k_grid=double(k_grid);

    % ������λ�ƻ�õ����ص�λ�ƾ���
    dk=bspline_transform(k_grid,Spacing,[m,n]);

    fprintf('---------------------------------------\n');   
    % �ƶ��Ժ�ͼ�� ...
    f_m=movepixels_2d(f_m,dk(:,:,1),dk(:,:,2));
    figure,
    subplot(1,3,1), imshow(f_m); title('����ͼ');
    subplot(1,3,2), imshow(f_f); title('�̶�ͼ');
    subplot(1,3,3), imshow(f_m-f_f); title('��ͼ');