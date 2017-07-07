% Group lasso example with random data
cd E:\GitHub\paper\lbfgs\lbfgs\minFunc_comprehension
addpath(genpath(pwd))
mexAll
%% Generate problem data
clc
clear all
close all
randn('seed', 0);
rand('seed',0);


m = 1500;       % amount of data
K = 200;        % number of blocks
partition = randi(50, [K 1]);

n = sum(partition); % number of features
p = 100/n;          % sparsity density  

% generate block sparse solution vector
x = zeros(n,1);
start_ind = 1;
cum_part = cumsum(partition);
for i = 1:K,
    x(start_ind:cum_part(i)) = 0;
    if( rand() < p)
        % fill nonzeros
        x(start_ind:cum_part(i)) = randn(partition(i),1);
    end
    start_ind = cum_part(i)+1;
end

% % lambda max
% start_ind = 1;
% for i = 1:K,
%     sel = start_ind:cum_part(i);
%     lambdas(i) = norm(A(:,sel)'*b);
%     start_ind = cum_part(i) + 1;
% end
% lambda_max = max(lambdas);
% 
% % regularization parameter
% lambda = 0.1*lambda_max; 
lamda=0.1;
%% Solve problem
rho_0=5.0;
alpha=1.0;
% f_m=randi(5,10,10);
% f_f=f_m;
% 
% f_f(10,10)=f_m(10,1);
% f_f(10,1)=f_m(10,10);
% 
% figure,
% subplot(1,2,1), imshow(f_m,[]); title('原浮动图');
% subplot(1,2,2), imshow(f_f,[]); title('原固定图');
% f_m=im2double(imread('fang.png'));
% f_f=im2double(imread('yuan.png'));
% 
% figure,
% subplot(1,2,1), imshow(f_m); title('原浮动图');
% subplot(1,2,2), imshow(f_f); title('原固定图');
f_m=[0,0,0,0,0, 0,0,0,0,0;
    0,0,0,0,0, 0,0,0,1,0;
    0,0,1,1,1, 1,1,0,0,0;
    0,1,1,0,1, 1,1,1,0,0;
    0,0,1,1,1, 1,1,1,0,0;
    
    0,0,1,1,1, 1,1,1,0,0;
    0,0,0,1,1, 1,1,1,0,0
    0,0,1,0,1, 1,1,1,0,0;
    0,1,0,0,0, 0,0,0,0,0;
    0,0,0,1,0, 0,0,0,0,0;];

f_f=[0,0,0,0,0, 0,0,0,0,0;
    0,0,0,0,0, 0,0,0,0,0;
    0,0,1,1,1, 1,1,1,0,0;
    0,0,1,1,1, 1,1,1,0,0;
    0,0,1,1,1, 1,1,1,0,0;
    
    0,0,1,1,1, 1,1,1,0,0;
    0,0,1,1,1, 1,1,1,0,0
    0,0,1,1,1, 1,1,1,0,0;
    0,0,0,0,0, 0,0,0,0,0;
    0,0,0,0,0, 0,0,0,0,0;];

figure,
subplot(1,2,1), imshow(f_m,[]); title('原浮动图');
subplot(1,2,2), imshow(f_f,[]); title('原固定图');

diary diary.txt
[x] = isopTV_new(f_m,f_f,lamda,partition, rho_0, alpha);
diary off
%% Reporting                                                                                                 


