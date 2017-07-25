function [k,f] = LBFGS_k(k0,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff,M_lbgfs )

maxFunEvals = M_lbgfs;

fprintf('Result after %d evaluations of limited-memory solvers on 2D rosenbrock:\n',maxFunEvals);

% fprintf('---------------------------------------\n');

options = [];% Ԫ��Ϊ��
options.display = 'none';  
options.maxFunEvals = maxFunEvals;
options.Method = 'lbfgs';
% ����   Ҫ����Ed_derivation2_2d(k,ffo,fmo,rho,Dissimilarity,spacing,ddk,ZU_diff)
if numel(size(ffo))<3
    [k,f] = minFunc(@Ed_2d,k0,options,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff);
else 
    [k,f] = minFunc(@Ed_3d,k0,options,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff);     
end
    [m,n]=size(fmo);   
    [kx,ky]=ndgrid(0:Spacing(1):m-1,0:Spacing(2):n-1); % ��������

    k_m=size(kx,1);  
    k_n=size(kx,2);
    num_control=k_m*k_n;     
    % Construct the k_grid

    % ����λ�ƾ���
    k_grid=reshape(k,k_m,k_n,2); % ��һƬΪTr,�ڶ�ƬΪTc
    k_grid=double(k_grid);

    % ������λ�ƻ�õ����ص�λ�ƾ���
    dk=bspline_transform(k_grid,Spacing,[m,n]);

    % �ƶ��Ժ�ͼ�� ...
    f_m=movepixels_2d(fmo,dk(:,:,1),dk(:,:,2));
    figure,
    subplot(1,2,1), imshow(f_m); title('����ͼ');
    subplot(1,2,2), imshow(ffo); title('�̶�ͼ');

end

