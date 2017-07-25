function [k,f] = LBFGS_k(k0,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff,M_lbgfs )

maxFunEvals = M_lbgfs;

fprintf('Result after %d evaluations of limited-memory solvers on 2D rosenbrock:\n',maxFunEvals);

% fprintf('---------------------------------------\n');

options = [];% 元素为空
options.display = 'none';  
options.maxFunEvals = maxFunEvals;
options.Method = 'lbfgs';
% 邻域   要计算Ed_derivation2_2d(k,ffo,fmo,rho,Dissimilarity,spacing,ddk,ZU_diff)
if numel(size(ffo))<3
    [k,f] = minFunc(@Ed_2d,k0,options,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff);
else 
    [k,f] = minFunc(@Ed_3d,k0,options,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff);     
end
    [m,n]=size(fmo);   
    [kx,ky]=ndgrid(0:Spacing(1):m-1,0:Spacing(2):n-1); % 网格坐标

    k_m=size(kx,1);  
    k_n=size(kx,2);
    num_control=k_m*k_n;     
    % Construct the k_grid

    % 网格位移矩阵
    k_grid=reshape(k,k_m,k_n,2); % 第一片为Tr,第二片为Tc
    k_grid=double(k_grid);

    % 由网格位移获得的像素点位移矩阵
    dk=bspline_transform(k_grid,Spacing,[m,n]);

    % 移动以后图像 ...
    f_m=movepixels_2d(fmo,dk(:,:,1),dk(:,:,2));
    figure,
    subplot(1,2,1), imshow(f_m); title('浮动图');
    subplot(1,2,2), imshow(ffo); title('固定图');

end

