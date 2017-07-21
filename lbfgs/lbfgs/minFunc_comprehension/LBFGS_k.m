function [k,f] = LBFGS_k(k0,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff,M_lbgfs )

maxFunEvals = M_lbgfs;

fprintf('Result after %d evaluations of limited-memory solvers on 2D rosenbrock:\n',maxFunEvals);

% fprintf('---------------------------------------\n');
% fprintf('starting point :\n');
% k0' % k0 nl行1列
% fprintf('---------------------------------------\n');

options = [];% 元素为空
options.display = 'none';  
options.maxFunEvals = maxFunEvals;
options.Method = 'lbfgs';
% 邻域   要计算Ed_derivation2_2d(k,ffo,fmo,rho,Dissimilarity,spacing,ddk,ZU_diff)
if length(size(ffo))<3
    [k,f] = minFunc(@Ed_2d,k0,options,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff);
else 
    [k,f] = minFunc(@Ed_3d,k0,options,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff);     
end
    
% fprintf('minFunc with limited-memory BFGS \n');
% k'
% fprintf('---------------------------------------\n');
end

