function [k] = LBFGS_k(k0,ffo,fmo,rho,Dissimilarity,spacing,ZU_diff )

maxFunEvals = 25;

fprintf('Result after %d evaluations of limited-memory solvers on 2D rosenbrock:\n',maxFunEvals);

fprintf('---------------------------------------\n');
fprintf('starting point :\n');
k0' % k0 nl行1列
fprintf('---------------------------------------\n');

options = [];% 元素为空
options.display = 'none';  
options.maxFunEvals = maxFunEvals;
options.Method = 'lbfgs';
% 邻域   要计算Ed_derivation2_2d(k,ffo,fmo,rho,Dissimilarity,spacing,ddk,ZU_diff)

k = minFunc(@Ed_2d,k0,options,ffo,fmo,rho,Dissimilarity,spacing,ZU_diff);

fprintf('minFunc with limited-memory BFGS \n');
k'
fprintf('---------------------------------------\n');
end

