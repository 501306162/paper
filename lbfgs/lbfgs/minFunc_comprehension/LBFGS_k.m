function [ output_args ] = LBFGS_k( input_args )

maxFunEvals = 25;

fprintf('Result after %d evaluations of limited-memory solvers on 2D rosenbrock:\n',maxFunEvals);

fprintf('---------------------------------------\n');
fprintf('x1 = %.4f, x2 = %.4f (starting point)\n',0,0);
fprintf('x1 = %.4f, x2 = %.4f (optimal solution)\n',1,1);
fprintf('---------------------------------------\n');

if exist('minimize') == 2
    % Minimize.m - conjugate gradient method
    x = minimize([0 0]', 'rosenbrock', -maxFunEvals);
    fprintf('x1 = %.4f, x2 = %.4f (minimize.m by C. Rasmussen)\n',x(1),x(2));
end
options = [];% 元素为空
options.display = 'none';  
options.maxFunEvals = maxFunEvals;
options.Method = 'lbfgs';
邻域   

        neig_zero=neighborhood_trans(ind,ffo_zero,spacing); % K1,K2,K3大小邻域...中的点对应在原图中标号 

x0=
x = minFunc(@Ed_derivation,x0',options);
fprintf('x1 = %.4f, x2 = %.4f (minFunc with limited-memory BFGS - default)\n',x(1),x(2));

fprintf('---------------------------------------\n');
end

