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
options = [];% ÔªËØÎª¿Õ
options.display = 'none';  
options.maxFunEvals = maxFunEvals;
options.Method = 'lbfgs';

x0=
x = minFunc(@Ed_derivation,x0',options);
fprintf('x1 = %.4f, x2 = %.4f (minFunc with limited-memory BFGS - default)\n',x(1),x(2));

fprintf('---------------------------------------\n');
end

