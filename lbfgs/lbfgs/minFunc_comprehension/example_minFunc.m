% Runs various limited-memory solvers on 2D rosenbrock function for 25
% function evaluations
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
%  Display - Level of display [ off | final | (iter) | full | excessive ]
options.maxFunEvals = maxFunEvals;
%   MaxFunEvals - Maximum number of function evaluations allowed (1000)

% Default: L-BFGS (default)
options.Method = 'lbfgs';
%   Method - [ sd | csd | bb | cg | scg | pcg | {lbfgs} | newton0 | pnewton0 |
%       qnewton | mnewton | newton | tensor ]
x = minFunc(@k_update,x0',options);
% function [f, df, ddf] = rosenbrock(x); ====>minFuncExamples文件夹中有三维的输出

% [x,f,exitflag,output] = minFunc(funObj,x0,options,varargin)
% input:
%   funObj - is a function handle
%   x0 - is a starting vector;
%   options - is a struct containing parameters (defaults are used for non-existent or blank fields)
%   varargin{:} - all other arguments are passed as additional arguments to funObj
%
% Outputs:
%   x is the minimum value found
%   f is the function value at the minimum found
%   exitflag returns an exit condition
%   output returns a structure with other information


fprintf('x1 = %.4f, x2 = %.4f (minFunc with limited-memory BFGS - default)\n',x(1),x(2));

fprintf('---------------------------------------\n');

x = minFunc_lbfgs(@rosenbrock,[0 0]',options);

fprintf('x1 = %.4f, x2 = %.4f (minFunc with limited-memory BFGS - translate)\n',x(1),x(2));

fprintf('---------------------------------------\n');



