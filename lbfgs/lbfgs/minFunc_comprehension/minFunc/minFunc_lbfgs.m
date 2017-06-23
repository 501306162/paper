function [x,f,exitflag] = minFunc_lbfgs(funObj,x0,options,varargin)
% [x,f,exitflag,output] = minFunc(funObj,x0,options,varargin)
%
% Unconstrained optimizer using a line search strategy
%
% Uses an interface very similar to fminunc
%

% It computes descent directions using one of ('Method'):
%   - 'lbfgs': Quasi-Newton with Limited-Memory BFGS Updating
%  (default: uses a predetermined nunber of previous steps to form a low-rank Hessian approximation)

% Several line search strategies are available for finding a step length satisfying
%   the termination criteria ('LS_type')
%   - 0 : A backtracking line-search based on the Armijo condition (default for 'bb')
%   - 1 : A bracekting line-search based on the strong Wolfe conditions (default for all other methods)
%   - 2 : The line-search from the Matlab Optimization Toolbox (requires Matlab's linesearch.m to be added to the path)
%
% For the Armijo line-search, several interpolation strategies are available ('LS_interp'):
%   - 0 : Step size halving
%   - 1 : Polynomial interpolation using new function values
%   - 2 : Polynomial interpolation using new function and gradient values (default)
%
% When (LS_interp = 1), the default setting of (LS_multi = 0) uses quadratic interpolation,
% while if (LS_multi = 1) it uses cubic interpolation if more than one point are available.
%
% When (LS_interp = 2), the default setting of (LS_multi = 0) uses cubic interpolation,
% while if (LS_multi = 1) it uses quartic or quintic interpolation if more than one point are available
%
% To use the non-monotonic Armijo condition, set the 'Fref' value to the number of previous function values to store
%
% For the Wolfe line-search, these interpolation strategies are available ('LS_interp'):
%   - 0 : Step Size Doubling and Bisection
%   - 1 : Cubic interpolation/extrapolation using new function and gradient values (default)
%   - 2 : Mixed quadratic/cubic interpolation/extrapolation
%
% Several strategies for choosing the initial step size are avaiable ('LS_init'):
%   - 0: Always try an initial step length of 1 (default for all except 'sd' and 'cg')
%       (t = 1)
%   - 1: Use a step similar to the previous step
%       (t = t_old*min(2,g'd/g_old'd_old)
%   - 2: Quadratic Initialization using previous function value and new
%   function value/gradient (use this if steps tend to be very long, default for 'sd' and 'cg')
%       (t = min(1,2*(f-f_old)/g))
%   - 3: The minimum between 1 and twice the previous step length
%       (t = min(1,2*t)
%   - 4: The scaled conjugate gradient step length (may accelerate
%   conjugate gradient methods, but requires a Hessian-vector product, default for 'scg')
%       (t = g'd/d'Hd)
%
% Inputs:
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
%
% Supported Input Options
%   Display - Level of display [ off | final | (iter) | full | excessive ]
%   MaxFunEvals - Maximum number of function evaluations allowed (1000)
%   MaxIter - Maximum number of iterations allowed (500)
%   optTol - Termination tolerance on the first-order optimality (1e-5)
%   progTol - Termination tolerance on progress in terms of function/parameter changes (1e-9)
%   Method - [ sd | csd | bb | cg | scg | pcg | {lbfgs} | newton0 | pnewton0 |
%       qnewton | mnewton | newton | tensor ]
%   c1 - Sufficient Decrease for Armijo condition (1e-4)
%   c2 - Curvature Decrease for Wolfe conditions (.2 for cg methods, .9 otherwise)
%   LS_init - Line Search Initialization - see above (2 for cg/sd, 4 for scg, 0 otherwise)
%   LS - Line Search type - see above (2 for bb, 4 otherwise)
%   Fref - Setting this to a positive integer greater than 1
%       will use non-monotone Armijo objective in the line search.
%       (20 for bb, 10 for csd, 1 for all others)
%   numDiff - [ 0 | 1 | 2] compute derivatives using user-supplied function (0),
%       numerically user forward-differencing (1), or numerically using central-differencing (2)
%       (default: 0) 
%       (this option has a different effect for 'newton', see below)
%   useComplex - if 1, use complex differentials if computing numerical derivatives
%       to get very accurate values (default: 0)
%   DerivativeCheck - if 'on', computes derivatives numerically at initial
%       point and compares to user-supplied derivative (default: 'off')
%   outputFcn - function to run after each iteration (default: []).  It
%       should have the following interface:
%       outputFcn(x,iterationType,i,funEvals,f,t,gtd,g,d,optCond,varargin{:});
%   useMex - where applicable, use mex files to speed things up (default: 1)
%
% Method-specific input options:
%   lbfgs:
%       Corr - number of corrections to store in memory (default: 100)
%           (higher numbers converge faster but use more memory)
%       Damped - use damped update (default: 0)

% Supported Output Options
%   iterations - number of iterations taken
%   funcCount - number of function evaluations
%   algorithm - algorithm used
%   firstorderopt - first-order optimality
%   message - exit message
%   trace.funccount - function evaluations after each iteration
%   trace.fval - function value after each iteration
% ===========================================================================================
% Get Parameters
[verbose,verboseI,debug,doPlot,maxFunEvals,maxIter,optTol,progTol,method,...
    corrections,c1,c2,LS_init,cgSolve,qnUpdate,cgUpdate,initialHessType,...
    HessianModify,Fref,useComplex,numDiff,LS_saveHessianComp,...
    Damped,HvFunc,bbType,cycle,...
    HessianIter,outputFcn,useMex,useNegCurv,precFunc,...
    LS_type,LS_interp,LS_multi,checkGrad] = ...
    minFunc_processInputOptions(options);

% Constants
SD = 0;
CSD = 1;
BB = 2;
CG = 3;
PCG = 4;
LBFGS = 5;
QNEWTON = 6;
NEWTON0 = 7;
NEWTON = 8;
TENSOR = 9;

% Initialize
p = length(x0);
d = zeros(p,1);
x = x0;
t = 1;
verbose = 1;
verboseI =1;
% If necessary, form numerical differentiation functions
funEvalMultiplier = 1;

numDiffType = numDiff;  %lbfgs 默认为0


% Evaluate Initial Point
if method < NEWTON   % 5
    [f,g] = funObj(x,varargin{:});  % 调用rosenblock函数...以x为输入...输出[f,df]
    computeHessian = 0;
end
funEvals = 1;

% Derivative Check   
% 0
if checkGrad   
	if numDiff
		fprintf('Can not do derivative checking when numDiff is 1\n');
		pause
	end
	derivativeCheck(funObj,x,1,numDiffType,varargin{:}); % Checks gradient
	if computeHessian
		derivativeCheck(funObj,x,2,numDiffType,varargin{:});
	end
end

% Output Log
% 0
if verboseI
    fprintf('%10s %10s %15s %15s %15s\n','Iteration','FunEvals','Step Length','Function Val','Opt Cond');
end

% Compute optimality of initial point
optCond = max(abs(g));   % max(abs(df))...梯度和  最大的

% Exit if initial point is optimal
if optCond <= optTol   % optTol=1e-5
    exitflag=1;
    msg = 'Optimality Condition below optTol';
%     if verbose     % verbose=0
%         fprintf('%s\n',msg);
%     end
    return; % 直接跳出函数
end

% Perform up to a maximum of 'maxIter' descent steps:
for i = 1:maxIter   % maxIter=500 

    % ****************** COMPUTE DESCENT DIRECTION *****************

    switch method
        case LBFGS % L-BFGS

            % Update the direction and step sizes
   % 0 不执行
			if ~Damped   
				if i == 1  % 第一次迭代 ..实现初始化操作
				% Initially use steepest descent direction
                    d = -g;  % g=df
					S = zeros(p,corrections);  % corrections=100 
					Y = zeros(p,corrections);  % p是输入起始向量的维度
					YS = zeros(corrections,1);
					lbfgs_start = 1;
					lbfgs_end = 0;
					Hdiag = 1;
                else      % 第2-500次迭代
					[S,Y,YS,lbfgs_start,lbfgs_end,Hdiag,skipped] = lbfgsAdd(g-g_old,t*d,S,Y,YS,lbfgs_start,lbfgs_end,Hdiag,useMex);
					if debug && skipped
						fprintf('Skipped L-BFGS updated\n');
					end
					if useMex
						d = lbfgsProdC(g,S,Y,YS,int32(lbfgs_start),int32(lbfgs_end),Hdiag);
					else
						d = lbfgsProd(g,S,Y,YS,lbfgs_start,lbfgs_end,Hdiag);
					end
				end
			end
			g_old = g;
    end

    % 当输入向量中不含0,不含非参数,不含无穷大值,为合法
    if ~isLegal(d)  
        fprintf('Step direction is illegal!\n');
        pause;
        return
    end

    % ****************** COMPUTE STEP LENGTH ************************

    % Directional Derivative
    gtd = g'*d;  % 第一次迭代计算的是 -df'*df== --*|  ...即负的梯度模值平方   -|df|^2

    % Check that progress can be made along direction
    if gtd > -progTol   % progTol=1e-9
        exitflag=2;
        msg = 'Directional Derivative below progTol';
        break;
    end

    % Select Initial Guess
    if i == 1    % 第一次迭代
        if method < NEWTON0   % method=5
            t = min(1,1/sum(abs(g)));   % 执行   
        else                       
            t = 1;                      % 不执行
        end
    else        % 第2-500次迭代
        if LS_init == 0
            % Newton step
            t = 1;
        elseif LS_init == 1
            % Close to previous step length
            t = t*min(2,(gtd_old)/(gtd));
        elseif LS_init == 2
            % Quadratic Initialization based on {f,g} and previous f
            t = min(1,2*(f-f_old)/(gtd));
        elseif LS_init == 3
            % Double previous step length
            t = min(1,t*2);
        elseif LS_init == 4
            % Scaled step length if possible
            if isempty(HvFunc)
                % No user-supplied Hessian-vector function,
                % use automatic differentiation
                dHd = d'*autoHv(d,x,g,0,funObj,varargin{:});
            else
                % Use user-supplid Hessian-vector function
                dHd = d'*HvFunc(d,x,varargin{:});
            end

            funEvals = funEvals + 1;
            if dHd > 0
                t = -gtd/(dHd);
            else
                t = min(1,2*(f-f_old)/(gtd));
            end
        end

        if t <= 0
            t = 1;
        end
    end
    f_old = f;   % f
    gtd_old = gtd; %df模值平方

    % Compute reference fr if using non-monotone objective
    if Fref == 1   % Fref=1 ,执行
        fr = f;
    else          %不执行
        if i == 1
            old_fvals = repmat(-inf,[Fref 1]);
        end

        if i <= Fref
            old_fvals(i) = f;
        else
            old_fvals = [old_fvals(2:end);f];
        end
        fr = max(old_fvals);
    end

    computeHessian = 0;
    % Line Search
    f_old = f;
    
    if LS_type == 1 % Find Point satisfying Wolfe conditions  % 执行
        if ~computeHessian   % computeHessian=0
            %  执行
              %  WolfeLineSearch 函数用于求取步长
              % Outputs:
              %   t: step length
              %   f_new: function value at x+t*d
              %   g_new: gradient value at x+t*d
              %   funEvals: number function evaluations performed by line search
           [t,f,g,LSfunEvals] = WolfeLineSearch(x,t,d,f,g,gtd,c1,c2,LS_interp,LS_multi,25,progTol,debug,doPlot,1,funObj,varargin{:});
        end
        funEvals = funEvals + LSfunEvals;
        x = x + t*d;
    end

	% Compute Optimality Condition
	optCond = max(abs(g));
	
    % Output iteration information
    if verboseI  % 0 不执行
        fprintf('%10d %10d %15.5e %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,t,f,optCond);
    end
	
    % Check Optimality Condition
    if optCond <= optTol  % optTol=1e-5 
        exitflag=1;
        msg = 'Optimality Condition below optTol';
        break;
    end

    % ******************* Check for lack of progress *******************

    if max(abs(t*d)) <= progTol % progTol=1e-9
        exitflag=2;
        msg = 'Step Size below progTol';
        break;
    end


    if abs(f-f_old) < progTol
        exitflag=2;
        msg = 'Function Value changing by less than progTol';
        break;
    end

    % ******** Check for going over iteration/evaluation limit *******************

    if funEvals*funEvalMultiplier >= maxFunEvals  % maxFunEvals=1000
        exitflag = 0;
        msg = 'Reached Maximum Number of Function Evaluations';
        break;
    end

    if i == maxIter
        exitflag = 0;
        msg='Reached Maximum Number of Iterations';
        break;
    end

end

if verbose  % verbose=0
    fprintf('%s\n',msg);
end

end

