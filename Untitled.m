function [d, history] = Untitled(A, b, mu, rho0, alpha)
% solves the following problem via ADMM:%
%   minimize   E_D(d(k);f_f,f_m)+lamda*eta*||Z||_2,1
%   s.t.       D(k)=Z;
% where A is a feature matrix and b is a response vector. The scalar m is
% the number of examples in the matrix A.
% ========================================================
% data set:
% f_m & f_f : images
% lamda     : regularization    0.1
% w         : LCC kernel bandwith   2.5
% M_pyr      : # of image pyramid levels  log2(m/8)
% M_ref      : # of grid refinements
% M_iter     : max.# of ADMM iterations
% M_lbgfs   : max.# of LBFGS iterations   5
% kesai_tol :optimization argument tolerance(10^-2)
% 
% ADMM internal parameters : ro_0=5;mu=10;tao=1.3;
% =========================================================
t_start = tic;
% Global constants and defaults
f_m=randi(10,10,10);
f_f=randi(10,10,10);
[m,n]=size(f_m);
lamda=0.1;
rho_0=5;
mu=10;
tao=1.3;

w=2.5;
M_pyr=ceil(log2(m/8));
M_ref=5;
M_iter=5;
M_lbgfs=5;
k(1,0)=0;

QUIET    = 0;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
kesai_tol= 1e-2;

% Data preprocessing

[m, n] = size(A);

% M_pyr sovler
for j=1:M_pyr+M_ref  
    if j<M_pyr
        %%
        f_mj=f_m(M_pyr-j)  %(Mpyr - j )-th Gaussian pyramid level of fm
        f_fj=f_f(M_pyr-j)  %(Mpyr - j )-th Gaussian pyramid level of ff
    else
        f_mj=f_mj_old;
        f_fj=f_fj_old;
    end
    
    num_control=numel(f_mj);
    % ADMM solver in j-th level 
    k = zeros(num_control+1,1);
    z = Dk(k); 
    u = zeros(num_control+1,1);
    rho=rho_0;

    for m_iter = 1:M_iter
        k_old=k;
        
        %% x-update use lbfgs
        k = update_x(A, b, u, z, rho);  ' undone!'
                
        kk=sum(k-k_old);% max {row(j)_sum}
        if kk<kesai_tol
            break;
        end
        %%

        %% z-update with relaxation && use group lasso
        z_old = z;
        x_hat = alpha*k + (1-alpha)*z_old;
       
        z=update_z(A, b, x_hat, u, rho);         ' undone!'
        %%
        % u-update
        u = u + (x_hat - z);
        % diagnostics, reporting, termination checks
        history.objval(m_iter)  = objective(A, b, mu, k, z);
        history.eps_pri(m_iter) = sqrt(n)*ABSTOL + RELTOL*max(norm(k), norm(z));
        history.eps_dual(m_iter)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
        history.r_norm(m_iter)=norm(Dk(k)-z);
        history.s_norm(m_iter)=norm(rho*D_conjugate(z , z_old));
        % varying penalty 
        if history.r_norm(m_iter)>=mu*history.s_norm(m_iter)
            rho=rho*tao;
            u=u/tao;
        end
        if history.s_norm(m_iter)>=mu*history.r_norm(m_iter)
            rho=rho/tao;
            u=u*tao;
        end
    end % end of M_iter

    if ~QUIET
        toc(t_start);
    end
      
end % end of M_pyr

end % end of function
function obj = objective(A, b, mu, x, z)
    m = size(A,1);
    obj = sum(log(1 + exp(-A*x(2:end) - b*x(1)))) + m*mu*norm(z,1);
end

function x = update_x(A, b, u, z, rho, x0)
    % solve the x update
    %   minimize [ -logistic(x_i) + (rho/2)||x_i - z^k + u^k||^2 ]
    % via Newton's method; for a single subsystem only.
    alpha = 0.1;
    BETA  = 0.5;
    TOLERANCE = 1e-5;
    MAX_ITER = 50;
    [m n] = size(A);
    I = eye(n+1);
    if exist('x0', 'var')
        x = x0;
    else
        x = zeros(n+1,1);
    end
    C = [-b -A];
    f = @(w) (sum(log(1 + exp(C*w))) + (rho/2)*norm(w - z + u).^2);
    for iter = 1:MAX_ITER
        fx = f(x);
        g = C'*(exp(C*x)./(1 + exp(C*x))) + rho*(x - z + u);
        H = C' * diag(exp(C*x)./(1 + exp(C*x)).^2) * C + rho*I;
        dx = -H\g;   % Newton step
        dfx = g'*dx; % Newton decrement
        if abs(dfx) < TOLERANCE
            break;
        end
        % backtracking
        t = 1;
        while f(x + t*dx) > fx + alpha*t*dfx
            t = BETA*t;
        end
        x = x + t*dx;
    end
end

%% z-updata
function z = shrinkage(a, kappa)
    z = max(0, a-kappa) - max(0, -a-kappa);
end
function z = update_z(A, b, x_hat, u, rho)    

        z = x_hat + u;
        z = shrinkage(z, (lamda*eta)/rho);
end

function z = Dk(k)

end
function s=D_conjugate(z , z_old)

end

function z = D(d,N)
     con_distance=1; % 控制点间距
     % forward difference :delta_i(dj[l])
     for j=1:N
         for i=1:N
             % displacement in j-th dimension :d(:,j);
                dis(r,c)=d(:,j);
             % difference of j-th dispalcement,in i-th direction   
     ' undone!'   维度不一致
                d_r=(dis(r+1,c)-dis(r,c))/con_distance;
                delta(1:length(d_r),j)=d_r(:);
                
                d_c=(dis(r,c+1)-dis(r,c))/con_distance;
                delta((length(d_r)+1):(length(d_r)+length(d_c)),j)=d_c(:);
                
         end
     end
end
%%