function [d, history] = isopTV_new(A, b, p,mu, rho0, alpha)
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
N=numel(size(f_m));
lamda=0.1;
rho_0=5;
mu=10;
tao=1.3;

w=2.5;
M_pyr=ceil(log2(m/8));
M_ref=5;
M_iter=5;
M_lbgfs=5;

Dissimilarity='LLC';
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
    
    [m, n] = size(f_mj);
    % ADMM solver in j-th level 
    Spacing=[4 4];
    km=length(0:Spacing(1):m);
    kn=length(0:Spacing(2):n);
    num_control=km*kn;
    
    k = zeros(num_control,1);
    k_grid=reshape(k,km,kn,2);       
    z = D(k_grid(:,:,1),k_grid(:,:,2)); 
    u = zeros(sizze(z));
    rho=rho_0;     
    
    % k   �� (n*L,1)ά
    % z\u ��(L,N^2)'ά
    for m_iter = 1:M_iter
        k_old=k;
        % x-update use lbfgs
        % LBFGS_k ����� (n*L,1)ά
        k=LBFGS_k(k_old,f_fj,f_mj,rho,Dissimilarity,Spacing,z-u );
        k_grid=reshape(k,km,kn,2);     
        % 'break' condition
        kk=norm((k-k_old),'inf');% max {row(j)_abs_sum}
        if kk>=kesai_tol
        else
            break;
        end
        % z-update with relaxation && use group lasso
        % update with k->D(k)
        z_old = z;
        Dk=D(k_grid(:,:,1),k_grid(:,:,2));
        x_hat = alpha.*Dk + (1-alpha).*z_old;
        z=update_z(lamda,eta, x_hat, u, rho,p,z_old);         'to check!'
        
        % u-update
        u = u + (x_hat - z);
        % diagnostics, reporting, termination checks
        history.objval(m_iter)  = objective(A, b, mu, k, z);
        history.eps_pri(m_iter) = sqrt(n)*ABSTOL + RELTOL*max(norm(k), norm(z));
        history.eps_dual(m_iter)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
        
        % varying penalty 
        history.r_norm(m_iter)=norm(Dk-z);
        history.s_norm(m_iter)=norm(rho*D_conjugate(z , z_old));
        
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



%% z-updata

%%
% z-update
function z_out = update_z(lamda,eta, x_hat, u, rho,p,z_old) 
    % cumulative partition
    cum_part = cumsum(p);
    start_ind = 1;
    for i = 1:length(p),
        sel = start_ind:cum_part(i);
        z_old(sel) = shrinkage(x_hat(sel) + u(sel), lamda*eta/rho);
        start_ind = cum_part(i) + 1;
    end
    z_out=z_old;
end
function z = shrinkage(a, kappa)
    z = max(0, a-kappa) - max(0, -a-kappa);
end
% function z = shrinkage(a, kappa)
%     z = pos(1 - kappa/norm(a))*a;
% end
% the function D(d) used in TV-Regularization
function Dd = D(dr,dc)  %2-D
      [dr_dc,dr_dr]=gradient(dr);
      [dc_dc,dc_dr]=gradient(dc);
      
      Dd=[dr_dr(:)';dr_dc(:)';dc_dr(:)';dc_dc(:)'];
end