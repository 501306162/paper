function [ki_up, dki_up] =Ed_derivation1(ff,fm,kx,rho,Dissimilarity_Metric)
% k-update 把k看成一个向量来计算,求出一个点的位移场更新  
% 对一个向量进行位移场更新...返回找到最小值的点...即最小 f 对应的 k..
%   ki<->xi  ...ki+1 <->xi+1
% This function returns the function value, partial derivatives
% and Hessian of the (general dimension)  function, given by:
%
%  rosenbrock:f(x) = sum_{i=1:D-1} 100*(x(i+1) - x(i)^2)^2 + (1-x(i))^2 
% where D is the dimension of x. The true minimum is 0 at x = (1 1 ... 1).
%  
%  F(kx)=... F(ky)=... F(kz)=... 

 最粗网格点
    Kmax=Mpyr*K;
    for j=1:Mpyr
        Kn=j.*K;
        f_m=fm(1:Kn(1):end,1:Kn(2):end);
        f_f=ff(1:Kn(1):end,1:Kn(2):end);
    end



%  x 表示 控制点 or 像素点 标号
    D = length(kx);
    % spacing : pixel (voxel)'s physical dimensions δn ,n=1,2...N
    v=prod(spacing);
    % the first term: 返回值对应原图像的维度
    switch Dissimilarity_Metric        
        case 'SSD'
             [Metric, dMetric] =SSD(f_f,f_m,v);
        case 'LCC'
             [Metric, dMetric] =LCC(f_f,f_m,v);
    end
    
    % the second term:
    [dfm_x,dfm_y,dfm_z]=Gradient(f_m);
  
    % 第三项是对第i个控制点求导还是对控制点的第i个方向求导?
    % the third term:
    % p和x都表示坐标位置
    dk=prod(max((1-abs(p-x)./K),0));
    
邻域   
    neighborhood_trans(ind,image,K); % K1,K2,K3大小邻域...中的点对应在原图中标号
    
    ki_up=Metric+(rho/2).*norm(D(p,N,f_m)-(Z-U),'fro');
    
    ed_up=dMetric.*dfm_x.*dk+dMetric.*dfm_y.*dk+dMetric.*dfm_z.*dk;
    
    dki_up=ed_up+;
end


%  以下函数输入均为图像,返回也为图像
function [df_x,df_y,df_z]=Gradient(f)
% 输入f为图像,返回也为图像
    if numel(size(f))<3
        df_z=0;
        [df_x,df_y]=gradient(f_m);
    else
        [df_x,df_y,df_z]=gradient(f_m);
    end
end

function [f, df] =SSD(f_f,f_m,v)    
% 输入图像为原图大小,返回也为图像 
%     if numel(size(f_f))>2
        ssd=sum(sum(sum((v/2)*((f_m-f_f).^2))));        
%     else
%         ssd=sum(sum((v/2)*((f_m-f_f).^2)));
%     end
    f=ssd;
    if nargout > 1
      df=(f_m-f_f).*v;
%       df=reshape(df,length(f_f(:)),1);
    end
end
function [f, df] =LCC(f_f,f_m,v)
%    sigma f_f:...感觉还需要reshape一下...f_f(x)...按每层维度
    windows=[15 15];
    w=2.5;
    Hw=fspecial('gaussian',windows,w);
    
    f_bar=imfilter(f_f,Hw,'conv'); 
    m_bar=imfilter(f_m,Hw,'conv');
    f2_bar=imfilter(f_f.^2,Hw,'conv');
    m2_bar=imfilter(f_m.^2,Hw,'conv');
    
    sigma_f=sqrt(f2_bar-f_bar.^2);
    sigma_m=sqrt(m2_bar-m_bar.^2);
    
    multiplication=imfilter(f_m.*f_f,Hw,'conv')-m_bar.*f_bar;
    
%     if numel(size(f_f))>2   % 其实可以直接用三次sum(3D)...这样也包含2D下的求和操作
        lcc=-sum(sum(sum((v)*(multiplication./(sigma_m.*sigma_f)))));
%     else
%         lcc=-sum(sum((v)*(multiplication./(sigma_m.*sigma_f))));
%     end   
    f=lcc;
    if nargout > 1
      df=(-sigma_m.*sigma_m).*((f_f-f_bar)-(f_m-m_bar).*multiplication./(sigma_m.^2));
%       df=reshape(df,length(f_f(:)),1);
    end
end

