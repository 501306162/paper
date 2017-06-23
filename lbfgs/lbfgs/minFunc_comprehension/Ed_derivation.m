function [ki_up, dki_up] =Ed_derivation(f_f,f_m,x,Dissimilarity_Metric)
% k-update 一次计算K的一个维度对应图像产生的所有值..计算2-3次,    
%  x 表示 控制点 or 像素点 标号
    D = length(x);
    % spacing : pixel (voxel)'s physical dimensions δn ,n=1,2...N
    v=prod(spacing);
    % the first term:
    switch Dissimilarity_Metric        
        case 'SSD'
             [Metric, dMetric] =SSD(f_f,f_m,v);
        case 'LCC'
             [Metric, dMetric] =LCC(f_f,f_m,v);
    end
    
    % the second term:
    [dfm_x,dfm_y,dfm_z]=gradient(f_m);
  
    % the third term:
    dk=prod(max((1-abs(p-x)./K),0));
    % 第三项是对第i个控制点求导还是对控制点的第i个方向求导?
    ki_up=Metric+(rho/2).*norm(D(p,N,f_m)-(Z-U),'fro');
    
    ed_up=dMetric.*dfm_x.*dk+dMetric.*dfm_y.*dk+dMetric.*dfm_z.*dk;
    
    dki_up=ed_up+;
end
function [f, df] =SSD(f_f,f_m,v)    
% 输入图像为原图大小
    ssd=sum(sum((v/2)*((f_m-f_f).^2)));
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
    
    lcc=-sum(sum((v)*(multiplication./(sigma_m.*sigma_f))));
    
    f=lcc;
    if nargout > 1
      df=(-sigma_m.*sigma_m).*((f_f-f_bar)-(f_m-m_bar).*multiplication./(sigma_m.^2));
%       df=reshape(df,length(f_f(:)),1);
    end
end

