function [f, df] =dissimilarity(f_f,f_m,x,Dissimilarity_Metric)
    %  x 表示 控制点 or 像素点 标号
    D = length(x);
    % spacing : pixel (voxel)'s physical dimensions δn ,n=1,2...N
    v=prod(spacing);
    switch Dissimilarity_Metric        
        case 'SSD'
             [f, df] =SSD(f_f,f_m,v);
        case 'LCC'
             [f, df] =LCC(f_f,f_m,v);
    end
  
end
function [f, df] =SSD(f_f,f_m,v)    
% 输入图像为原图大小
    ssd=sum(sum((v/2)*((f_m-f_f).^2)));
    f=ssd;
    if nargout > 1
      df=(f_m-f_f).*v;
      df=reshape(df,length(f_f(:)),1);
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
      df=reshape(df,length(f_f(:)),1);
    end
end