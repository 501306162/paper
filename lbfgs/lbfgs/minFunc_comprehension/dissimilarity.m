function [f, df] =dissimilarity(f_f,f_m,x,Dissimilarity_Metric)
    %  x ��ʾ ���Ƶ� or ���ص� ���
    D = length(x);
    % spacing : pixel (voxel)'s physical dimensions ��n ,n=1,2...N
    v=prod(spacing);
    switch Dissimilarity_Metric        
        case 'SSD'
             [f, df] =SSD(f_f,f_m,v);
        case 'LCC'
             [f, df] =LCC(f_f,f_m,v);
    end
  
end
function [f, df] =SSD(f_f,f_m,v)    
% ����ͼ��Ϊԭͼ��С
    ssd=sum(sum((v/2)*((f_m-f_f).^2)));
    f=ssd;
    if nargout > 1
      df=(f_m-f_f).*v;
      df=reshape(df,length(f_f(:)),1);
    end
end
function [f, df] =LCC(f_f,f_m,v)
%    sigma f_f:...�о�����Ҫreshapeһ��...f_f(x)...��ÿ��ά��
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