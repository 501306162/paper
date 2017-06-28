function [ki_up, dki_up] =Ed_derivation2_2d(ffo,fmo,k,rho,Dissimilarity,spacing,ddk,ZU_diff)
% k-update 把k看成一个向量来计算,求出一个点的位移场更新  
% 对一个向量进行位移场更新...返回找到最小值的点...即最小 f 对应的 k..
%  
%     输入最好是一维列向量...输出也是一维列向量


    % spacing : pixel (voxel)'s physical dimensions δn ,n=1,2...N
    v=prod(spacing);
    kr=reshape(k(1:numel(fmo)),size(fmo));
    kc=reshape(k(numel(fmo)+1:end),size(fmo));
    [m,n]=size(ffo);

    % 移动以后图像 ...根据demons算法里的程序进行修改
    f_m=movepixels_2d(fmo,kr,kc);
   
    % 第一项操作  移动以后的  网格图像...
    switch Dissimilarity        
        case 'SSD'
             [Metric, dMetric] =SSD(ffo,f_m,v);  % <-此时图像为金字塔某一层的图像 ,v 为该层像素大小
        case 'LCC'
             [Metric, dMetric] =LCC(ffo,f_m,v);
    end
    
    % 第二项    移动以后的 也是网格图像
    [dfm_y,dfm_x,dfm_z]=Gradient(f_m);  
    
    % 第三项    是网格坐标和原图坐标    
    % 要计算网格点在原图的位置和之前像素点的位置  
    % p和x都表示坐标位置   ...一个网格像素的空间内才有dk
    % 感觉该项可以直接传入....可在外层函数内计算
    %     dk_x=max((1-abs(px_diff(1))./K(1)),0);
    %     dk_y=max((1-abs(px_diff(2))./K(2)),0);
    %     ddk=dk_x.*dk_y;
    %     dk_z=max((1-abs(px_diff)./K(3)),0);

    ki_up=Metric+(rho/2).*norm(D_new(kr,kc)-ZU_diff,'fro');
    
    k_up_1r=dMetric.*dfm_x.*ddk;
    k_up_1c=dMetric.*dfm_y.*ddk;
%     k_up_1z=dMetric.*dfm_z.*ddk;
    

    % D_conju_new  返回值第一行为第一维度,第二行对应第二维度,都为L列
    % D_new   返回值为[d11;d12;d21;d22]...其中dii为L维列向量...故ZU_diff也应该是这个维度
 怎么调整？    
    conju_in=D_new(kr,kc)-ZU_diff;
    conju_in_r=
    conju_in_c=
    k_up_2=rho.*D_conju_new(conju_in_r,conju_in_c); 
    
    
end

%  以下函数输入均为图像,返回也为图像
function [df_y,df_x,df_z]=Gradient(f)
% 输入f为图像,返回也为图像
    if numel(size(f))<3
        df_z=0;
        [df_y,df_x]=gradient(f);
    else
        [df_y,df_x,df_z]=gradient(f);
    end
end

function [f, df] =SSD(f_f,f_m,v)    
% 输入图像为原图大小,返回也为图像 
    ssd=sum(sum(sum((v/2)*((f_m-f_f).^2))));        
    f=ssd;
    if nargout > 1
      df=(f_m-f_f).*v;
   end
end
function [f, df] =LCC(f_f,f_m,v)
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
    
    lcc=-sum(sum(sum((v)*(multiplication./(sigma_m.*sigma_f))))); 
    f=lcc;
    if nargout > 1
      df=(-sigma_m.*sigma_m).*((f_f-f_bar)-(f_m-m_bar).*multiplication./(sigma_m.^2));
    end
end

