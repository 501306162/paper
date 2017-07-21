function [ki_up, dki_up] =Ed_2d(k,ffo,fmo,rho,Dissimilarity,Spacing,delta,ZU_diff)
% k-update 把k看成一个向量来计算,求出一个点的位移场更新  
% 对一个向量进行位移场更新...返回找到最小值的点...即最小 f 对应的 k..
%  
%     输入最好是一维列向量...输出也是一维列向量
%   k    是   n*L行 1列...一共3维/
%   Dk、u、z 是 n2列 L行

    [m,n]=size(fmo);
    [kx,ky]=ndgrid(0:Spacing(1):m-1,0:Spacing(2):n-1); % 网格坐标
    [x,y]=ndgrid(0:m-1,0:n-1); % 像素坐标    
    
    v=prod(delta);
    k_m=size(kx,1);  
    k_n=size(kx,2);
    num_control=k_m*k_n;     
    % Construct the k_grid

    % 网格位移矩阵
    k_grid=reshape(k,k_m,k_n,2); % 第一片为Tr,第二片为Tc
    k_grid=double(k_grid);
%  fprintf('minFunc in  BFGS :k\n');
% k(:)'
% fprintf('---------------------------------------\n');   
    % 由网格位移获得的像素点位移矩阵
    dk=bspline_transform(k_grid,Spacing,[m,n]);
%  fprintf('minFunc in  BFGS :dk\n');
% dk(:)'
% fprintf('---------------------------------------\n');   
    % 移动以后图像 ...
    f_m=movepixels_2d(fmo,dk(:,:,1),dk(:,:,2));
    figure,
    subplot(1,2,1), imshow(f_m); title('浮动图');
    subplot(1,2,2), imshow(ffo); title('固定图');
    
    % 第一项操作  返回值与图像同维度..
    switch Dissimilarity        
        case 'SSD'
             [Metric, dMetric] =SSD(ffo,f_m,v);  % <-此时图像为金字塔某一层的图像 ,v 为该层像素大小
        case 'LCC'
             [Metric, dMetric] =LCC(ffo,f_m,v);
    end
    
    % 第二项  也是图像维度
    [dfm_c,dfm_r,dfm_z]=Gradient(f_m);  
    
    % 函数值 一个数
    ki_up=Metric+(rho/2).*norm(D_new(k_grid(:,:,1),k_grid(:,:,2))-ZU_diff,'fro');
    
    % 函数导数值 
    % 第三项  返回当前控制点的邻域标号,以及ddk值       

    k_up_1r=zeros(1,num_control);
    k_up_1c=zeros(1,num_control);
    
    for i=1:num_control
        k_p=find(x==kx(i)&y==ky(i));    % 找到控制点在像素矩阵中的位置
        
%         返回了空k_p值
        ddk=neighborhood_trans(k_p,[m n],Spacing);       
        neigh_ind=ddk(:,1);       
        
        k_up_1r(i)=sum(dMetric(neigh_ind).*dfm_r(neigh_ind).*ddk(:,2));
        k_up_1c(i)=sum(dMetric(neigh_ind).*dfm_c(neigh_ind).*ddk(:,2));        
    end   

    k_up_1=[k_up_1r;k_up_1c]; % 2行 num_control列
   
    temp=D_new(k_grid(:,:,1),k_grid(:,:,2))-ZU_diff;  
    k_up_2=rho.*D_conju(temp,k_m,k_n); 
    
    % 应该输出维度为 (n*L,1) 
    dki_up=(k_up_1+k_up_2)';
    dki_up=dki_up(:);
end

%  以下函数输入均为图像,返回也为图像
function [df_c,df_r,df_z]=Gradient(f)
% 输入f为图像,返回也为图像
    if numel(size(f))<3
        df_z=0;
        [df_c,df_r]=gradient(f);
    else
        [df_c,df_r,df_z]=gradient(f);
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

