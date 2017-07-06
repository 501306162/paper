function [ki_up, dki_up] =Ed_2d(k,ffo,fmo,rho,Dissimilarity,Spacing,ZU_diff)
% k-update ��k����һ������������,���һ�����λ�Ƴ�����  
% ��һ����������λ�Ƴ�����...�����ҵ���Сֵ�ĵ�...����С f ��Ӧ�� k..
%  
%     ���������һά������...���Ҳ��һά������
%   k    ��   n*L�� 1��...һ��3ά/
%   Dk��u��z �� n2�� L��

    [m,n]=size(fmo);
    [kx,ky]=ndgrid(0:Spacing(1):m,0:Spacing(2):n); % ��������
    [x,y]=ndgrid(0:m-1,0:n-1); % ��������    
    
    v=prod(Spacing);
    k_m=size(kx,1);  
    k_n=size(kx,2);
    num_control=k_m*k_n;     
    % Construct the k_grid

    % ����λ�ƾ���
    k_grid=reshape(k,k_m,k_n,2); % ��һƬΪTr,�ڶ�ƬΪTc
    k_grid=double(k_grid);
    
    % ������λ�ƻ�õ����ص�λ�ƾ����
    dk=bspline_transform(k_grid,Spacing,[m,n]);

    % �ƶ��Ժ�ͼ�� ...
    f_m=movepixels_2d(fmo,dk(:,:,1),dk(:,:,2));
   
    % ��һ�����  ����ֵ��ͼ��ͬά��..
    switch Dissimilarity        
        case 'SSD'
             [Metric, dMetric] =SSD(ffo,f_m,v);  % <-��ʱͼ��Ϊ������ĳһ���ͼ�� ,v Ϊ�ò����ش�С
        case 'LCC'
             [Metric, dMetric] =LCC(ffo,f_m,v);
    end
    
    % �ڶ���  Ҳ��ͼ��ά��
    [dfm_c,dfm_r,dfm_z]=Gradient(f_m);  
    
    % ����ֵ һ����
    ki_up=Metric+(rho/2).*norm(D_new(k_grid(:,:,1),k_grid(:,:,2))-ZU_diff,'fro');
    
    % ��������ֵ 
    % ������  ���ص�ǰ���Ƶ��������,�Լ�ddkֵ       

    k_up_1r=zeros(1,num_control);
    k_up_1c=zeros(1,num_control);
    
    for i=1:num_control
        k_p=find(x==kx(i)&y==ky(i));    % �ҵ����Ƶ������ؾ����е�λ��
        ddk=neighborhood_trans(k_p,[m n],Spacing);       
        neigh_ind=ddk(:,1);       
        
        k_up_1r(i)=sum(dMetric(neigh_ind).*dfm_r(neigh_ind).*ddk(:,2));
        k_up_1c(i)=sum(dMetric(neigh_ind).*dfm_c(neigh_ind).*ddk(:,2));        
    end   

    k_up_1=[k_up_1r;k_up_1c]; % 2�� num_control��
   
    temp=D_new(k_grid(:,:,1),k_grid(:,:,2))-ZU_diff;  
    k_up_2=rho.*D_conju(temp,k_m,k_n); 
    
    % Ӧ�����ά��Ϊ (n*L,1) 
    dki_up=(k_up_1+k_up_2)';
    dki_up=dki_up(:);
end

%  ���º��������Ϊͼ��,����ҲΪͼ��
function [df_c,df_r,df_z]=Gradient(f)
% ����fΪͼ��,����ҲΪͼ��
    if numel(size(f))<3
        df_z=0;
        [df_c,df_r]=gradient(f);
    else
        [df_c,df_r,df_z]=gradient(f);
    end
end

function [f, df] =SSD(f_f,f_m,v)    
% ����ͼ��Ϊԭͼ��С,����ҲΪͼ�� 
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

