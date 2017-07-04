function [ki_up, dki_up] =Ed_derivation2_2d_new(k,ffo,fmo,rho,Dissimilarity,spacing,ZU_diff)
% k-update ��k����һ������������,���һ�����λ�Ƴ�����  
% ��һ����������λ�Ƴ�����...�����ҵ���Сֵ�ĵ�...����С f ��Ӧ�� k..
%  
%     ���������һά������...���Ҳ��һά������
%   k    �� 1��  n*L��...һ��3ά/
%   Dk��u��z �� n2�� L��
    [m,n]=size(fmo);
    control_row=1:spacing(1):m;
    control_col=1:spacing(2):n;
    control_row=repmat(control_row',1,numel(control_col));
    control_col=repmat(control_col,size(control_row,1),1);
    
    control_ind=sub2ind(size(fmo),control_row(:),control_col(:));
    for i=1:numel(control_ind)
     neigh_ddk=neighborhood_trans(control_ind(i),fmo,spacing);
    end
    
    num_control=numel(control_ind); 

    % spacing : pixel (voxel)'s physical dimensions ��n ,n=1,2...N
    v=prod(spacing);
    kr=reshape(k(1:num_control),m,n);
    kc=reshape(k(num_control+1:end),m,n);


    % �ƶ��Ժ�ͼ�� ...����demons�㷨��ĳ�������޸�
    f_m=movepixels_2d(fmo,kr,kc);
   
    % ��һ�����  �ƶ��Ժ��  ����ͼ��...
    switch Dissimilarity        
        case 'SSD'
             [Metric, dMetric] =SSD(ffo,f_m,v);  % <-��ʱͼ��Ϊ������ĳһ���ͼ�� ,v Ϊ�ò����ش�С
        case 'LCC'
             [Metric, dMetric] =LCC(ffo,f_m,v);
    end
    
    % �ڶ���    �ƶ��Ժ�� Ҳ������ͼ��
    [dfm_c,dfm_r,dfm_z]=Gradient(f_m);  
    
    % ������    �����������ԭͼ����    
    % Ҫ�����������ԭͼ��λ�ú�֮ǰ���ص��λ��  
    % p��x����ʾ����λ��   ...һ���������صĿռ��ڲ���dk
    % �о��������ֱ�Ӵ���....������㺯���ڼ���
    %     dk_x=max((1-abs(px_diff(1))./K(1)),0);
    %     dk_y=max((1-abs(px_diff(2))./K(2)),0);
    %     ddk=dk_x.*dk_y;
    %     dk_z=max((1-abs(px_diff)./K(3)),0);

    
    ki_up=Metric+(rho/2).*norm(D_new(kr,kc)-ZU_diff,'fro');
    
    k_up_1r=dMetric.*dfm_r.*ddk;
    k_up_1c=dMetric.*dfm_c.*ddk;
%     k_up_1z=dMetric.*dfm_z.*ddk;
    k_up_1=[k_up_1r(:)';k_up_1c(:)'];

    % D_conju_new  ����Ϊdr,dc��Ӧͼ�����...����ֵ��һ��Ϊ��һά��,�ڶ��ж�Ӧ�ڶ�ά��,��ΪL��
    % D_new   ����ֵΪ[d11;d12;d21;d22]...����diiΪLά������...��ZU_diffҲӦ�������ά��
    
    temp=D_new(kr,kc)-ZU_diff;  
    k_up_2=rho.*D_conju_new2(temp,num_control,m,n); 
    
    % Ӧ�����ά��Ϊ (1,n*L) 
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

