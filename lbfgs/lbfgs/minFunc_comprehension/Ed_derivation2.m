function [ki_up, dki_up] =Ed_derivation2(ffo,fmo,k,rho,Dissimilarity)
% k-update ��k����һ������������,���һ�����λ�Ƴ�����  
% ��һ����������λ�Ƴ�����...�����ҵ���Сֵ�ĵ�...����С f ��Ӧ�� k..
%   ki<->xi  ...ki+1 <->xi+1
% This function returns the function value, partial derivatives
% and Hessian of the (general dimension)  function, given by:
%
%  rosenbrock:f(x) = sum_{i=1:D-1} 100*(x(i+1) - x(i)^2)^2 + (1-x(i))^2 
% where D is the dimension of x. The true minimum is 0 at x = (1 1 ... 1).
%  
%  F(kx)=... F(ky)=... F(kz)=... 

%  ��������
%     Kmax=Mpyr*K;
%     for j=1:Mpyr
%         Kn=j.*K;
%         f_m=fm(1:Kn(1):end,1:Kn(2):end);
%         f_f=ff(1:Kn(1):end,1:Kn(2):end);
%     end

    
%  x ��ʾ ���Ƶ� or ���ص� ���
    % spacing : pixel (voxel)'s physical dimensions ��n ,n=1,2...N
    v=prod(spacing);
    kx=k(:,:1);
    ky=k(:,:2);

�ƶ��Ժ�ͼ�� ...���Ը���demons�㷨��ĳ�������޸�
    if(numel(k)>(2*numel(ffo)+1)) % 2d
        f_m=movepixels_2d(double(fmo(:,:)),double(kx),double(ky));
    else    
        f_m=movepixels_3d(fmo(:,:),kx,ky); % 3d
    end   

��һ�����  �ƶ��Ժ��  ����ͼ��...
    % the first term: ����ֵ��Ӧԭͼ���ά��
    switch Dissimilarity        
        case 'SSD'
             [Metric, dMetric] =SSD(f_f,f_m,v);  % <-��ʱͼ��Ϊ������ĳһ���ͼ�� ,v Ϊ�ò����ش�С
        case 'LCC'
             [Metric, dMetric] =LCC(f_f,f_m,v);
    end
�ڶ���    �ƶ��Ժ�� Ҳ������ͼ��
    % the second term:
%     [Fx,Fy]=gradient(x)������FxΪ��ˮƽ�����ϵ��ݶȣ�FyΪ�䴹ֱ�����ϵ��ݶȣ�Fx�ĵ�һ��Ԫ��Ϊԭ����ڶ������һ��Ԫ��֮�Fx�ĵڶ���Ԫ��Ϊԭ������������һ��Ԫ��֮�����2
    [dfm_y,dfm_x,dfm_z]=Gradient(f_m);  % <-��ʱͼ��Ϊ������ĳһ���ͼ��
������    �����������ԭͼ����    
Ҫ�����������ԭͼ��λ�ú�֮ǰ���ص��λ��  
    % �������ǶԵ�i�����Ƶ��󵼻��ǶԿ��Ƶ�ĵ�i��������?
    % the third term:
    % p��x����ʾ����λ��   ...һ���������صĿռ��ڲ���dk
    dk=prod(max((1-abs(p-x)./K),0));
    
����   
    neighborhood_trans(ind,image,K); % K1,K2,K3��С����...�еĵ��Ӧ��ԭͼ�б�� 
    
    ki_up=Metric+(rho/2).*norm(D(p,N,f_m)-(Z-U),'fro');
    
    ed_up=dMetric.*dfm_x.*dk+dMetric.*dfm_y.*dk+dMetric.*dfm_z.*dk;
    
    dki_up=ed_up+;
end


%  ���º��������Ϊͼ��,����ҲΪͼ��
function [df_y,df_x,df_z]=Gradient(f)
% ����fΪͼ��,����ҲΪͼ��
    if numel(size(f))<3
        df_z=0;
        [df_y,df_x]=gradient(f);
    else
        [df_y,df_x,df_z]=gradient(f);
    end
end

function [f, df] =SSD(f_f,f_m,v)    
% ����ͼ��Ϊԭͼ��С,����ҲΪͼ�� 
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
    
%     if numel(size(f_f))>2   % ��ʵ����ֱ��������sum(3D)...����Ҳ����2D�µ���Ͳ���
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

