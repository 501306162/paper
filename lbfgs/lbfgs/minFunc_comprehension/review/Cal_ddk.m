function ddk = Cal_ddk( m,n,Spacing )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% kx,ky�ǿ��Ƶ�����
O_m=0:Spacing(1):m;
O_n=0:Spacing(2):n;
[kx,ky]=ndgrid(O_m,O_n);

% ͼ�������
[x,y]=ndgrid(0:m-1,0:n-1);

% �ҵ����Ƶ������ؾ����е�λ��
k_p=zeros(size(kx));
ddk=zeros(numel(kx),1);
for i=1:numel(k_p)
    k_p(i)=find(x==kx(i)&y==ky(i));
    % ddk=[neighbors_ind(:),ddk_1.*ddk_2]; ��øÿ��Ƶ��������еı���Լ�����ddk��ֵ
    ddk=neighborhood_trans(ind,[m n],Spacing);
end


% x_dif=u_index;
% y_dif=v_index;
% ddk_x=1-abs(x_dif)./Spacing(1);
% ddk_y=1-abs(y_dif)./Spacing(2);
% ddk=max(ddk_x,0).*max(ddk_y,0);

end

