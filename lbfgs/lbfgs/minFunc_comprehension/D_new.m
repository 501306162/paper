function Dd = D_new(dr,dc)  %2-D
%      con_distance(1:N)=1; % 控制点间距
     % 按矩阵方式存储。。。只有其相邻点不为0,即+1/delta or -1/delta;
     % 所以delta 会是L*L维
     % forward difference :delta_i(dj[l]) in 2-D
%      Dd=zeros(N*N,numel(f_m));
%      k=1;
%      for j=1:N
%      % displacement in j-th dimension :d(:,j);
%         dis=reshape(d(:,j),size(f_m));
%         [md,nd]=size(dis);
%      % forward difference of j-th dispalcement,in i-th direction   
%      % no.1 direction
%         r_1=dis(2:md,:);
%         r_0=dis(1:md-1,:);
%         d_r=(r_1-r_0)/con_distance(1);
%         d_r(md,:)=0;
%      % no.2 direction   
%         c_1=dis(:,2:nd);
%         c_0=dis(:,1:nd-1);
%         d_c=(c_1-c_0)/con_distance(2);
%         d_c(:,nd)=0;
%         Dd(k:k+1,:)=[d_r(:)';d_c(:)'];
%         k=k+2;
%      end
%      Dd=Dd(:);


% 输入为两行L列--->变为图像维度
% 返回值第1、2行为第一维度,第3、4行对应第二维度,都为L列
%        dr=d(1,:);
%        dc=d(2,:);
%        dr=reshape(dr,m,n);
%        dc=reshape(dc,m,n);
      [dr_dc,dr_dr]=gradient(dr);
      [dc_dc,dc_dr]=gradient(dc);
      
      Dd=[dr_dr(:)';dr_dc(:)';dc_dr(:)';dc_dc(:)'];
end