function Dd = D_conju_new(dr,dc)  %2-D
%      con_distance(1:N)=1; % ���Ƶ���
     % ������ʽ�洢������ֻ�������ڵ㲻Ϊ0,��+1/delta or -1/delta;
     % ����delta ����L*Lά
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


% ����Ϊ����L��--->��Ϊͼ��ά��
% ����ֵ��һ��Ϊ��һά��,�ڶ��ж�Ӧ�ڶ�ά��,��ΪL��
%        dr=d(1,:);
%        dc=d(2,:);
%        dr=reshape(dr,m,n);
%        dc=reshape(dc,m,n);
      [d_rc,d_rr]=gradient(dr);
      [d_rrc,~]=gradient(d_rr);
      [~,d_rcr]=gradient(d_rc);
      
      [d_cc,d_cr]=gradient(dc);      
      [d_crc,~]=gradient(d_cr);
      [~,d_ccr]=gradient(d_cc);

      theta_1=d_rrc+d_rcr;
      theta_2=d_crc+d_ccr;
    
     Dd=[theta_1(:)';theta_2(:)']; 
      
        
end