function Dd = D(d,N,f_m)  %2-D
     con_distance(1:N)=1; % 控制点间距
     % 按矩阵方式存储。。。只有其相邻点不为0,即+1/delta or -1/delta;
     % 所以delta 会是L*L维
     % forward difference :delta_i(dj[l]) in 2-D
     Dd=zeros(N*N,numel(f_m));
     k=1;
     for j=1:N
     % displacement in j-th dimension :d(:,j);
        dis=reshape(d(:,j),size(f_m));
        [md,nd]=size(dis);
     % forward difference of j-th dispalcement,in i-th direction   
     % no.1 direction
        r_1=dis(2:md,:);
        r_0=dis(1:md-1,:);
        d_r=(r_1-r_0)/con_distance(1);
        d_r(md,:)=0;
     % no.2 direction   
        c_1=dis(:,2:nd);
        c_0=dis(:,1:nd-1);
        d_c=(c_1-c_0)/con_distance(2);
        d_c(:,nd)=0;
        Dd(k:k+1,:)=[d_r(:)';d_c(:)'];
        k=k+2;
     end
%      Dd=Dd(:);
end