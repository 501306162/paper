function Dd = D_conju(temp,m,n)  %2-D

% ����Ϊ����L��--->��Ϊͼ��ά��
% ����ֵ��һ��Ϊ��һά��,�ڶ��ж�Ӧ�ڶ�ά��,��ΪL��
       num_control=m*n;
       
       r_dr=temp(1,1:num_control);
       r_dc=temp(2,1:num_control);
       c_dr=temp(3,1:num_control);
       c_dc=temp(4,1:num_control);
       
       r_dr=reshape(r_dr,m,n);
       r_dc=reshape(r_dc,m,n);
       c_dr=reshape(c_dr,m,n);
       c_dc=reshape(c_dc,m,n);

%       [r_drc,r_drr]=gradient(r_dr);
%       [r_dcc,r_dcr]=gradient(r_dc);
%       [c_drc,c_drr]=gradient(c_dr);
%       [c_dcc,c_dcr]=gradient(c_dc);

      [r_drc,~]=gradient(r_dr);
      [~,r_dcr]=gradient(r_dc);
      [c_drc,~]=gradient(c_dr);
      [~,c_dcr]=gradient(c_dc);
      
      theta_1=r_drc+r_dcr;
      theta_2=c_drc+c_dcr;
    
     Dd=[theta_1(:)';theta_2(:)']; 
      
        
end