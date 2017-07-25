function Dd = D_conju(temp,m,n,q)  %2-D
  if(nargin<4)  %2d
% 输入为两行L列--->变为图像维度
% 返回值第一行为第一维度,第二行对应第二维度,都为L列
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
     
  else
       num_control=m*n*q;
       
       r_dr=temp(1,1:num_control);
       r_dc=temp(2,1:num_control);
       r_dz=temp(3,1:num_control);
       c_dr=temp(4,1:num_control);
       c_dc=temp(5,1:num_control);
       c_dz=temp(6,1:num_control);   
       z_dr=temp(7,1:num_control);
       z_dc=temp(8,1:num_control);
       z_dz=temp(9,1:num_control);  
       
       r_dr=reshape(r_dr,m,n,q);
       r_dc=reshape(r_dc,m,n,q);
       r_dz=reshape(r_dz,m,n,q);
       
       c_dr=reshape(c_dr,m,n,q);
       c_dc=reshape(c_dc,m,n,q);
       c_dz=reshape(c_dz,m,n,q);
       
       z_dr=reshape(z_dr,m,n,q);
       z_dc=reshape(z_dc,m,n,q);
       z_dz=reshape(z_dz,m,n,q);

      [r_drc,r_drr,r_drz]=gradient(r_dr);
      [r_dcc,r_dcr,r_dcz]=gradient(r_dc);     
      [r_dzc,r_dzr,r_dzz]=gradient(r_dz);
      
      [c_drc,c_drr,c_drz]=gradient(c_dr);
      [c_dcc,c_dcr,c_dcz]=gradient(c_dc);
      [c_dzc,c_dzr,c_dzz]=gradient(c_dz);
      
      [z_drc,z_drr,z_drz]=gradient(z_dr);
      [z_dcc,z_dcr,z_dcz]=gradient(z_dc);
      [z_dzc,z_dzr,z_dzz]=gradient(z_dz);
      
      theta_1=r_drc+r_dcr   +r_drz+r_dzr   +r_dcz+r_dzc;
      theta_2=c_drc+c_dcr   +c_drz+c_dzr   +c_dcz+c_dzc;
      theta_3=z_drc+z_dcr   +z_drz+z_dzr   +z_dcz+z_dzc;  
      
     Dd=[theta_1(:)';theta_2(:)';theta_3(:)']; 
  end
        
end