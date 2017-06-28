function Iout=movepixels_2d(Iin_o,Tx,Ty)
% This function movepixels, will translate the pixels of an image
%  according to x and y translation images (bilinear interpolated). 
% 
%  Iout = movepixels_2d_double(I,Tx,Ty,mode);
%
% Inputs;
%   Tx, Ty: The transformation images, describing the
%             (backwards) translation of every pixel in x and y direction.
%   mode:   linear interpolation and outside pixels set to nearest pixel
%
% Outputs,
%   Iout : The transformed image
%
% Function is written by D.Kroon University of Twente (February 2009)
Iin=double(Iin_o(:,:));
Tx=double(Tx);Ty=double(Ty);
% Make all x,y indices
[x,y]=ndgrid(0:size(Iin,1)-1,0:size(Iin,2)-1);    % 0--->
                                                  % |
                                                  % \/
% Calculate the Transformed coordinates
Tlocalx = x+Tx;
Tlocaly = y+Ty;

% All the neighborh pixels involved in linear interpolation.
xBas0=floor(Tlocalx); 
yBas0=floor(Tlocaly);
xBas1=xBas0+1;           
yBas1=yBas0+1;

% Linear interpolation constants (percentages)
xCom=Tlocalx-xBas0; 
yCom=Tlocaly-yBas0;
perc0=(1-xCom).*(1-yCom);
perc1=(1-xCom).*yCom;
perc2=xCom.*(1-yCom);
perc3=xCom.*yCom;

% limit indexes to boundaries
check_xBas0=(xBas0<0)|(xBas0>(size(Iin,1)-1));
check_yBas0=(yBas0<0)|(yBas0>(size(Iin,2)-1));
xBas0(check_xBas0)=0; 
yBas0(check_yBas0)=0; 
check_xBas1=(xBas1<0)|(xBas1>(size(Iin,1)-1));
check_yBas1=(yBas1<0)|(yBas1>(size(Iin,2)-1));
xBas1(check_xBas1)=0; 
yBas1(check_yBas1)=0; 


Iout=zeros(size(Iin));
  
Iin_one=Iin(:,:);
% Get the intensities
intensity_xyz0=Iin_one(1+xBas0+yBas0*size(Iin,1));
intensity_xyz1=Iin_one(1+xBas0+yBas1*size(Iin,1)); 
intensity_xyz2=Iin_one(1+xBas1+yBas0*size(Iin,1));
intensity_xyz3=Iin_one(1+xBas1+yBas1*size(Iin,1));

Iout_one=intensity_xyz0.*perc0+intensity_xyz1.*perc1+intensity_xyz2.*perc2+intensity_xyz3.*perc3;
Iout(:,:)=reshape(Iout_one, [size(Iin,1) size(Iin,2)]);

if(~isa(Iin_o,'double')&&~isa(Iin_o,'single'))
    if(isa(Iin_o,'uint8')), Iout=uint8(Iout); end
    if(isa(Iin_o,'uint16')), Iout=uint16(Iout); end
    if(isa(Iin_o,'uint32')), Iout=uint32(Iout); end
    if(isa(Iin_o,'int8')),   Iout=int8(Iout); end
    if(isa(Iin_o,'int16')), Iout=int16(Iout); end
    if(isa(Iin_o,'int32')), Iout=int32(Iout); end
end
end
    
    


 
    