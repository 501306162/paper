function Iout=movepixels_3d(Iin_o,Tx,Ty,Tz)
% This function movepixels, will translate the pixels of an image
%  according to x & y & z translation images (bilinear interpolated). 
% 
%  Iout = movepixels_3d_double(I,Tx,Ty,Tz,mode);
%
% Inputs;
%   Tx, Ty, Tz: The transformation images, describing the
%             (backwards) translation of every pixel in x ,y ,z direction.
%   mode:   linear interpolation and outside pixels set to nearest pixel
%
% Outputs,
%   Iout : The transformed image
%
% Function is written by D.Kroon University of Twente (February 2009)
Iin=double(Iin_o(:,:,:));
Tx=double(Tx);Ty=double(Ty);Tz=double(Tz(:));
% Make all x,y,z indices
[x,y,z]=ndgrid(0:size(Iin,1)-1,0:size(Iin,2)-1,0:size(Iin,3)-1);

% Calculate the Transformed coordinates
Tlocalx = x+Tx;
Tlocaly = y+Ty;
Tlocalz = z+Tz;

% All the neighborh pixels involved in linear interpolation.
xBas0=floor(Tlocalx); 
yBas0=floor(Tlocaly);
zBas0=floor(Tlocalz);
xBas1=xBas0+1;           
yBas1=yBas0+1;
zBas1=zBas0+1;

% Linear interpolation constants (percentages)
xCom=Tlocalx-xBas0; 
yCom=Tlocaly-yBas0;
zCom=Tlocalz-zBas0;
perc0=(1-xCom).*(1-yCom).*(1-zCom);
perc1=(1-xCom).*yCom.*(1-zCom);
perc2=xCom.*(1-yCom).*(1-zCom);
perc3=(1-xCom).*(1-yCom).*zCom;

perc4=xCom.*yCom.*(1-zCom);
perc5=(1-xCom).*yCom.*zCom;
perc6=xCom.*(1-yCom).*zCom;
perc7=xCom.*yCom.*zCom;
% limit indexes to boundaries
check_xBas0=(xBas0<0)|(xBas0>(size(Iin,1)-1));
check_yBas0=(yBas0<0)|(yBas0>(size(Iin,2)-1));
check_zBas0=(zBas0<0)|(zBas0>(size(Iin,3)-1));
xBas0(check_xBas0)=0; 
yBas0(check_yBas0)=0; 
zBas0(check_zBas0)=0; 
check_xBas1=(xBas1<0)|(xBas1>(size(Iin,1)-1));
check_yBas1=(yBas1<0)|(yBas1>(size(Iin,2)-1));
check_zBas1=(zBas1<0)|(zBas1>(size(Iin,3)-1));
xBas1(check_xBas1)=0; 
yBas1(check_yBas1)=0; 
zBas1(check_zBas1)=0; 

Iout=zeros(size(Iin));
  
Iin_one=Iin(:,:);
% Get the intensities
intensity_xyz0=Iin_one(1+xBas0+yBas0*size(Iin,1)+zBas0*size(Iin,1)*size(Iin,2));
intensity_xyz1=Iin_one(1+xBas0+yBas1*size(Iin,1)+zBas0*size(Iin,1)*size(Iin,2)); 

intensity_xyz2=Iin_one(1+xBas1+yBas0*size(Iin,1)+zBas0*size(Iin,1)*size(Iin,2));
intensity_xyz3=Iin_one(1+xBas0+yBas0*size(Iin,1)+zBas1*size(Iin,1)*size(Iin,2));

intensity_xyz4=Iin_one(1+xBas1+yBas1*size(Iin,1)+zBas0*size(Iin,1)*size(Iin,2));
intensity_xyz5=Iin_one(1+xBas0+yBas1*size(Iin,1)+zBas1*size(Iin,1)*size(Iin,2)); 

intensity_xyz6=Iin_one(1+xBas1+yBas0*size(Iin,1)+zBas1*size(Iin,1)*size(Iin,2));
intensity_xyz7=Iin_one(1+xBas1+yBas1*size(Iin,1)+zBas1*size(Iin,1)*size(Iin,2));

Iout_one=intensity_xyz0.*perc0+intensity_xyz1.*perc1+intensity_xyz2.*perc2+intensity_xyz3.*perc3+...
         intensity_xyz4.*perc4+intensity_xyz5.*perc5+intensity_xyz6.*perc6+intensity_xyz7.*perc7;
Iout(:,:,:)=reshape(Iout_one, [size(Iin,1) size(Iin,2) size(Iin,3)]);

if(~isa(Iin_o,'double')&&~isa(Iin_o,'single'))
    if(isa(Iin_o,'uint8')), Iout=uint8(Iout); end
    if(isa(Iin_o,'uint16')), Iout=uint16(Iout); end
    if(isa(Iin_o,'uint32')), Iout=uint32(Iout); end
    if(isa(Iin_o,'int8')),   Iout=int8(Iout); end
    if(isa(Iin_o,'int16')), Iout=int16(Iout); end
    if(isa(Iin_o,'int32')), Iout=int32(Iout); end
end
end

    
    


 
    