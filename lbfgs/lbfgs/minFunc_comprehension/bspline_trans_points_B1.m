function [Tlocal]=bspline_trans_points_B1(O_trans,Spacing,X,check)
% If Spacing and points are not an integer, we can not use fast look up tables,
% but have to calculate the bspline coefficients for every point
% X is [x(:),y(:)]=ndgrid(0:size(Iin,1)-1,0:size(Iin,2)-1);
%  the output Tlocal is the image dimensional deformation 
if(nargin<4)
    check = any(mod(Spacing,1)>0)|any(mod(X(:),1)>0);
end
switch size(X,2)
    case 2,
        if(check)
            [Tlocal]=bspline_transform_slow_2d(O_trans,Spacing,X);
        else
            [Tlocal]=bspline_transform_fast_2d(O_trans,Spacing,X);
        end
    case 3,
        [Tlocal]=bspline_transform_slow_3d(O_trans,Spacing,X);
end
end
function [Tlocal]=bspline_transform_fast_2d(O_trans,Spacing,X)
% Make row vectors of input coordinates
x2=X(:,1); y2=X(:,2); % 像素点坐标0~size(I)-1

% Make polynomial look up tables
Bu=zeros(2,Spacing(1));
Bv=zeros(2,Spacing(2));

x=0:Spacing(1)-1;
u=(x/Spacing(1))-floor(x/Spacing(1));
Bu(0*Spacing(1)+x+1) = (1-u);
Bu(1*Spacing(1)+x+1) = u;


y=0:Spacing(2)-1;
v=(y/Spacing(2))-floor(y/Spacing(2));
Bv(0*Spacing(2)+y+1) = (1-v);
Bv(1*Spacing(2)+y+1) = v;


% Calculate the indexes need to loop up the B-spline values.
u_index=mod(x2,Spacing(1));  %图像像素点在每个块中的标号
v_index=mod(y2,Spacing(2));

i=floor(x2/Spacing(1)); % (first row outside image against boundary artefacts)
j=floor(y2/Spacing(2));  % 块的标号
%   k0 
% ...__x0__| __|__...
% ...__|_0_|_1_|_2_...
% ...__| __| __|__...
%      |   |   | 
% This part calculates the coordinates of the pixel
% which will be transformed to the current x,y pixel.
% 
% % x2,y2是像素点坐标...IndexO1是像素点对应控制点标号,即坐标
% x_dif=u_index;
% y_dif=v_index;
% ddk_x=1-abs(x_dif)./Spacing(1);
% ddk_y=1-abs(y_dif)./Spacing(2);
% ddk=max(ddk_x,0).*max(ddk_y,0);

Ox=O_trans(:,:,1); % 网格点的位移...
Oy=O_trans(:,:,2);

Tlocalx=0; Tlocaly=0;

a=zeros(size(X,1),2);
b=zeros(size(X,1),2);

IndexO1l=zeros(size(X,1),2);  % the image's x index 
IndexO2l=zeros(size(X,1),2);
Check_bound1=false(size(X,1),2);
Check_bound2=false(size(X,1),2);

for r=0:1,
    a(:,r+1)=Bu(r*Spacing(1)+u_index(:)+1);
    b(:,r+1)=Bv(r*Spacing(2)+v_index(:)+1);

    IndexO1l(:,r+1)=(i+r);
    IndexO2l(:,r+1)=(j+r);   % 像素块的标号
    Check_bound1(:,r+1)=(IndexO1l(:,r+1)<0)|(IndexO1l(:,r+1)>(size(O_trans,1)-1));  % 判断每个点的控制点坐标是否超出网格坐标
    Check_bound2(:,r+1)=(IndexO2l(:,r+1)<0)|(IndexO2l(:,r+1)>(size(O_trans,2)-1));
end

    
for l=0:1,
    for m=0:1,
        IndexO1= IndexO1l(:,l+1);
        IndexO2= IndexO2l(:,m+1);
        Check_bound=Check_bound1(:,l+1)|Check_bound2(:,m+1);
        IndexO1(Check_bound)=1;   % check为1时,即是边界,变值为1
        IndexO2(Check_bound)=1;
        Check_bound_inv=double(~Check_bound);

        ab=a(:,l+1).*b(:,m+1);
        
%         c=Ox(IndexO1(:)+IndexO2(:)*size(Ox,1)+1);  % i,j 号网格...像素点四周网格点的x坐标
        c=Ox(IndexO1(:)+IndexO2(:)*size(Ox,1)+1);  % i,j 号网格块...像素点四周网格点的x方向的位移
      
%         Tlocalx=Tlocalx+Check_bound_inv(:).*ab.*c; % 获得每个像素点新的x坐标...即经过B样条处理的坐标
        Tlocalx=Tlocalx+Check_bound_inv(:).*ab.*c; % 获得每个像素点新的x方向的位移...即经过B样条处理的位移量

        c=Oy(IndexO1(:)+IndexO2(:)*size(Oy,1)+1); % 同理
        Tlocaly=Tlocaly+Check_bound_inv(:).*ab.*c;

    end
end
Tlocal(:,1)=Tlocalx(:);
Tlocal(:,2)=Tlocaly(:);

end

function Tlocal=bspline_transform_slow_2d(O_trans,Spacing,X)
% --->一阶B样条处理
% Make row vectors of input coordinates
x2=X(:,1); y2=X(:,2); % 图像坐标

% This code calculates for every coordinate in X, the indices of all
% b-spline knots which have influence on the transformation value of
% this point
[m,l]=ndgrid(0:1,0:1); m=m(:)'; l=l(:)';

ixs=floor(x2/Spacing(1));  % 图像点所属的块标号 。。。1 列向量 ...spacing 取4
                            % x2:  0 1 2 3 ,4 5 6 ...     ->  0,0,0,0 * 1,1,1,1...
iys=floor(y2/Spacing(2));   % y2:  0 0 0 0 ,0 0 0 ...     ->  0,0,0,0 * 0,0,0,0... 竖着读标号矩阵

 % 像素点周围4*4范围的控制点在控制矩阵中的坐标                           
ix=repmat(ixs,[1 4])+repmat(m,[length(x2) 1]); ix=ix(:);  
iy=repmat(iys,[1 4])+repmat(l,[length(y2) 1]); iy=iy(:);   

% 以 （x2，y2）=（0，0）为例： （0到第一个SPACING点范围内的值都和下方一致）。。。即，将spacing范围内所有点限制到 0号 控制点上
% ix=0,1,2,3...0,1,2,3...0,1,2,3...0,1,2,3...        
% iy=0,0,0,0...1,1,1,1...2,2,2,2...3,3,3,3...
% ... 0 __ d__ h__l___
% ... a |__e|__i|__m|...      以00为原点的，a->n,p顺序的控制点
% ... b |__f|__j|__n|...
% ... c |__g|__k|__p|...
%     1 |   |   |   |
%以 （x2，y2）=（4，0）为例： 
% ix=1,2,3,4...1,2,3,4...1,2,3,4...1,2,3,4...        
% iy=0,0,0,0...1,1,1,1...2,2,2,2...3,3,3,3...

% Size of the b-spline grid
s=size(O_trans);   % 网格大小
  
% Points outside the bspline grid are set to the upper corner
Check_bound=(ix<0)|(ix>(s(1)-1))|(iy<0)|(iy>(s(2)-1));   % 超出位置设为1 其余0
ix(Check_bound)=1; iy(Check_bound)=1;
Check_bound_inv=double(~Check_bound);       % 超出为0 其余1

% Look up the b-spline knot values in neighborhood of the points in (x2,y2)
Cx=O_trans(ix+iy*s(1)+1).*Check_bound_inv;  % 获得未超范围像素对应的控制点的行位移值
Cx=reshape(Cx,[length(x2) 4]);

Cy=O_trans(ix+iy*s(1)+s(1)*s(2)+1).*Check_bound_inv;  % 列位移值
Cy=reshape(Cy,[length(x2) 4]); 

% Calculate the b-spline interpolation constants u,v in the center cell
% range between 0 and 1
v  = (x2-ixs*Spacing(1))/Spacing(1);  % x2-ixs*Spacing(1) 是像素点距离起始网格控制点的距离...除以分母,表示占该网格的比例
u  = (y2-iys*Spacing(2))/Spacing(2);

% Get the b-spline coefficients in a matrix W, which contains
% the influence of all knots on the points in (x2,y2)
W=bspline_coefficients_B1(v,u);

% Calculate the transformation of the points in (x2,y2) by the b-spline grid
Tlocal(:,1)=sum(W.*Cx,2);
Tlocal(:,2)=sum(W.*Cy,2);
end

function Tlocal=bspline_transform_slow_3d(O_trans,Spacing,X)
% Make row vectors of input coordinates
x2=X(:,1); y2=X(:,2); z2=X(:,3);

% This code calculates for every coordinate in X, the indices of all
% b-spline knots which have influence on the transformation value of
% this point
[m,l,k]=ndgrid(0:1,0:1,0:1); m=m(:)'; l=l(:)'; k=k(:)';
ixs=floor(x2/Spacing(1));
iys=floor(y2/Spacing(2));
izs=floor(z2/Spacing(3));
ix=repmat(floor(ixs),[1 8])+repmat(m,[length(x2) 1]); ix=ix(:);
iy=repmat(floor(iys),[1 8])+repmat(l,[length(y2) 1]); iy=iy(:);
iz=repmat(floor(izs),[1 8])+repmat(k,[length(z2) 1]); iz=iz(:);

% Size of the b-spline grid
s=size(O_trans);

% Points outside the bspline grid are set to the upper corner
Check_bound=(ix<0)|(ix>(s(1)-1))|(iy<0)|(iy>(s(2)-1))|(iz<0)|(iz>(s(3)-1));
ix(Check_bound)=1; iy(Check_bound)=1; iz(Check_bound)=1;
Check_bound_inv=double(~Check_bound);

% Look up the b-spline knot values in neighborhood of the points in (x2,y2)
Cx=O_trans(ix+iy*s(1) +iz*s(1)*s(2) +                    1).*Check_bound_inv;
Cy=O_trans(ix+iy*s(1) +iz*s(1)*s(2) + s(1)*s(2)*s(3)   + 1).*Check_bound_inv;
Cz=O_trans(ix+iy*s(1) +iz*s(1)*s(2) + s(1)*s(2)*s(3)*2 + 1).*Check_bound_inv;
Cx=reshape(Cx,[length(x2) 8]);
Cy=reshape(Cy,[length(x2) 8]);
Cz=reshape(Cz,[length(x2) 8]);

% Calculate the b-spline interpolation constants u,v in the center cell
% range between 0 and 1
v  = (x2-ixs*Spacing(1))/Spacing(1);
u  = (y2-iys*Spacing(2))/Spacing(2);
w  = (z2-izs*Spacing(3))/Spacing(3);

% Get the b-spline coefficients in a matrix W, which contains
% the influence of all knots on the points in (x2,y2)
W=bspline_coefficients_B1(v,u,w);

% Calculate the transformation of the points in (x2,y2) by the b-spline grid
Tlocal(:,1)=sum(W.*Cx,2);
Tlocal(:,2)=sum(W.*Cy,2);
Tlocal(:,3)=sum(W.*Cz,2);
end