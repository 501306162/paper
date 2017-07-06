function [ output_args ] = k_initial_2d( Spacing, sizes,f_mj)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    [k_grid]=make_init_grid(Spacing,sizes);

    % Convert all values tot type double
     k_grid=double(k_grid); 
    
    [I_tran,Trans]=bspline_transform(k_grid,double(f_mj),Spacing);
end


function k_grid=make_init_grid(Spacing,sizeI)
    % Determine grid spacing
    dx=Spacing(1); dy=Spacing(2);

    % Calculate te grid coordinates (make the grid)
    [X,Y]=ndgrid(-dx:dx:(sizeI(1)+(dx*2)),-dy:dy:(sizeI(2)+(dy*2)));
    k_grid=ones(size(X,1),size(X,2),2);
    k_grid(:,:,1)=X;
    k_grid(:,:,2)=Y;
end

function [fm,T]=bspline_transform(O,fm,Spacing,mode)
    if(sum(Spacing-floor(Spacing))>0), error('Spacing must be a integer'); end
    if(~exist('mode','var')), mode=0; end
    if(~isa(fm,'double')), fm=im2double(fm); end
    
    if(nargout > 1 )
        [fm,Tx,Ty]=bspline_transform_2d_double(double(O(:,:,1)),double(O(:,:,2)),fm,double(Spacing(1)),double(Spacing(2)),double(mode));
        T(:,:,1)=Tx; 
        T(:,:,2)=Ty;
    else
        fm=bspline_transform_2d_double(double(O(:,:,1)),double(O(:,:,2)),fm,double(Spacing(1)),double(Spacing(2)),double(mode));
        T=0;
    end

end

function [Iout,Tx,Ty]=bspline_transform_2d_double(Ox,Oy,Iin,dx,dy,mode)
    % Make all x,y indices
    [x,y]=ndgrid(0:size(Iin,1)-1,0:size(Iin,2)-1);

    % Calulate the transformation of all image coordinates by the b-spline grid
    O_trans(:,:,1)=Ox; O_trans(:,:,2)=Oy;
    Tlocal=bspline_trans_points_double_1st(O_trans,[dx dy],[x(:) y(:)],false);

    switch(mode)
        case 0
            Interpolation='bilinear';
            Boundary='replicate';
        case 1
            Interpolation='bilinear';
            Boundary='zero';
        case 2
            Interpolation='bicubic';
            Boundary='replicate';
        otherwise
            Interpolation='bicubic';
            Boundary='zero';
    end

    if(nargout>1)
        % Store transformation fields
        Tx=reshape(Tlocal(:,1),[size(Iin,1) size(Iin,2)])-x;  
        Ty=reshape(Tlocal(:,2),[size(Iin,1) size(Iin,2)])-y; 
    end
    Iout=image_interpolation(Iin,Tlocal(:,1),Tlocal(:,2),Interpolation,Boundary);
end