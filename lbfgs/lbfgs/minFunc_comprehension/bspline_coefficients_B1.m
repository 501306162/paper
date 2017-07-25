function W=bspline_coefficients_B1(u,v,w)
switch(nargin)
    case 1
        W = bspline_coefficients_1d(u);
    case 2
        W = bspline_coefficients_2d(u,v);
    case 3
        W = bspline_coefficients_3d(u,v,w);
end
end    
function W=bspline_coefficients_1d(u)
W(:,1) = (1-u);
W(:,2) = u;
end

function W = bspline_coefficients_2d(u,v)
Bu=bspline_coefficients_1d(u);
Bv=bspline_coefficients_1d(v);

% Calculate influences of all neighborh b-spline knots
W = [Bu(:,1).*Bv(:,1),Bu(:,2).*Bv(:,1),...
     Bu(:,1).*Bv(:,2),Bu(:,2).*Bv(:,2)];
end

function W = bspline_coefficients_3d(u,v,w)
Bu=bspline_coefficients_1d(u);
Bv=bspline_coefficients_1d(v);
Bw=bspline_coefficients_1d(w); 

% Calculate influences of all neighborh b-spline knots
W = [Bu(:,1).*Bv(:,1).*Bw(:,1) ,Bu(:,2).*Bv(:,1).*Bw(:,1) , ...
     Bu(:,1).*Bv(:,2).*Bw(:,1) ,Bu(:,2).*Bv(:,2).*Bw(:,1) , ...
     Bu(:,1).*Bv(:,1).*Bw(:,2) ,Bu(:,2).*Bv(:,1).*Bw(:,2) , ...
     Bu(:,1).*Bv(:,2).*Bw(:,2) ,Bu(:,2).*Bv(:,2).*Bw(:,2) ];
end