function [tar]=TAR(x1,y1,x2,y2,x3,y3)

% 获得特征点i,与边缘点k  对应得到的TAR签名值

    A=[x1,y1,1; x2,y2,1;x3,y3,1;];
    tar=det(A)./2;
end
% 
% function tar=TAR(x1,y1,x2,y2,x3,y3)
%     A=[x1,y1,1; x2,y2,1;x3,y3,1;];
%     tar=det(A)./2;
% end