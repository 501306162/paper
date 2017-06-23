function [RegistrationParameter]=Ant_K(M,I1,I2)

%%%%%%%%%%%%%%%%%初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ant=10;   %蚂蚁数量
Times=100;   %蚂蚁移动次数
Rou=0.5;   %信息素挥发系数
P0=0.2;   %转移概率常数
Lower_1=0;   %设置搜索范围
Upper_1=1;

%}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Ant
    X(i,1)=(Lower_1+(Upper_1-Lower_1)*rand);
                    %随机设置蚂蚁的初值位置
    Tau(i)=F(X(i,1),M,I1,I2);  %计算信息量
end

for T=1:Times
    lamda=1/T;
    [Tau_Best(T),BestIndex]=max(Tau);
    for i=1:Ant
        P(T,i)=(Tau(BestIndex)-Tau(i))/Tau(BestIndex);
                                  %计算状态转移概率
    end
    for i=1:Ant
        if P(T,i)<P0  %局部搜索
            temp1=X(i,1)+(2*rand-1)*lamda;
        else      %全局搜索
             temp1=X(i,1)+(Upper_1-Lower_1)*(rand-0.5);
        end
        %%%%%%%%%%%%%%%越界处理%%%%%%%%%%%%%%%%%%%%%%
        if temp1<Lower_1
            temp1=Lower_1;
        end
        if temp1>Upper_1
            temp1=Upper_1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if F(temp1,M,I1,I2)>F(X(i,1),M,I1,I2)
           %判断蚂蚁是否移动
           X(i,1)=temp1;
        end
    end
    for i=1:Ant
        Tau(i)=(1-Rou)*Tau(i)+F(X(i,1),M,I1,I2);
          %更新信息量
    end
end

[max_value,max_index]=max(Tau);
maxX=X(max_index,1);
RegistrationParameter(1)=maxX;  %y
end

function [F]=F(k,M,I1,I2)  %目标函数

F=-mse(MM(k,M,I1,I2),I2);

end   
function [M]=MM(k,M,I1,I2)  %目标函数    

        Hsmooth=fspecial('gaussian',[60 60],10);

        Tx=zeros(size(M)); Ty=zeros(size(M));
        Idiff=M-I2;
        [My,Mx] = gradient(M);
%         Default demon force, (Thirion 1998)
        Ux = -(Idiff.*Mx)./((Mx.^2+My.^2)+k^2);
        Uy = -(Idiff.*My)./((Mx.^2+My.^2)+k^2);
  
 
        % Extended demon force. With forces from the gradients from both
% %         moving as static image. (Cachier 1999, He Wang 2005)
%         [My,Mx] = gradient(M);
%         Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+(Idiff).^2))+(Mx./((Mx.^2+My.^2)+(Idiff).^2)));
%         Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+(Idiff).^2))+(My./((Mx.^2+My.^2)+(Idiff).^2)));
% %  
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;  %判断是否是非数值参数

        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);

        % Add the new transformation field to the total transformation field.
        Tx=Tx+Uxs;
        Ty=Ty+Uys;
        M=movepixels(I1,Tx,Ty); 
end   