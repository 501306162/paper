function mse=myMSE( X,Y )

X=double(X);  
Y=double(Y); 
MSE=0;
[r,c]=size(X);
for i=1:r
    for j=1:c
        MSE=MSE+(X(i,j)-Y(i,j)).*(X(i,j)-Y(i,j));
    end
end
mse=MSE/(r*c);
end

