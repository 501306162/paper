function Idown = down_sample(I,N) 
[row,col] = size(I); 
drow = round(row/N); 
dcol = round(col/N); 
Idown = zeros(drow,dcol); 
p =1; 
q =1; 
for i = 1:N:row 
    for j = 1:N:col 
         Idown(p,q) = I(i,j); 
         q = q+1; 
    end 
    q =1; 
    p = p+1; 
end 
end 