function [mse_result] = mse(f1, f2)
%其中，结果值越大，代表与原图相差越大
% MSE = E( (img－Eimg)^2 ) 
%     = SUM((img-Eimg)^2)/(M*N);
%计算f1和f2均方根误差
e = double(f1) - double(f2);
[m, n] = size(e);
mse_results = sqrt(sum(e.^2)/(m*n));       %得到的是一组结果集，最终结果要使用平均值
mse_result = mean2(mse_results);
end
