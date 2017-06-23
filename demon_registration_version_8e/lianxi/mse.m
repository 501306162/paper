function [mse_result] = mse(f1, f2)
%���У����ֵԽ�󣬴�����ԭͼ���Խ��
% MSE = E( (img��Eimg)^2 ) 
%     = SUM((img-Eimg)^2)/(M*N);
%����f1��f2���������
e = double(f1) - double(f2);
[m, n] = size(e);
mse_results = sqrt(sum(e.^2)/(m*n));       %�õ�����һ�����������ս��Ҫʹ��ƽ��ֵ
mse_result = mean2(mse_results);
end
