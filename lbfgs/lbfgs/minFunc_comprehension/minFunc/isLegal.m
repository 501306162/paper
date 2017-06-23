function [legal] = isLegal(v)
legal = sum(any(imag(v(:))))==0 & sum(isnan(v(:)))==0 & sum(isinf(v(:)))==0;

% B=any(A)，如果A是向量，如果向量里有非0的数，则返回1（true），否则返回0
%           如果A是矩阵，则把矩阵的列当做向量来处理，函数返回每个列向量的逻辑值

% 该函数即表示:
% 当输入向量中不含0,不含非参数,不含无穷大值,为合法,返回1