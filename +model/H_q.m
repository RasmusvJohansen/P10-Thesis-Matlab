function H_q = H_q(q,param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Mq = cell(1,param.n);
for i = 1:param.n 
    Mq{i} = param.model.M{i}*q(i);
end
H_q = blkdiag(Mq{:});
end