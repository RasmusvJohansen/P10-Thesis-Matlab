function H_xi = H_xi(xi,param)
%Calculate H_xi from the state vector xi
%   Detailed explanation goes here'

%calculates M1*Z1 ... Mn*Zn and store in a cell to write it into a block
%diagonal matrice
MZ = cell(1,param.n);
for i = 1:param.n 
    MZ{i} = param.model.M{i}*xi(3*(i-1)+1:i*3);% fetch 1:3, then 4:6, ... variables of xi
end
H_xi = blkdiag(MZ{:});
end