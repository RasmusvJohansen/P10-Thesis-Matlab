function Ks = DesignProcedure1(param)
%DESIGNPROCEDURE1 Summary of this function goes here
%   Returns a control gain based on design procedure 1 (Repport)
%k = zeros(1,3,param.n);

yalmip('clear')

Q = cell(1,param.n);
y = cell(1,param.n);
R = cell(1,param.n);

for i=1:param.n
    Q{i} = sdpvar(3);
    y{i} = sdpvar(1,3,'full');

end
cali_Q = blkdiag(Q{:})
Y = blkdiag(y{:})
%Lyapuanov
constraints = [param.model.A*cali_Q+cali_Q*param.model.A.'+param.model.B_Bar*Y+Y.'*param.model.B_Bar.' <=0]
constraints = [constraints, cali_Q >= 0]

%Robust stability
%constraints = []
options = sdpsettings('verbose',0,'solver','mosek');
sol = optimize(constraints,[],options);
check(constraints)

K = double(Y)/double(cali_Q);
eig(param.model.A+param.model.B_Bar*K)
Ks = zeros(1,3,param.n);
for i=1:param.n
    Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
end
Ks
end

