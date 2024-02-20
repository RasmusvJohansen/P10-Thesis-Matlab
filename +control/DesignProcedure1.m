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
    R{i} = sdpvar(3);
end
cali_Q = blkdiag(Q{:});
Y = blkdiag(y{:});
cali_R = blkdiag(R{:});
%Lyapuanov
constraints = [param.model.A*cali_Q+cali_Q*param.model.A.'+param.model.B_Bar * inv(param.ctrl.Lambda_Bar)*Y + Y.'* inv(param.ctrl.Lambda_Bar)*param.model.B_Bar.' <=0];
constraints = [constraints, cali_Q >= 0];
constraints
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

clear constraints
tau = sdpvar(1);
%Lyapuanov
LPV = [cali_R * param.model.A + param.model.A.' * cali_R + tau*K.'*K, cali_R*param.model.B;param.model.B.'*cali_R, -tau*eye(4)]

constraints = LPV <= 0;
constraints = [constraints, cali_R >= 0, tau >=0];
constraints

%constraints = []
options = sdpsettings('verbose',0,'solver','mosek');
sol = optimize(constraints,[],options);
check(constraints)
double(tau)
eig(double(cali_R))
eig(param.model.A+param.model.B*K)
end

