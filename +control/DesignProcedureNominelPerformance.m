function [Ks] = DesignProcedureNominelPerformance(param,info)
%DESIGNPROCEDURENOMINELPERFORMANCE Summary of this function goes here
%   Detailed explanation goes here


Q = cell(1,param.n);
y = cell(1,param.n);
R = cell(1,param.n);
P = cell(1,param.n);
for i=1:param.n
    Q{i} = sdpvar(3);
    y{i} = sdpvar(1,3,'full');
    R{i} = sdpvar(3);
    P{i} = sdpvar(3);
    
end
Cali_Q = blkdiag(Q{:});
Cali_R = blkdiag(R{:});
Cali_Y = blkdiag(y{:});
Cali_phi = sdpvar(1);   
gamma = sdpvar(1);
B1 = -blkdiag([0;0;1],[0;0;1],[0;0;1],[0;0;1]);
C1 = blkdiag([0 0 1],[0 0 1],[0 0 1],[0 0 1]);
alpha = 0.0001;
gamma2 = 1


%Nominel performance 
NP = [Cali_Q*param.model.A.' + param.model.A*Cali_Q + Cali_Y.'*param.model.B_Bar.'+param.model.B_Bar*Cali_Y + C1.'*C1, B1;
      B1.', -eye(param.n)*gamma];

constraints = [NP <= 0];

%limit on integral state 
% x_bound = [Cali_Q/(eye(12)*gamma2); Cali_Q*[0 0 1 0 0 1 0 0 1 0 0 1].'; 
%            [0 0 1 0 0 1 0 0 1 0 0 1].'*Cali_Q 1]
% constraints = [constraints, x_bound >= 0]



%Robust stability 
LPV = [param.model.A*Cali_Q + Cali_Q*param.model.A.' + Cali_phi*param.model.B*param.model.B.' Cali_Y.';
        Cali_Y -Cali_phi*eye(param.n)];
constraints = [constraints, LPV <= 0]

constraints = [constraints, Cali_phi >= 0, 1>= Cali_phi, Cali_Q >= 0];

options = sdpsettings('verbose',0,'solver','mosek');
sol = optimize(constraints,[gamma, -Cali_phi],options)

disp("Eig of Q")
eig(value(Cali_Q))
disp("K")
K = value(Cali_Y)/value(Cali_Q)
disp("tau")
tau = 1/value(Cali_phi)
disp("gamma")
gam = sqrt(value(gamma))

disp("Eig of NP")
eig(value(NP))

disp("eig of LPV")
% eig(value(LPV))
Ks = zeros(1,3,param.n);
    for i=1:param.n
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end
end

