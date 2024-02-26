function Ks = DesignProcedure2(param,info)
%DESIGNPROCEDURE2 Summary of this function goes here
% Detailed explanation goes here

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
Cali_Y = blkdiag(y{:});
Cali_phi = sdpvar(1);   
%Lyapuanov 
Lyap = param.model.A*Cali_Q +Cali_Q*param.model.A.' + param.model.B_Bar*Cali_Y + Cali_Y.'*param.model.B_Bar.';
constraints = [Lyap <= 0]

%Robust stability 
LPV = [param.model.A*Cali_Q + Cali_Q*param.model.A.' + Cali_phi*param.model.B*param.model.B.' Cali_Y.';
        Cali_Y -Cali_phi*eye(param.n)];

constraints = [constraints, LPV <= 0]

constraints = [constraints, Cali_phi >= 0, Cali_Q >= 0];

options = sdpsettings('verbose',0,'solver','mosek');
sol = optimize(constraints,[],options)
    if info
        check(constraints);
    end
    
    K = value(Cali_Y)/value(Cali_Q)
    tau = 1/value(Cali_phi)
    Ks = zeros(1,3,param.n);
    for i=1:param.n
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end

end

