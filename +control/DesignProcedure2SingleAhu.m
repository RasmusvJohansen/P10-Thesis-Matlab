function Ks = DesignProcedure2SingleAhu(param,info)
%DESIGNPROCEDURE2 Summary of this function goes here
% Detailed explanation goes here

Q = cell(1,1);
y = cell(1,1);
R = cell(1,1);
P = cell(1,1);
for i=1:1
    Q{i} = sdpvar(3);
    y{i} = sdpvar(1,3,'full');
    R{i} = sdpvar(3);
    P{i} = sdpvar(3);
    
end
Cali_Q = blkdiag(Q{:});
Cali_Y = blkdiag(y{:});
Cali_phi = sdpvar(1);   
%Lyapuanov 
Lyap = param.model.A(1:3,1:3)*Cali_Q +Cali_Q*param.model.A(1:3,1:3).' + param.model.B_Bar(1:3,1)*Cali_Y + Cali_Y.'*param.model.B_Bar(1:3,1).';
constraints = [Lyap <= 0]

%Robust stability 
LPV = [param.model.A(1:3,1:3)*Cali_Q + Cali_Q*param.model.A(1:3,1:3).' + Cali_phi*param.model.B_Bar(1:3,1)*param.model.B_Bar(1:3,1).' Cali_Y.';
        Cali_Y -Cali_phi*eye(1)];

constraints = [constraints, LPV <= 0]

constraints = [constraints, Cali_phi >= 0, Cali_Q >= 0];

options = sdpsettings('verbose',0,'solver','mosek');
sol = optimize(constraints,[Cali_phi],options)
    if info
        check(constraints);
    end
    value(Cali_phi)
    K = value(Cali_Y)/value(Cali_Q)
    tau = 1/value(Cali_phi)
    Ks = zeros(1,3,1);
    for i=1:1
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end

end

