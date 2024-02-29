function Ks = DesignProcedureDScalingDiag(param,info)
%Calculates the feedback gain for diagoanl input uncertenties
%   Detailed explanation goes here
 
     yalmip('clear')
    
   Q = cell(1,param.n);
    y = cell(1,param.n);
    R = cell(1,param.n);
    P = cell(1,param.n);
    l = cell(1,param.n)
    for i=1:param.n
        Q{i} = sdpvar(3);
        y{i} = sdpvar(1,3,'full');
        R{i} = sdpvar(3);
        P{i} = sdpvar(3);
        l{i} = sdpvar(1);
    end
    cali_Q = blkdiag(Q{:});
    Y = blkdiag(y{:});
    %cali_R = sdpvar(12,12);
    cali_R = blkdiag(R{:});
    cali_P = blkdiag(P{:});
    L = blkdiag(l{:});
    %Lyapuanov
    constraints = [param.model.A*cali_P+cali_P*param.model.A.'+param.model.B_Bar *Y + Y.'* param.model.B_Bar.' <=0];
    constraints = [constraints,  cali_P >= 0];
    
    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,[],options)
    if info
        check(constraints);
    end
    
    K = value(Y)/value(cali_P)
    Ks = zeros(1,3,param.n);
    for i=1:param.n
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end
    C = [eye(3); eye(3); eye(3);eye(3)].'
    clear constraints
    tau = sdpvar(1);
    %Robust stability
    
    LPV = [param.model.A.'*cali_R+ cali_R*param.model.A cali_R*param.model.B_Bar*L eye(12)*L;
           s*param.model.B_Bar'*cali_R -L zeros(4,4);
           S*eye(12) zeros(4,4) L ];

    constraints = LPV <= 0;
    constraints = [constraints, cali_R >= 0, tau >=0, L>= 0];
    
    
    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,[],options);
    if info
        check(constraints);
    end
    tau = value(tau)
    value(L)
    disp("eig of L")
    eig(value(L))
end