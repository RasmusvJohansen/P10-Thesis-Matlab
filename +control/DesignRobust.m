function Ks = DesignRobust(param,info)
%DESIGNROBUST Summary of this function goes here
%   Detailed explanation goes here

   yalmip('clear')
    
   Q = cell(1,param.n);
    y = cell(1,param.n);
    R = cell(1,param.n);
    P = cell(1,param.n);
    Kg = cell(1,param.n);
    for i=1:param.n
        Q{i} = sdpvar(3);
        y{i} = sdpvar(1,3,'full');
        R{i} = sdpvar(3);
        P{i} = sdpvar(3);
        Kg{i} = sdpvar(3,3);
    end
    cali_Q = blkdiag(Q{:});
    Y = blkdiag(y{:});
    %cali_R = sdpvar(12,12);
    cali_R = blkdiag(R{:});
    cali_P = blkdiag(P{:});
    
    tau = sdpvar(1);
    tbk = blkdiag(Kg{:})
    LPV = [cali_R * param.model.A + param.model.A.' * cali_R, cali_R*param.model.B_Bar; param.model.B_Bar.'*cali_R, eye(4)];
    constraints = LPV <= 0;
    Robustnes =  [-tbk, zeros(12,4); zeros(4,12), eye(4)]
    constraints = [constraints, Robustnes <= 0]
    
    
    constraints = [constraints, cali_R >= 0,tau >=1];
    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,[],options)
    if info
        check(constraints);
    end
    value(tbk)
    value(tau)
    k = zeros(12,12);
    val_Big_K = value(tbk)
    
    
    for i = 1:12
        k(i,i) = sqrt((val_Big_K(i,i)));
    end
    k






     K = value(Y)/value(cali_Q);
    Ks = zeros(1,3,param.n);
    for i=1:param.n
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end
   

Ks
end

