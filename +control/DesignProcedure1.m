function Ks = DesignProcedure1(param, info)
%DESIGNPROCEDURE1 Summary of this function goes here
%   Returns a control gain based on design procedure 1 (Repport)
%k = zeros(1,3,param.n);
    
    yalmip('clear')
    
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
    cali_Q = blkdiag(Q{:});
    Y = blkdiag(y{:});
    %cali_R = sdpvar(12,12);
    cali_R = blkdiag(R{:});
    cali_P = blkdiag(P{:});
    %Lyapuanov
    constraints = [param.model.A*cali_Q+cali_Q*param.model.A.'+param.model.B_Bar *Y + Y.'* param.model.B_Bar.' <=0];
    constraints = [constraints,  cali_Q >= 0];
    
    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,[],options)
    if info
        check(constraints);
    end
    
    K = value(Y)/value(cali_Q);
    Ks = zeros(1,3,param.n);
    for i=1:param.n
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end
    
    clear constraints
    tau = sdpvar(1);
    %Robust stability
    LPV = [cali_R * param.model.A + param.model.A.' * cali_R + tau*(K.'*K), cali_R*param.model.B; param.model.B.'*cali_R, -tau*eye(4)];
    
    constraints = LPV <= 0;
    constraints = [constraints, cali_R >= 0, tau >=0];
    
    
    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,[],options);
    if info
        check(constraints);
    end
    
    disp("LPV")
    ttau =1;
    LPV = [value(cali_R) * param.model.A + param.model.A.' * value(cali_R) + ttau*(K.'*K), value(cali_R)*param.model.B; param.model.B.'*value(cali_R), -ttau*eye(4)];
    disp("eig LPV")
    eig(LPV)
    fprintf("Evaluation of design procedure 1 \n")
    if(eig(value(cali_Q))<=0)
        fprintf("cali_Q is not positive definite!!!\n")
    else
        fprintf("cali_Q is positive definite\n")
    end
    eig(value(cali_Q)).'

    if(eig(value(cali_R))<=0)
        fprintf("cali_R is not positive definite!!!\n")
    else
        fprintf("cali_R is positive definite\n")
    end
    eig(value(cali_R)).'
    
    fprintf("Tau is found to be:\n")
    value(tau)
    
    if(eig(param.model.A+param.model.B_Bar*K) > 0)
        fprintf("One or more eigenvalues in the decoupled system is positive!!! \n")
    else
        fprintf("Eigenvalues of the decoupled systems are all negative \n")
    end
    
    eig(param.model.A+param.model.B_Bar*K).'
    
    if(eig(param.model.A+param.model.B*K) > 0)
        fprintf("One or more eigenvalues in the coupled system is positive!!! \n")
    else
        fprintf("Eigenvalues of the coupled systems are all negative \n")
    end
    eig(param.model.A+param.model.B*K).'
    
    
    disp("eigenvalues of A")
    eig(param.model.A)
    disp("eigenvlaues of Bk")
    eig(param.model.B*K)
   
    %Lyapuanov på specefict gain
    clear constraints
    eps1 = 1; 
    eps2 = 0; 
    eps3 = 0;
    eps4 = 0;
    epsilon =[eps1 0 0 0; 0 eps2 0 0; 0 0 eps3 0; 0 0 0 eps4]

    constraints = [param.model.A.'*cali_P + K.'*epsilon.'*param.model.B.'*cali_P + cali_P*param.model.A + cali_P*param.model.B*epsilon*K <=0];
    constraints = [constraints, cali_P >= 0];

    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,[],options);
    if info
        check(constraints);
    end
    eig(value(cali_P))
    eig(param.model.A.'*value(cali_P) + K.'*epsilon.'*param.model.B.'*value(cali_P) + value(cali_P)*param.model.A + value(cali_P)*param.model.B*epsilon*K)

  
end

