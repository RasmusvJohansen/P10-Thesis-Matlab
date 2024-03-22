function Ks = DesignProcedureFunky(param, info)
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
    cali_R = blkdiag(R{:});
    cali_P = blkdiag(P{:});
    %Lyapuanov
    constraints = [param.model.A*cali_Q+cali_Q*param.model.A.'+param.model.B_Bar *Y + Y.'* param.model.B_Bar.' <=0];
    constraints = [constraints, cali_Q >= 0];
    
    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,[],options);
    if info
        check(constraints);
    end
    
    K = double(Y)/double(cali_Q);
    Ks = zeros(1,3,param.n);
    for i=1:param.n
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end
    
    clear constraints
    tau1 = sdpvar(1);
    tau2 = sdpvar(1);
    %epsilon = sdpvar(1);
    te = sdpvar(1)
    %Lyapuanov
    LPV = [cali_R * param.model.A + param.model.A.' * cali_R + tau1*(K.'*K) - te*(K.'*K), cali_R*param.model.B;param.model.B.'*cali_R, -tau1*eye(4)-tau2*eye(4)];
    
    constraints = LPV <= 0;
    constraints = [constraints, cali_R >= 0, tau1 >=0,tau2 >=0,te];
    
    
    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,te,options);
    if info
        check(constraints);
    end
    
    fprintf("Evaluation of design procedure 1 \n")
    if(eig(double(cali_Q))<=0)
        fprintf("cali_Q is not positive definite!!!\n")
    else
        fprintf("cali_Q is positive definite\n")
    end
    eig(double(cali_Q)).'

    if(eig(double(cali_R))<=0)
        fprintf("cali_R is not positive definite!!!\n")
    else
        fprintf("cali_R is positive definite\n")
    end
    eig(double(cali_R)).'
    
    fprintf("Tau is found to be:\n")
    double(tau1)
    double(tau2)
    epsilon = sqrt(double(te)/double(tau2))
    
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

    %Lyapuanov på specefict gain
    clear constraints
    epsilon = 1

    constraints = [param.model.A.'*cali_P + K.'*epsilon.'*param.model.B.'*cali_P + cali_P*param.model.A + cali_P*param.model.B*epsilon*K <=0];
    constraints = [constraints, cali_P >= 0];
    
    %constraints = []
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints,[],options);
    if info
        check(constraints);
    end
    eig(param.model.A.'*double(cali_P + K.'*epsilon.'*param.model.B.'*cali_P + cali_P*param.model.A + cali_P*param.model.B*epsilon*K)
    eig(double(cali_P))

    
    
end

