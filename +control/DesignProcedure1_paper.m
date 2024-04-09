function param = DesignProcedure1_paper(param, info)
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
    m = 1;
    res = zeros(1,m)
    eps = linspace(0,1,m);
    for i = 1:m
        clear constraints
        
        eps1 = eps(i); 
        eps2 = eps(i); 
        eps3 = eps(i);
        eps4 = eps(i);
        epsilon =[eps1 0 0 0; 0 eps2 0 0; 0 0 eps3 0; 0 0 0 eps4];
    
        constraints = [param.model.A.'*cali_P + K.'*epsilon.'*param.model.B_Bar.'*cali_P + cali_P*param.model.A + cali_P*param.model.B_Bar*epsilon*K <=0];
        constraints = [constraints, cali_P >= 0];
    
        %constraints = []
        options = sdpsettings('verbose',0,'solver','mosek');
        sol = optimize(constraints,[],options);
        if info
            check(constraints);
        end
        if(min(eig(value(cali_P))) < 0)
            eig(value(cali_P));
            res(i) = eps(i);
        end
            %eig(param.model.A.'*value(cali_P) + K.'*epsilon.'*param.model.B_Bar.'*value(cali_P) + value(cali_P)*param.model.A + value(cali_P)*param.model.B_Bar*epsilon*K)
        
       

    end
    %saves the block version of K this is used further in some simulations
    param.ctrl.K = K;
    %saves a 3d matrix version of Ks this is used in the simulation of the
    %nonlinear system.
    param.ctrl.Ks = Ks;

end

