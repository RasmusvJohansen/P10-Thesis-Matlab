function [Gamma] = CalculatePerformanceLPV(param,listOfUncertainties)
    m = length(listOfUncertainties);
    Uncertainties = zeros(m,m,m^2);
    if(m == 1)% check if the inputs share uncertainties or if they have individual uncertainties
        comb = combinations(cell2mat(listOfUncertainties(1)));
    else
        comb = combinations(cell2mat(listOfUncertainties(1)),cell2mat(listOfUncertainties(2)),cell2mat(listOfUncertainties(3)),cell2mat(listOfUncertainties(4)));
    end
    combDouble = table2array(comb);
    for i=1:length(combDouble)
        Uncertainties(:,:,i) = diag(combDouble(i,:));
    end

    B = @(i) param.model.B*(eye(4) + Uncertainties(i) * param.model.Dwd);
    A_CL = @(k) [param.model.A+B(k)* param.ctrl.K, param.model.Bw*param.model.Cww, zeros(12,4);
                zeros(4,12), param.model.Aww, zeros(4,4);
                param.model.Bwz*param.model.Cz, -param.model.Bwz*param.model.Cww, param.model.Awz];
    
    B_CL = [param.model.Bw*param.model.Dww; param.model.Bww; -param.model.Bwz*param.model.Dww];
    C_CL = [param.model.Dwz*param.model.Cz, -param.model.Dwz*param.model.Cww, param.model.Cwz];
    D_CL = [-param.model.Dwz*param.model.Dww];
    
    yalmip('clear')
    p = cell(1,param.n);
    for i=1:param.n
        p{i} = sdpvar(3);
        ksof{i} = sdpvar(1,2,'full');
    end
    P = sdpvar(20);
    P(1:12,1:12) = blkdiag(p{1:end});
    gamma = sdpvar(1);
    constraints = [];
    for i=1:m
        constraints = [constraints, [A_CL(i).'*P + P*A_CL(i), P*B_CL, C_CL.';
                                     B_CL.'*P, -gamma*eye(4), D_CL.';
                                     C_CL,  D_CL, -gamma*eye(4) ]<=0];
    end
    constraints = [constraints, P>=0, gamma >=0];
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints, gamma, options);
    Gamma = value(gamma)

    % Calculate Hinfnorm from the TF
    A_CL_Nominal = [param.model.A+param.model.B* param.ctrl.K, param.model.Bw*param.model.Cww, zeros(12,4);
                zeros(4,12), param.model.Aww, zeros(4,4);
                param.model.Bwz*param.model.Cz, -param.model.Bwz*param.model.Cww, param.model.Awz];
    s=tf('s');
    N = C_CL/(s*eye(20)-A_CL_Nominal)*B_CL+D_CL;
    ActualGamma = hinfnorm(N)
end