function [Gamma] = CalculateFWMM(param,listOfUncertainties)
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


    A = param.model.A + (0)*eye(12);
    A_bar = [param.model.A, param.model.Bw*param.model.Cww, zeros(12,4),zeros(12,4),zeros(12,4);
             zeros(4,12),param.model.Aww, zeros(4),zeros(4),zeros(4);
             param.model.Bwz*param.model.Cz, -param.model.Bwz*param.model.DWref*param.model.Cww, param.model.Awz, -param.model.Bwz*param.model.CWref,zeros(4);
             zeros(4,12),param.model.BWref*param.model.Cww, zeros(4),param.model.AWref,zeros(4);
             zeros(4,12), zeros(4),zeros(4),zeros(4),param.model.AWu]; 
    
    B_bar = @(k) [B(k);
                 zeros(4,4);
                 zeros(4,4);
                 zeros(4,4);
                 param.model.BWu];
    
    B_DS = [param.model.Bw*param.model.Dww;
             param.model.Bww;
             -param.model.Bwz*param.model.DWref*param.model.Dww;
             param.model.BWref*param.model.Dww;
             zeros(4,4)];
    
    C1_bar = [param.model.Dwz*param.model.Cz, -param.model.Dwz*param.model.DWref*param.model.Cww, param.model.Cwz,-param.model.Dwz*param.model.CWref,zeros(4);
          zeros(4,12),zeros(4),zeros(4),zeros(4),param.model.CWu];
    
    C2_bar = [zeros(4,4);
          param.model.DWu];
    
    D_DS = [-param.model.Dwz*param.model.DWref*param.model.Dww;
            zeros(4,4);];


    A_CL = @(k) A_bar + B_bar(k)*[param.ctrl.K,zeros(4,16)];
    
    B_CL = B_DS;
    C_CL = C1_bar+C2_bar*[param.ctrl.K,zeros(4,16)];
    D_CL = D_DS;
    
    yalmip('clear')
    p = cell(1,param.n);
    for i=1:param.n
        p{i} = sdpvar(3);
    end
    P = sdpvar(28);
    P(1:12,1:12) = blkdiag(p{1:end});
    gamma = sdpvar(1);
    constraints = [];
    for i=1:m
        constraints = [constraints, [A_CL(i).'*P + P*A_CL(i), P*B_CL, C_CL.';
                                     B_CL.'*P, -gamma*eye(4), D_CL.';
                                     C_CL,  D_CL, -gamma*eye(8) ]<=0];
    end
    constraints = [constraints, P>=0, gamma >=0];
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints, gamma, options);
    Gamma = value(gamma);
end