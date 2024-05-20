function [param] = DesignProcedureSOF_FWMM(param, info, listOfUncertainties, alpha)
    % Generates a diagonal matrix for each vertex given the uncertainties.
    % Firstly every combination of uncertainties is generated, then
    % structured into a 3D array of diagonal matrices.
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

    A = param.model.A + (alpha)*eye(12);
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
    
    D_Bar = [-param.model.Dwz*param.model.DWref*param.model.Dww;
            zeros(4,4);];


    % Step 1-3
    q = cell(param.n);
    y = cell(param.n);
    %fill in p^-1
    for i = 1:param.n
        q{i} = sdpvar(3);
        y{i} = sdpvar(1,3,'full');
    end
    Q_hat = sdpvar(28,28);
    Y_hat = [sdpvar(4,28,'full')];
    %insert Q into Q_hat
    Q_hat(1:12,1:12) = blkdiag(q{:});
    %insert Y into y_hat
    Y_hat(1:4,1:12) = blkdiag(y{:});
    gamma = sdpvar(1);

    constraints = [param.model.A * Q_hat(1:12,1:12) + Q_hat(1:12,1:12) * param.model.A.' + param.model.B_Bar * Y_hat(1:4,1:12) + Y_hat(1:4,1:12).' * param.model.B_Bar.' <=0];
    for i=1:length(Uncertainties)
        constraints = [constraints, [Q_hat*A_bar.'+A_bar*Q_hat+Y_hat.'*B_bar(i).'+B_bar(i)*Y_hat, B_DS, Q_hat*C1_bar.'+Y_hat.'*C2_bar.';
                         B_DS.', -gamma*eye(4),D_Bar.';
                         C1_bar*Q_hat+C2_bar*Y_hat, D_Bar, -gamma*eye(8)]<=0];
    end
    constraints = [constraints,gamma >= 0, Q_hat >=0];
    
    options = sdpsettings('verbose',0,'solver','mosek');
    optimize(constraints,[gamma],options);
    gammaSF = value(gamma);
    y_out = value(Y_hat);
    q_out = value(Q_hat);
    
    % step 4
    Ksf = y_out(1:4,1:12)/q_out(1:12,1:12);
    
    %step 5
    yalmip('clear')
    p = cell(1,param.n);
    for i=1:param.n
        p{i} = sdpvar(3);
    end
    gamma = sdpvar(1);
    P = sdpvar(28);
    P(1:12,1:12) = blkdiag(p{:});
    sigma = sdpvar(1);
    constraints = [A.'*P(1:12,1:12) + P(1:12,1:12)*A - sigma*(param.model.CySOF.'*param.model.CySOF) <= 0];
    for i=1:length(Uncertainties)
        constraints = [constraints, [(A_bar + B_bar(i)*Ksf*[eye(12),zeros(12,16)]).'*P + P*(A_bar + B_bar(i)*Ksf*[eye(12),zeros(12,16)]), P * B_DS, (C1_bar+C2_bar*Ksf*[eye(12),zeros(12,16)]).';
                         B_DS.' * P, -gamma*eye(4),D_Bar.';
                         (C1_bar+C2_bar*Ksf*[eye(12),zeros(12,16)]), D_Bar, -gamma*eye(8)]<=0];
    end
    constraints = [constraints, P>=0, sigma >= 0];
    options = sdpsettings('verbose',0,'solver','mosek');
    optimize(constraints, [gamma], options);
    gammaSOF_P = value(gamma); 
    P_min = value(P);
    
    % Step 6
    yalmip('clear')
    k = cell(1,param.n);
    for j=1:param.n
        k{j} = sdpvar(1,2,'full');
    end
    K = blkdiag(k{:});
    gamma = sdpvar(1);
        
    kvec = [K(:,1);K(:,2);K(:,3);K(:,4);K(:,5);K(:,6);K(:,7);K(:,8)];
    constraints = [gamma, kvec.'; kvec, eye(size(kvec,1))]>=0;
    for i=1:length(Uncertainties)
        constraints = [constraints, [(A + B(i)*K*param.model.CySOF).'*P_min(1:12,1:12) + P_min(1:12,1:12)*(A + B(i)*K*param.model.CySOF)]<=0];
    end
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints, gamma, options);
    
    Ksof = value(K);
    Ks = zeros(1,3,param.n);
    for j=1:param.n
        Ks(1,2:3,j) = Ksof(j,2*(j-1)+1:2*j);
    end
    %saves a 3d matrix version of Ks this is used in the simulation of the nonlinear system.
    param.ctrl.Ks = Ks;

    % Evaluation of the results from the design procedure
    fprintf("Evaluation of design procedure for three-step SOF \n")
    fprintf("Gamma for SF is found to be: %f \n", value(gammaSF))
    fprintf("Gamma for P SOF is found to be: %f \n", value(gammaSOF_P))
    if(eig(param.model.A+param.model.B_Bar*Ksof*param.model.CySOF) > 0)
        fprintf("One or more eigenvalues in the decoupled system is positive!!! \n")
    else
        fprintf("Eigenvalues of the decoupled systems are all negative \n")
    end
    eig(param.model.A+param.model.B_Bar*Ksof*param.model.CySOF).'
    
    if(eig(param.model.A+param.model.B*Ksof*param.model.CySOF) > 0)
        fprintf("One or more eigenvalues in the coupled system is positive!!! \n")
    else
        fprintf("Eigenvalues of the coupled systems are all negative \n")
    end
    eig(param.model.A+param.model.B*Ksof*param.model.CySOF).'
    %saves the block version of K this is used further in some simulations
    param.ctrl.K = Ksof*param.model.CySOF;
    %saves a 3d matrix version of Ks this is used in the simulation of the
    %nonlinear system.
    param.ctrl.Ks = Ks;
end

