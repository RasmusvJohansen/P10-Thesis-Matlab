function [param] = DesignProcedureSOF_ILMI_Hinf(param, info, listOfUncertainties)
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
    Q = eye(20);
    tol = 1e-2;
    tol_alpha = 0.005;
    n = param.n;
    yalmip('clear')
    gamma = sdpvar(1);
    A_1 = [param.model.A, param.model.Bw*param.model.Cww, zeros(12,4);
            zeros(4,12), param.model.Aww, zeros(4,4);
            param.model.Bwz*param.model.Cz, -param.model.Bwz*param.model.Cww, param.model.Awz];
    B_1 = @(k) [B(k); zeros(4,4); zeros(4,4)];
    C_1 = [param.model.CySOF, zeros(8,4), zeros(8,4)];
    B_CL = [param.model.Bw*param.model.Dww; param.model.Bww; -param.model.Bwz*param.model.Dww];
    C_CL = [param.model.Dwz*param.model.Cz, -param.model.Dwz*param.model.Cww, param.model.Cwz];
    D_CL = [-param.model.Dwz*param.model.Dww];
    A_bar = [A_1, B_CL, zeros(20,4);
            zeros(4,20), -gamma/2*eye(4), zeros(4,4);
            C_CL, D_CL, -gamma/2*eye(4)];
    B_bar = @(k) [B_1(k); zeros(4,4); zeros(4,4)];
    C_bar = [C_1, zeros(8,4), zeros(8,4)];


    % Step 1
    % [P_bar, ~, P] = constructVariables(param);
    p = cell(1,param.n);
    for i=1:param.n
        p{i} = sdpvar(3);
    end
    P = sdpvar(20);
    P(1:12,1:12) = blkdiag(p{1:end})
    % A_bar
    % P_bar
    % Q
    % B_bar(1)
    constraints = [];
    for k=1:length(Uncertainties)
        constraints = [constraints, [A_1.' * P + P * A_1 + Q, P * B_1(k); B_1(k).' * P, eye(n)] >= 0];
    end
    % 
    % constraints = [A.' * P + P * A + Q,P * B;
    %                 B.' * P, eye(n)] >= 0;
    constraints = [constraints, P >= 0];
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints, -trace(P), options)
    X(:,:,1) = blkdiag(value(P),eye(8))
    i = 1;
    
    while 1
        % Step 2
        yalmip('clear')
        [P_bar, Ksof] = constructVariables(param);
        gamma = sdpvar(1);
         A_bar = [A_1, B_CL, zeros(20,4);zeros(4,20), -gamma/2*eye(4), zeros(4,4);C_CL, D_CL, -gamma/2*eye(4)];
        
        alpha = sdpvar(1);
        constraints = [];
        for k=1:length(Uncertainties)
            % constraints = [constraints, [A.'*P + P*A - X(:,:,i)*(B(k)*B(k).')*P - P*(B(k)*B(k).')*X(:,:,i) + X(:,:,i)*(B(k)*B(k).')*X(:,:,i) - alpha*P , (B(k).'*P + Ksof*C).';
                            % (B(k).'*P + Ksof*C), -eye(n)] <= 0];
            % constraints = [constraints, [P_bar*B_bar(k)*Ksof*C_bar + (P_bar*B_bar(k)*Ksof*C_bar).' + A_bar.'*P_bar + P_bar*A_bar]<=0];
            constraints = [constraints, [A_bar.'*P_bar + P_bar*A_bar - X(:,:,i)*(B_bar(k)*B_bar(k).')*P_bar - P_bar*(B_bar(k)*B_bar(k).')*X(:,:,i) + X(:,:,i)*(B_bar(k)*B_bar(k).')*X(:,:,i) - alpha*P_bar , (B_bar(k).'*P_bar + Ksof*C_bar).';
                            (B_bar(k).'*P_bar + Ksof*C_bar), -eye(n)] <= 0];
        end
        constraints = [constraints, P_bar >= 0];
        options = sdpsettings('verbose',0,'solver','mosek');
        sol = bisection(constraints, alpha, options);
        alpha_min = value(alpha)
        % Step 3
        if(alpha_min <= tol_alpha)
            break
        end
        % Step 4
        yalmip('clear')
        [P_bar, Ksof] = constructVariables(param);
        gamma = sdpvar(1);
        A_bar = [A_1, B_CL, zeros(20,4);zeros(4,20), -gamma/2*eye(4), zeros(4,4);C_CL, D_CL, -gamma/2*eye(4)];
        constraints = [];
        for k=1:length(Uncertainties)
            % constraints = [constraints, [A.'*P + P*A - X(:,:,i)*(B(k)*B(k).')*P - P*(B(k)*B(k).')*X(:,:,i) + X(:,:,i)*(B(k)*B(k).')*X(:,:,i) - alpha_min*P , (B(k).'*P + Ksof*C).';
                            % (B(k).'*P + Ksof*C), -eye(n)] <= 0];
        constraints = [constraints, [A_bar.'*P_bar + P_bar*A_bar - X(:,:,i)*(B_bar(k)*B_bar(k).')*P_bar - P_bar*(B_bar(k)*B_bar(k).')*X(:,:,i) + X(:,:,i)*(B_bar(k)*B_bar(k).')*X(:,:,i) - alpha_min*P_bar , (B_bar(k).'*P_bar + Ksof*C_bar).';
                            (B_bar(k).'*P_bar + Ksof*C_bar), -eye(n)] <= 0];
        end
        constraints = [constraints, P_bar >= 0];
        options = sdpsettings('verbose',0,'solver','mosek');
        sol = optimize(constraints, trace(P_bar), options);
        P_min(:,:,i) = value(P_bar);
        % Step 5
        if (norm(X(:,:,i) - P_min(:,:,i)) < tol)
            % Step 6
            break
        else
            i = i + 1;
            X(:,:,i) = P_min(:,:,i-1);
        end
    end
    
    Ksof = value(Ksof);
    Ks = zeros(1,3,param.n);
    for j=1:param.n
        Ks(1,2:3,j) = Ksof(j,2*(j-1)+1:2*j);
    end
    
    % Evaluation of the results from the design procedure
    fprintf("Evaluation of design procedure for ILMI SOF \n")
    fprintf("Gamma is found to be: %f \n", value(gamma))
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

function [P_bar, Ksof, P] = constructVariables(param)
    ksof = cell(1,param.n);
    p = cell(1,param.n);
    for i=1:param.n
        p{i} = sdpvar(3);
        ksof{i} = sdpvar(1,2,'full');
    end
    P = sdpvar(20);
    P(1:12,1:12) = blkdiag(p{1:end});
    P_bar = blkdiag(P,eye(4),eye(4));
    Ksof = blkdiag(ksof{:});
end