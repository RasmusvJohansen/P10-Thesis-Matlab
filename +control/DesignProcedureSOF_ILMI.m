function [param] = DesignProcedureSOF_ILMI(param, info, listOfUncertainties)
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
    Q = eye(12);
    tol = 1e-2;
    tol_alpha = 0.001;
    n = param.n;
    % Step 1
    p = cell(1,param.n);
    for i=1:param.n
        p{i} = sdpvar(3);
    end
    P = blkdiag(p{1:end});

    constraints = [];
    for k=1:length(Uncertainties)
        constraints = [constraints, [param.model.A.' * P + P * param.model.A + Q, P * B(k); B(k).' * P, eye(n)] >= 0];
    end
    % 
    % constraints = [A.' * P + P * A + Q,P * B;
    %                 B.' * P, eye(n)] >= 0;
    constraints = [constraints, P >= 0];
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints, -trace(P), options);
    X(:,:,1) = value(P);
    i = 1;
    
    while 1
        % Step 2
        yalmip('clear')
        [P, Ksof] = constructVariables(param);
        alpha = sdpvar(1);
        constraints = [];
        for k=1:length(Uncertainties)
            constraints = [constraints, [param.model.A.'*P + P*param.model.A - X(:,:,i)*(B(k)*B(k).')*P - P*(B(k)*B(k).')*X(:,:,i) + X(:,:,i)*(B(k)*B(k).')*X(:,:,i) - alpha*P , (B(k).'*P + Ksof*param.model.CySOF).';
                            (B(k).'*P + Ksof*param.model.CySOF), -eye(n)] <= 0];
        end
        constraints = [constraints, P >= 0];
        options = sdpsettings('verbose',0,'solver','mosek');
        sol = bisection(constraints, alpha, options);
        alpha_min(i) = value(alpha);
        alpha_min(i)
        % Step 3
        if(alpha_min(i) <= tol_alpha)
            break
        end
        if (isnan(alpha_min(i)))
            fprintf("Alpha became NaN after %i itr. Last value was %f \n", i, alpha_min(i-1));
            i = i-1;
            alpha_min(i) = alpha_min(i) + 0.001; % 0.01
        end
        % Step 4
        yalmip('clear')
        [P, Ksof] = constructVariables(param);
        constraints = [];
        for k=1:length(Uncertainties)
            constraints = [constraints, [param.model.A.'*P + P*param.model.A - X(:,:,i)*(B(k)*B(k).')*P - P*(B(k)*B(k).')*X(:,:,i) + X(:,:,i)*(B(k)*B(k).')*X(:,:,i) - alpha_min*P , (B(k).'*P + Ksof*param.model.CySOF).';
                            (B(k).'*P + Ksof*param.model.CySOF), -eye(n)] <= 0];
        end
        constraints = [constraints, P >= 0];
        options = sdpsettings('verbose',0,'solver','mosek');
        sol = optimize(constraints, trace(P), options);
        P_min(:,:,i) = value(P);
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
    fprintf("Solved the problem in %i iterations\n", i)
    fprintf("The solution had alpha: %f \n", alpha_min)
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

function [P, Ksof] = constructVariables(param)
    ksof = cell(1,param.n);
    p = cell(1,param.n);
    for i=1:param.n
        p{i} = sdpvar(3);
        ksof{i} = sdpvar(1,2,'full');
    end
    P = blkdiag(p{1:end});
    Ksof = blkdiag(ksof{:});
end