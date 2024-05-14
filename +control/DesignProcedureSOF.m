function [param] = DesignProcedureSOF(param, info, listOfUncertainties, alpha)
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

    % Step 1-3
    A = param.model.A + (alpha)*eye(12);

    yalmip('clear')
    q = cell(1,param.n);
    y = cell(1,param.n);
    for i=1:param.n
        q{i} = sdpvar(3);
        y{i} = sdpvar(1,3,'full');
    end
    Q = blkdiag(q{:});
    Y = blkdiag(y{:});


    % Setup constraint
    constraints = [];
    for i=1:length(Uncertainties)
        constraints = [constraints, A * Q + Q * A.' + B(i) * Y + Y.' * B(i).' <= 0];
    end
    constraints = [constraints,  Q - eye(12)>= 0];
    % constraints = [constraints,  Q >= 0];
    
    options = sdpsettings('verbose',0,'solver','mosek');
    if info
        optimize(constraints, [], options)
        check(constraints)
    else
        optimize(constraints, [], options);
    end
    
    % step 4
    Ksf = value(Y)/value(Q);
    
    %step 5
    yalmip('clear')
    p = cell(1,param.n);
    for i=1:param.n
        p{i} = sdpvar(3);
    end
    P = blkdiag(p{:});
    sigma = sdpvar(1);
    constraints = [A.'*P + P*A - sigma*(param.model.CySOF.'*param.model.CySOF) <= 0];
    for i=1:length(Uncertainties)
        constraints = [constraints, (A + B(i)*Ksf).'*P + P*(A + B(i)*Ksf) <= 0];
    end
    constraints = [constraints, P - eye(12) >= 0, sigma >= 0];
    options = sdpsettings('verbose',0,'solver','mosek');
    if info
        sol = optimize(constraints, [], options)
        check(constraints)
    else
        sol = optimize(constraints, [], options);
    end
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
    % constraints = [(A + B(2)*K*param.model.CySOF).'*P_min + P_min*(A + B(2)*K*param.model.CySOF)]<=0;
    % constraints = [constraints, [(A + param.model.B*K*param.model.CySOF).'*P_min + P_min*(A + param.model.B*K*param.model.CySOF)]<=0];
    for i=1:length(Uncertainties)
        constraints = [constraints, [(A + B(i)*K*param.model.CySOF).'*P_min + P_min*(A + B(i)*K*param.model.CySOF)]<=0];
    end
    
    options = sdpsettings('verbose',0,'solver','mosek');
    if info
        sol = optimize(constraints, gamma, options)
        check(constraints)
    else
        sol = optimize(constraints, gamma, options);
    end
    Ksof = value(K);
    
    Ks = zeros(1,3,param.n);
    for j=1:param.n
        Ks(1,2:3,j) = Ksof(j,2*(j-1)+1:2*j);
    end

    % Evaluation of the results from the design procedure
    fprintf("Evaluation of design procedure for three-step SOF \n")
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

