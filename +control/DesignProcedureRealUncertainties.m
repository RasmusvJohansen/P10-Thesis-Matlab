function [param] = DesignProcedureRealUncertainties(param, info, listOfUncertainties)
% Synthesises a controller based on design procedure 1 in the report, and
% returns it as an 3D-array with each 2D controller ocyping a slot in the
% 3rd dimension. listOfUncertainties is an array contining the
% uncertainties and their range

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

    % Initializes yalmip and optimisation variables.
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
    constraints = [param.model.A * Q + Q * param.model.A.' + param.model.B_Bar * Y + Y.' * param.model.B_Bar.' <=0];
    for i=1:length(Uncertainties)
        constraints = [constraints, param.model.A * Q + Q * param.model.A.' + param.model.B * (eye(m) + Uncertainties(i)) * Y + Y.' * (eye(m) + Uncertainties(i))' * param.model.B.' <= 0];
    end
    constraints = [constraints,  Q >= 0];
    options = sdpsettings('verbose',0,'solver','mosek');
    
    if info
        sol = optimize(constraints, [], options)
        check(constraints)
    else
        sol = optimize(constraints, [], options);
    end

    % Recovers the controller from the variables
    K = value(Y)/value(Q);
    Ks = zeros(1,3,param.n);
    for i=1:param.n
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end
    
    % Evaluation of the results from the design procedure
    fprintf("Evaluation of design procedure for real uncertainties \n")
    if(eig(value(Q))<=0)
        fprintf("Q is not positive definite!!!\n")
    else
        fprintf("Q is positive definite\n")
    end
    eig(value(Q)).'

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
    %saves the block version of K this is used further in some simulations
    param.ctrl.K = K;
    %saves a 3d matrix version of Ks this is used in the simulation of the
    %nonlinear system.
    param.ctrl.Ks = Ks;
end