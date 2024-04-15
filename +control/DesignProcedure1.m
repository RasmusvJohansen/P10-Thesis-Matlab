function [param] = DesignProcedure1(param,info)
% Synthesises a controller based on design procedure 1 in the report, and
% returns it as an 3D-array with each 2D controller ocyping a slot in the
% 3rd dimension.

    % Initializes yalmip and optimisation variables.
    yalmip('clear')
    q = cell(1,param.n);
    y = cell(1,param.n);
    r = cell(1,param.n);
    for i=1:param.n
        q{i} = sdpvar(3);
        y{i} = sdpvar(1,3,'full');
        r{i} = sdpvar(3); 
    end
    Q = blkdiag(q{:});
    Y = blkdiag(y{:});
    R = blkdiag(r{:});
    gamma = sdpvar(1);

    % Setup and solve for stability in the decoupled system
    constraints = [param.model.A*Q+Q*param.model.A.'+param.model.B_Bar *Y + Y.'* param.model.B_Bar.' <=0];
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
    
    % Clear constraints and set up for robust stability of the coupled
    % system
    clear constraints
    constraints = [R*(param.model.A + param.model.B * K * eye(12)) + (param.model.A + param.model.B * K * eye(12)).' * R, R * param.model.B, eye(12)*K.';
                    param.model.B.' * R, -gamma * eye(4), zeros(4,4);
                    K*eye(12), zeros(4,4), -gamma * eye(4)] <= 0;
    constraints = [constraints,  R >= 0, gamma >= 0];
    options = sdpsettings('verbose',0,'solver','mosek');
    if info
        sol = optimize(constraints, gamma, options)
        check(constraints)
    else
        sol = optimize(constraints, gamma, options);
    end

    % Evaluation of the results from the design procedure
    fprintf("Evaluation of design procedure 1 \n")
    if(eig(value(Q))<=0)
        fprintf("Q is not positive definite!!!\n")
    else
        fprintf("Q is positive definite\n")
    end
    eig(value(Q)).'

    if(eig(value(R))<=0)
        fprintf("R is not positive definite!!!\n")
    else
        fprintf("R is positive definite\n")
    end
    eig(value(R)).'
    
    fprintf("Gamma is found to be:\n")
    value(gamma)
    
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