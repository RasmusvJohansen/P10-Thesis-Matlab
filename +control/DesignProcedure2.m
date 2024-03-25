function [Ks] = DesignProcedure2(param,info,isGammaOne)
% Synthesises a controller based on design procedure 1 in the report, and
% returns it as an 3D-array with each 2D controller ocyping a slot in the
% 3rd dimension. isGammaOne is used to determine if gamma=1 or if it is
% an optimisation variable

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

    % Is gamma = 1 or an optimisation variable
    if isGammaOne
        gamma = 1;
    else
        gamma = sdpvar(1);
    end

    % Setup constraint
    constraints = [param.model.A * Q + Q * param.model.A.' + param.model.B_Bar * Y + Y.' * param.model.B_Bar.' <=0];
    constraints = [param.model.A * Q + Q * param.model.A.' + param.model.B * Y + Y.' * param.model.B.', param.model.B, Y.';
                    param.model.B.', -gamma * eye(4), zeros(4,4);
                    Y, zeros(4,4), -gamma * eye(4)] <= 0;
    constraints = [constraints,  Q >= 0, gamma >= 0];
    options = sdpsettings('verbose',0,'solver','mosek');
    sol = optimize(constraints, gamma, options);
    if info
        check(constraints);
    end

    % Recovers the controller from the variables
    K = value(Y)/value(Q);
    Ks = zeros(1,3,param.n);
    for i=1:param.n
        Ks(:,:,i) = K(i,3*(i-1)+1:3*i);
    end
    
    % Evaluation of the results from the design procedure
    fprintf("Evaluation of design procedure 2 \n")
    if(eig(value(Q))<=0)
        fprintf("Q is not positive definite!!!\n")
    else
        fprintf("Q is positive definite\n")
    end
    eig(value(Q)).'

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
    
end