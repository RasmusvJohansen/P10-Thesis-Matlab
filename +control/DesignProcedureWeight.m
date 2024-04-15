function [param] = DesignProcedureWeight(param,info,isGammaOne)
% Synthesises a controller based on design procedure 2 with weight in the report, and
% returns it as an 3D-array with each 2D controller ocyping a slot in the
% 3rd dimension. isGammaOne is used to determine if gamma=1 or if it is
% an optimisation variable

    % Initializes yalmip and optimisation variables.
    yalmip('clear')
    q = cell(1,param.n + 1);
    y = cell(1,param.n);
    for i=1:param.n
        q{i} = sdpvar(3);
        y{i} = sdpvar(1,3,'full');
    end
    q{end} = sdpvar(4);
    Q_small = blkdiag(q{1:end-1});
    Q = blkdiag(Q_small,q{end});
    Y_small = blkdiag(y{:});
    Y = blkdiag(Y_small, sdpvar(4));

    % Is gamma = 1 or an optimisation variable
    if isGammaOne
        gamma = 1;
    else
        gamma = sdpvar(1);
    end
    
    % % Structure the matrices
    % A_split = [param.model.A, zeros(12,4); zeros(4,12), param.model.Awd];
    % B_split = [param.model.B, zeros(12,4); param.model.Bwd, zeros(4,4)];
    % % K_split = [K, zeros(4); zeros(4), eye(4)]
    % B_mw = [param.model.B ; zeros(4)];
    % C_split = [param.model.Dwd, param.model.Cwd];
    % D_mw = zeros(4);


    % Structure the matrices
    A_split = [param.model.A, zeros(12,4); zeros(4,12), param.model.Awd];
    B_split = [param.model.B, param.model.B * param.model.Cwd; zeros(4,4), zeros(4,4)];
    % K_split = [K, zeros(4); zeros(4), eye(4)]
    B_mw = [param.model.B * param.model.Dwd; param.model.Bwd];
    C_split = [eye(4), zeros(4)];
    D_mw = zeros(4);

    % Setup constraint
    constraints = param.model.A * Q_small + Q_small * param.model.A.' + param.model.B_Bar * Y_small + Y_small.' * param.model.B_Bar.' <=0;
    constraints = [constraints, [A_split * Q + Q * A_split.' + B_split * Y + Y.' * B_split.', B_mw, Y.' * C_split.';
                    B_mw.', -gamma * eye(4), D_mw;
                    C_split * Y, D_mw.', -gamma * eye(4)] <= 0];
    constraints = [constraints,  Q >= 0, gamma >= 0];
    options = sdpsettings('verbose',0,'solver','mosek');
    if info
        sol = optimize(constraints, gamma, options)
        check(constraints)
    else
        sol = optimize(constraints, gamma, options);
    end
    
    % Recovers the controller from the variables
    K = value(Y_small)/value(Q_small);
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
    %saves the block version of K this is used further in some simulations
    param.ctrl.K = K;
    %saves a 3d matrix version of Ks this is used in the simulation of the
    %nonlinear system.
    param.ctrl.Ks = Ks;
end