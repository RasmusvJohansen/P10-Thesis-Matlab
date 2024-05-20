function [param] = DesingProcedureFWMM_Polytopic(param,info,listOfUncertainties)
    % Design procedure for frequency weighted model matching with polytopic
    % uncertainties
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
    A_1 = [param.model.A, param.model.Bw*param.model.Cww, zeros(12,4),zeros(12,4),zeros(12,4);
             zeros(4,12),param.model.Aww, zeros(4),zeros(4),zeros(4);
             param.model.Bwz*param.model.Cz, -param.model.Bwz*param.model.DWref*param.model.Cww, param.model.Awz, -param.model.Bwz*param.model.CWref,zeros(4);
             zeros(4,12),param.model.BWref*param.model.Cww, zeros(4),param.model.AWref,zeros(4);
             zeros(4,12), zeros(4),zeros(4),zeros(4),param.model.AWu]; 
    B_2 = @(k) [B(k);
                 zeros(4,4);
                 zeros(4,4);
                 zeros(4,4);
                 param.model.BWu];
    B_1 = [param.model.Bw*param.model.Dww;
             param.model.Bww;
             -param.model.Bwz*param.model.DWref*param.model.Dww;
             param.model.BWref*param.model.Dww;
             zeros(4,4)];
    C_1 = [param.model.Dwz*param.model.Cz, -param.model.Dwz*param.model.DWref*param.model.Cww, param.model.Cwz,-param.model.Dwz*param.model.CWref,zeros(4);
          zeros(4,12),zeros(4),zeros(4),zeros(4),param.model.CWu];
    D_12 = [zeros(4,4);
          param.model.DWu];
    D_11 = [-param.model.Dwz*param.model.DWref*param.model.Dww;
            zeros(4,4);];
    
    n = param.n;
    C_2 = [param.model.CySOF, zeros(8,4), zeros(8,4), zeros(8,4), zeros(8,4)];
    
    B_bar = @(k) [B_2(k); zeros(4,4); D_12];
    C_bar = [C_2, zeros(8,8), zeros(8,4)];

    q = cell(param.n);
    y = cell(param.n);
    for i = 1:param.n
        q{i} = sdpvar(3);
        y{i} = sdpvar(1,3,'full');
    end
    Q_hat = sdpvar(28,28);
    Y_hat = [sdpvar(4,28,'full')];
    Q_hat(1:12,1:12) = blkdiag(q{:});
    Y_hat(1:4,1:12) = blkdiag(y{:});
    gammasf = sdpvar(1);
    constraints = [];
    for i=1:length(Uncertainties)
        constraints = [constraints, [Q_hat*A_1.'+A_1*Q_hat+Y_hat.'*B_2(i).'+B_2(i)*Y_hat, B_1, Q_hat*C_1.'+Y_hat.'*D_12.';
                         B_1.', -gammasf*eye(4),D_11.';
                         C_1*Q_hat+D_12*Y_hat, D_11, -gammasf*eye(8)]<=0];
    end
    constraints = [constraints,gammasf >= 0, Q_hat >=0];
    
    options = sdpsettings('verbose',0,'solver','mosek');
    optimize(constraints,[gammasf],options);
    k = value(Y_hat(1:4,1:12))/(value(Q_hat(1:12,1:12)));
       Ks = zeros(1,3,param.n);
      for i=1:param.n
            Ks(:,:,i) = k(i,3*(i-1)+1:3*i);
      end
    param.ctrl.Ks = Ks;
    param.ctrl.K=k;

    % Evaluation of the results from the design procedure
    fprintf("Evaluation of design procedure for FWMM uncertainties\n")
    fprintf("Gamma is found to be: %f \n", value(gammasf))
    if(info)
        disp("Eig of Q_hat")
        eig(value(Q_hat))
    end
end

