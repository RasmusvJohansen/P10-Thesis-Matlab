function param = DesignProcedureFWMM(param,info)
    % Design procedure for frequency weighted model matching
    yalmip('clear')
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
    
    %constraints
    constriants = [];
    
    %Define the lyapunaov constraint for the decoupled system 
    Lyap = param.model.A * Q_hat(1:param.n*3,1:param.n*3) + Q_hat(1:param.n*3,1:param.n*3) * param.model.A.' + param.model.B_Bar * Y_hat(1:param.n,1:param.n*3) + Y_hat(1:param.n,1:param.n*3).' * param.model.B_Bar.';
    constriants = [constriants, Lyap <= 0, Q_hat(1:param.n*3,1:param.n*3) >= 0];
    
    A_bar = [param.model.A, param.model.Bw*param.model.Cww, zeros(12,4),zeros(12,4),zeros(12,4);
             zeros(4,12),param.model.Aww, zeros(4),zeros(4),zeros(4);
             param.model.Bwz*param.model.Cz, -param.model.Bwz*param.model.DWref*param.model.Cww, param.model.Awz, -param.model.Bwz*param.model.CWref,zeros(4);
             zeros(4,12),param.model.BWref*param.model.Cww, zeros(4),param.model.AWref,zeros(4);
             zeros(4,12), zeros(4),zeros(4),zeros(4),param.model.AWu]; 
    
    B_bar = [param.model.B;
             zeros(4,4);
             zeros(4,4);
             zeros(4,4);
             param.model.BWu];
    
    B_DS = [param.model.Bw*param.model.Dww, param.model.B;
             param.model.Bww, zeros(4);
             -param.model.Bwz*param.model.DWref*param.model.Dww, zeros(4);
             param.model.BWref*param.model.Dww,zeros(4);
             zeros(4,8)];
    
    C1_bar = [param.model.Dwz*param.model.Cz, -param.model.Dwz*param.model.DWref*param.model.Cww, param.model.Cwz,-param.model.Dwz*param.model.CWref,zeros(4);
          zeros(4,12),zeros(4),zeros(4),zeros(4),param.model.CWu;
          zeros(4,28)];
    
    C2_bar = [zeros(4,4);
          param.model.DWu;
          param.model.Dwd];
    
    D_DS = [-param.model.Dwz*param.model.DWref*param.model.Dww,zeros(4);
            zeros(4,8);
            zeros(4,8)];
    
    
    DesiredPerformance = [Q_hat*A_bar.'+A_bar*Q_hat+Y_hat.'*B_bar.'+B_bar*Y_hat, B_DS, Q_hat*C1_bar.'+Y_hat.'*C2_bar.';
                         B_DS.', -gamma*eye(8),D_DS.';
                         C1_bar*Q_hat+C2_bar*Y_hat, D_DS, -gamma*eye(12)];
    constriants = [constriants, DesiredPerformance <= 0,gamma >= 0, Q_hat >=0];
    
    options = sdpsettings('verbose',0,'solver','mosek');    
    if(info)
        sol = optimize(constriants,[gamma],options)
        check(constriants)
    else
        sol = optimize(constriants,[gamma],options);
    end
    
    y_out = value(Y_hat);
    q_out = value(Q_hat);
    k = y_out(1:4,1:12)/q_out(1:12,1:12);
    
    Ks = zeros(1,3,param.n);
      for i=1:param.n
            Ks(:,:,i) = k(i,3*(i-1)+1:3*i);
      end
    param.ctrl.Ks = Ks;
    param.ctrl.K = k;

    % Evaluation of the results from the design procedure
    fprintf("Evaluation of design procedure for Frequency weighted model matching \n")
    fprintf("Gamma is found to be: %f \n", value(gamma))
    if(info)
        disp("Eig of Q_hat")
        eig(value(Q_hat))
    end
end

