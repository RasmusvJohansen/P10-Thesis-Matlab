function [dx,q] = coupledDynamics(t,x,param,step_time)    
    % States
    T_w = x(1:4);
    T_a = x(5:8);
    
    % Inputs
    T_c = param.thermo.T_c;
    T_amb = param.thermo.T_A;
    Q = param.ctrl.Q;

    % Initialize variables used in loops
    dT_w = zeros(4,1);
    dT_a = zeros(4,1);
    dT = zeros(4,1);
    omega = zeros(4,1);

    % Controller
    for i=1:4
        if nargin > 3
            if t>step_time % Perform step if time
                T_ref = param.ctrl.T_ref + 1;
            else
                T_ref = param.ctrl.T_ref;
            end
        else
            T_ref = param.ctrl.T_ref;
        end
        dT(i) = T_a(i) - T_ref;
    
        % Apply linearisation and repack x to only contain the current
        % loops variables
        xlin = [T_w(i) - param.ctrl.T_wOP(i); T_a(i) - param.ctrl.T_aOP(i); x(8+i)];
        % Apply the control law for the coupleSystem (25)+(36,5)
        omega(i) = param.ctrl.Ks(:,:,i) * xlin + param.ctrl.omega_OP(i);
    end
    % Hydraulic network 
    % The flow could be calculated from the inverse funktion g(omega), 
    % however this becomes a troublesome function. Instead an the flow is 
    % found through simulation of the hydraulics at the given omega)
        % The simulation should be long enought to reach steadystate,
        % the desired flow can then be taken from the last element
    [~, q_sim] = ode15s(@(tt,xx)model.calculateFlow(tt,xx,omega,param),[0 60],zeros(4,1));
    q = q_sim(end,:).';
    [dq] = model.calculateFlow(60,q,omega,param);
    if max(abs(dq)) > 0.1
        fprintf('%f, ',dq)
        fprintf('\n')
    end
    % Thermodynamics(6)
    for i=1:4
        dT_w(i) = (param.thermo.C_w * param.thermo.rho_w * q(i) * (T_c - T_w(i)) - param.thermo.B(i) * (T_w(i) - T_a(i)))/(param.thermo.C_w * param.thermo.rho_w * param.thermo.V_w(i));
        dT_a(i) = (param.thermo.C_a * param.thermo.rho_a * Q(i) * (T_amb - T_a(i)) + param.thermo.B(i) * (T_w(i) - T_a(i)))/(param.thermo.C_a * param.thermo.rho_a * param.thermo.V_a(i));
    end

    dx = [dT_w; dT_a; dT];
end

