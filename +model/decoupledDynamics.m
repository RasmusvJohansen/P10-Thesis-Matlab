function [dx,q] = decoupledDynamics(t,x,i,param,step_time)    
    % States
    T_w = x(1);
    T_a = x(2);
    
    % Inputs
    T_c = param.thermo.T_c;
    T_amb = param.thermo.T_A;
    Q = param.ctrl.Q(i);

    % Controller
    if nargin > 4
        if t>step_time % Perform step if time
            T_ref = param.ctrl.T_ref + 1;
        else
            T_ref = param.ctrl.T_ref;
        end
    else
        T_ref = param.ctrl.T_ref;
    end
    dT = T_a - T_ref;
    % Apply linearisation
    x = [T_w - param.ctrl.T_wOP(i); T_a - param.ctrl.T_aOP(i); x(3)];
    % (14) with the linearisation offset replaced by omega_s_OP calculated
    % for the OP
    omega = param.ctrl.Ks(:,:,i) * x + param.ctrl.omega_s_OP(i);

    % Hydraulic network (13)
    q = omega/param.ctrl.alpha(i);
    q = max(q,0);

    % Thermodynamics(6)
    dT_w = (param.thermo.C_w * param.thermo.rho_w * q * (T_c - T_w) - param.thermo.B(i) * (T_w - T_a))/(param.thermo.C_w * param.thermo.rho_w * param.thermo.V_w(i));
    dT_a = (param.thermo.C_a * param.thermo.rho_a * Q * (T_amb - T_a) + param.thermo.B(i) * (T_w - T_a))/(param.thermo.C_a * param.thermo.rho_a * param.thermo.V_a(i));
    
    dx = [dT_w; dT_a; dT];
end

