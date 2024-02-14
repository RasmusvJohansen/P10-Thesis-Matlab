function [dx,q] = decoupledDynamics(t,x,u,i,K,param,step_time)
% DECOUPLEDDYNAMICS Summary of this function goes here
%   Detailed explanation goes here
    
    % States
    T_w = x(1);
    T_a = x(2);
    
    % Inputs
    T_c = u(1);
    T_amb = u(2);
    Q = u(3);

    % Controller(14)
    if t>step_time % Perform step if time
        T_ref = param.ctrl.T_ref+1;
    else
        T_ref = param.ctrl.T_ref;
    end
    dT = T_a - T_ref;
    % Apply linearisation
    x = [T_w - param.ctrl.T_wOP(i); T_a - param.ctrl.T_aOP(i); x(3)];
    omega = param.ctrl.alpha(i) * K(:,:,i) * x + param.ctrl.alpha(i) * param.ctrl.q(i);

    % Hydraulic network (13)
    q = omega/param.ctrl.alpha(i);
    q = max(q,0);

    % Thermodynamics(6)
    dT_w = (param.thermo.C_w * param.thermo.rho_w * q * (T_c - T_w) - param.thermo.B(i) * (T_w - T_a))/(param.thermo.C_w * param.thermo.rho_w * param.thermo.V_w(i));
    dT_a = (param.thermo.C_a * param.thermo.rho_a * Q * (T_amb - T_a) + param.thermo.B(i) * (T_w - T_a))/(param.thermo.C_a * param.thermo.rho_a * param.thermo.V_a(i));
    

    dx = [dT_w; dT_a; dT];
end

