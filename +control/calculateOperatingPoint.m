function [T_wOP, q_OP, alpha, omega_s_OP, omega_OP] = calculateOperatingPoint(T_aOP,i)
%CALCULATEOPERATINGPOINT Summary of this function goes here
%   Detailed explanation goes here
    model.Parameters

    % Thermodynamic operating point
    T_wOP = (-param.thermo.C_a*param.thermo.rho_a * param.ctrl.Q(i) * (param.thermo.T_A - T_aOP) + param.thermo.B(i)*T_aOP)/(param.thermo.B(i));
    q_OP = (param.thermo.B(i)*(T_wOP - T_aOP))/(param.thermo.C_w*param.thermo.rho_w*(param.thermo.T_c - T_wOP));
    
    % Hydraulic operating point
        % Decoupled system, meaning only i'th flow is considered
        alpha = sqrt((param.pump.a(i)+param.pipe.r(i)+param.pipe.R_c + 2*sum(param.pipe.R(1:i)))/param.pump.b(i));
        omega_s_squared = (param.pump.a(1)+param.pipe.r(i))/param.pump.b(i) * param.ctrl.q(i)^2 + param.pipe.R_c/param.pump.b(i)*param.ctrl.q(i)^2;
        for k=1:i  % Loops the pipe sections only considering the i'th pumps flow
            omega_s_squared = omega_s_squared + 2*param.pipe.R(k)/param.pump.b(i)*param.ctrl.q(i)^2;
        end
        omega_s_OP = sqrt(omega_s_squared);

        % Coupled system, meaning every flow is considered
        omega_squared = (param.pump.a(1)+param.pipe.r(i))/param.pump.b(i) * param.ctrl.q(i)^2 + param.pipe.R_c/param.pump.b(i)*sum(param.ctrl.q)^2;
        for k=1:i % Loops the pipe sections summing over the flows of the k'th and remaining pumps
            omega_squared = omega_squared + 2*param.pipe.R(k)/param.pump.b(i)*sum(param.ctrl.q(k:end))^2;
        end
        omega_OP = sqrt(omega_squared);


    
end

