function param = calculateOperatingPoint(param)
% Calculates the operating point for both the hydraulics, thermodynamics
% and control parameters. First the water temperature is calculated, and
% then used in the calculation of the flow at the operating point. Next the
% flow is used to calculate the control parameter alpha and the two angular
% velocity offset for the pumps at the linearisation point.
    
    fprintf("Calculate Operating Point for T_a = %f [K]",param.ctrl.T_ref)
    T_aOP = param.ctrl.T_ref;
    for i=1:param.n
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


        param.ctrl.alpha(i) = alpha;
        param.ctrl.T_wOP(i) = T_wOP;
        param.ctrl.T_aOP(i) = T_aOP;
        param.ctrl.omega_OP(i) = omega_OP;
        param.ctrl.omega_s_OP(i) = omega_s_OP;
        param.ctrl.q_OP(i) = q_OP;
        param.ctrl.xi_OP(3*(i-1)+1:3*i,1) = [T_wOP; param.ctrl.T_aOP(i); 0];
    end
    param.ctrl.Lambda_Bar = diag(param.ctrl.alpha);

end

