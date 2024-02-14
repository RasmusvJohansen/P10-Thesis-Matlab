function dq = calculateFlow(t,q,omega,param)
    % ODE function that calculates the flow of the individual pumps based 
    % on the inputted omega. The flow is calculated based on the pressure
    % drop of the hydraulic network defined by (2) and (3). This is then
    % combined using KVL, which means that if the pressure drop of the pump
    % and pipes aren't equal, then the pressure must be wrong, and hence a
    % change in pressure must be applied.

    % Pump
    % Hvorfor max? Tog det fra carsten og john. det virker ikke uden, men kan ikke se ideen bag
    deltaP_pump = -param.pump.a.' .* abs(q).*q + param.pump.b.' .* max(omega,0).^2;
    deltaP_pipe = zeros(4,1);
    for i=1:4
        % WAHE and Chiller
        deltaP_pipe(i) = param.pipe.r(i) * abs(q(i))*q(i) + param.pipe.R_c*abs(sum(q))*sum(q);
        % Pipes
        for k=1:i
            deltaP_pipe(i) = deltaP_pipe(i) + 2*param.pipe.R(k) * abs(sum(q(k:end))) * sum(q(k:end));
        end
    end
    
    % KVL
    % Same her, skal have dette tjek ellers fejler ode
    dq = zeros(4,1);
    for k=1:4
        if q(k) >= 0
            dq(k) = deltaP_pump(k) - deltaP_pipe(k);
        else
            dq(k) = -q(k);
        end
    end
end

