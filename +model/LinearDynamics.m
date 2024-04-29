function [dx,q] = LinearDynamics(t,x,param,decoupled,step_time)
    % States
    T_w = x(1:4);
    T_a = x(5:8);
    delta = x(9:12);
    if nargin > 4
        if t>step_time % Perform step if time
            T_ref = param.ctrl.T_ref + 1;
        else
            T_ref = param.ctrl.T_ref;
        end
    else
        T_ref = param.ctrl.T_ref;
    end
    % x is [T_w;T_a;delta] but it should be [T_w1;T_a1;delta1;T_w2;T_a2;delta2,...]
    xLin = zeros(12,1);
    for i=1:4
        xLin((i-1)*3+1:i*3,1) = [T_w(i) - param.ctrl.T_wOP(i); T_a(i) - param.ctrl.T_aOP(i); delta(i)];
    end
    q = param.ctrl.K*xLin;

    if(decoupled == 1)
        dx_temp = (param.model.A + param.model.B_Bar * param.ctrl.K)*xLin - blkdiag([0;0;1],[0;0;1],[0;0;1],[0;0;1])*ones(4,1)*(T_ref - param.ctrl.T_ref);
    else
        dx_temp = (param.model.A + param.model.B * param.ctrl.K)*xLin - blkdiag([0;0;1],[0;0;1],[0;0;1],[0;0;1])*ones(4,1)*(T_ref - param.ctrl.T_ref);
    end
    dx = [dx_temp(1);dx_temp(4);dx_temp(7);dx_temp(10);dx_temp(2);dx_temp(5);dx_temp(8);dx_temp(11);dx_temp(3);dx_temp(6);dx_temp(9);dx_temp(12)];
end

