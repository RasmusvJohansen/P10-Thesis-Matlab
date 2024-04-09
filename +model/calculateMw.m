function [Mw] = calculateMw(param)
%CALCULATEMW Summary of this function goes here
%   Detailed explanation goes here
s = tf('s')
Am = [param.model.A+param.model.B*param.ctrl.K*param.model.Cy, zeros(12,4);
      param.model.Bwd*param.ctrl.K*param.model.Cy, param.model.Awd];

Bm = [param.model.B;
      zeros(4,4)]; 

Cm = [param.model.Dwd*param.ctrl.K*param.model.Cy, param.model.Cwd];

Dm = [zeros(4,4)];

Mw = Cm/(s*eye(16)-Am)*Bm+Dm;
end

