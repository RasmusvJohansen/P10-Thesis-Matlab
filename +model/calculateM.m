function [M] = calculateM(param)
%Calculate the interconnection M.
s = tf('s')

M = param.ctrl.K * param.model.Cy/(s * eye(12) - (param.model.A + param.model.B * param.ctrl.K * param.model.Cy) )*param.model.B; 
end

