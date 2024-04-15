function param = calculateWDelta(param,unityGain)
%CALCULATEFILTER calculates the A,B,C,D matrices of the filter and add them
%to param struct 
%   Detailed explanation goes here

%apply Wdelta
if(unityGain == 0)    
    unc =  [param.model.B(1,:)-param.model.B_Bar(1,:);
            param.model.B(4,:)-param.model.B_Bar(4,:);
            param.model.B(7,:)-param.model.B_Bar(7,:);
            param.model.B(10,:)-param.model.B_Bar(10,:)];
    Wdelta = unc*tf([1],[1]);
    WdeltaTF2SS = ss(Wdelta);
    
    WdeltaSS.A = zeros(4);
    WdeltaSS.B = zeros(4);
    WdeltaSS.C = zeros(4);
    WdeltaSS.D = WdeltaTF2SS.D;

%Use a passthrough wdelta
else
% Create a unity gain feedforward filter
    Wdelta = eye(4)*tf([1],[1]);
    WdeltaTF2SS = ss(Wdelta);
    WdeltaSS.A = zeros(4);
    WdeltaSS.B = zeros(4);
    WdeltaSS.C = zeros(4);
    WdeltaSS.D = WdeltaTF2SS.D;

end
param.model.Awd = WdeltaSS.A;
param.model.Bwd = WdeltaSS.B;
param.model.Cwd = WdeltaSS.C;
param.model.Dwd = WdeltaSS.D;

end

