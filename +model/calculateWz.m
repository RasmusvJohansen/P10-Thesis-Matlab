function param = calculateWz(param,filter,feedforward)
%CALCULATEWZ Summary of this function goes here
%   Detailed explanation goes here Filter can either be a tf of the filter
%   which should be applied on Z. If feedforward is wanted the feedforward
%   parameter can be set to 1 to force feedforward if 0 the supplied filter
%   is used. 
    s=tf('s');
    if(feedforward == 1)
     param.model.Awz = zeros(4);
     param.model.Bwz = zeros(4);
     param.model.Cwz = zeros(4);
     param.model.Dwz = 1*eye(4);
    else
    
    WZ = filter;
    WZSS = ss(WZ);
    param.model.Awz = WZSS.A;
    param.model.Bwz = WZSS.B;
    param.model.Cwz = WZSS.C;
    param.model.Dwz = WZSS.D;
    end
end

