function param= calculateWw(param,filter,feedforward)
%CALCULATEWZ Summary of this function goes here
%   Detailed explanation goes here Filter can either be a tf of the filter
%   which should be applied on w. If feedforward is wanted the feedforward
%   parameter can be set to 1 to force feedforward if 0 the supplied filter
%   is used. 
    s=tf('s');
    if(feedforward == 1)
     param.model.Aww = zeros(4);
     param.model.Bww = zeros(4);
     param.model.Cww = zeros(4);
     param.model.Dww = eye(4);
    else
    
    Ww = filter;
    WwSS = ss(Ww);
    param.model.Aww = WwSS.A;
    param.model.Bww = WwSS.B;
    param.model.Cww = WwSS.C;
    param.model.Dww = WwSS.D;
    end
end

