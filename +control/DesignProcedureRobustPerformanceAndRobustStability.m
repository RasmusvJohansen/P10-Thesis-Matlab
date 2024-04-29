function param = DesignProcedureRobustPerformanceAndRobustStability(param,info)
%DESIGNPROCEDUREROBUSTPERFORMANCEANDROBUSTSTABILITY Summary of this function goes here
%   Detailed explanation goes here

s = tf('s');
yalmip('clear')

q = cell(param.n);
y = cell(param.n);
%fill in p^-1
for i = 1:param.n
    q{i} = sdpvar(3);
    y{i} = sdpvar(1,3,'full');
end
% fill in the rest of p_bar
for i = 1:3
    q{param.n+i} = sdpvar(4);
    y{param.n+i} = sdpvar(4);
end
Q_hat = blkdiag(q{:});
Y_hat = blkdiag(y{:});
gammaRS = sdpvar(1); %(1);
gammaPerformance = sdpvar(1);

%constraints
constriants = [];

%Define the lyapunaov constraint for the decoupled system 
Lyap = param.model.A * Q_hat(1:param.n*3,1:param.n*3) + Q_hat(1:param.n*3,1:param.n*3) * param.model.A.' + param.model.B_Bar * Y_hat(1:param.n,1:param.n*3) + Y_hat(1:param.n,1:param.n*3).' * param.model.B_Bar.';
constriants = [constriants, Lyap <= 0 Q_hat(1:param.n*3,1:param.n*3) >= 0];

%Define the paramerts for robust stability
% Structure the matrices
%     A_split = [param.model.A, zeros(12,4); zeros(4,12), param.model.Awd];
%     B_split = [param.model.B, param.model.B * param.model.Cwd; zeros(4,4), zeros(4,4)];
%     % K_split = [K, zeros(4); zeros(4), eye(4)]
%     B_mw = [param.model.B * param.model.Dwd; param.model.Bwd];
%     C_split = [eye(4), zeros(4)];
%     D_mw = zeros(4);
% 
% 
% %Define the Robust stability constraints
% RobustStability = [A_split * Q_hat(1:param.n*3+4,1:param.n*3+4) + Q_hat(1:param.n*3+4,1:param.n*3+4) * A_split.' + B_split * Y_hat(1:param.n+4,1:param.n*3+4) + Y_hat(1:param.n+4,1:param.n*3+4).' * B_split.', B_mw, Y_hat(1:param.n+4,1:param.n*3+4).' * C_split.';
%                     B_mw.', -gammaRS * eye(4), D_mw;
%                     C_split * Y_hat(1:param.n+4,1:param.n*3+4), D_mw.', -gammaRS * eye(4)];
% constriants = [constriants, RobustStability <= 0,Q_hat(1:param.n*3+4,1:param.n*3+4) >=0,gammaRS >= 0, gammaRS <=1];


%Define paramerts for robust performance
Cz = -blkdiag([0 1 0],[0 1 0],[0 1 0],[0 1 0]);
Bw = -blkdiag([0;0;1],[0;0;1],[0;0;1],[0;0;1]);
Dz = eye(4);
A_hat = [param.model.A, zeros(12,4) Bw*param.model.Cww, zeros(12,4);
         zeros(4,12),param.model.Awd, zeros(4),zeros(4);
         zeros(4,12),zeros(4),param.model.Aww, zeros(4);
         param.model.Bwz*Cz, zeros(4),param.model.Bwz*Dz*param.model.Cww, param.model.Awz]; 

B_hat = [param.model.B, zeros(12,4),zeros(12,4),zeros(12,4);
         param.model.Bwd, zeros(4),zeros(4),zeros(4);
         zeros(4,16)
         zeros(4,16)];

B_bar = [param.model.B, Bw*param.model.Dww;
         zeros(4,4),zeros(4,4);
         zeros(4,4),param.model.Bww;
         zeros(4,4),param.model.Bwz*Dz*param.model.Dww];

C1 = [zeros(4,12), param.model.Cwd,zeros(4),zeros(4);
      param.model.Dwz*Cz , zeros(4), param.model.Dwz*Dz*param.model.Cww, param.model.Cwz];
C2 = [param.model.Dwd, zeros(4),zeros(4),zeros(4);
      zeros(4),zeros(4),zeros(4),zeros(4)];

D_Bar = [zeros(4),zeros(4);
         zeros(4), param.model.Bwz*Dz*param.model.Dww];

%Define the robust performance constraint

RoubstPerformance = [Q_hat*A_hat.'+A_hat*Q_hat+Y_hat.'*B_hat.'+B_hat*Y_hat, B_bar, Q_hat*C1.'+Y_hat.'*C2.';
                     B_bar.', -gammaPerformance*eye(8),D_Bar.';
                     C1*Q_hat+C2*Y_hat, D_Bar, -gammaPerformance*eye(8)];

constriants = [constriants, RoubstPerformance <= 0, Q_hat >=0,gammaPerformance >= 0];

% if info 
%         sol = optimize(constraints, [gammap], options)
%         check(constraints)
%     else
options = sdpsettings('verbose',0,'solver','mosek');    
if(info)
    sol = optimize(constriants, [gammaPerformance], options)
    check(constriants)
else
    sol = optimize(constriants, [gammaPerformance], options);

end

disp("gamma per")
value(gammaPerformance)

disp("gamma robust stab")
value(gammaRS)

Kfull = value(Y_hat)/value(Q_hat)    
k = Kfull(1:4,1:12)
Ks = zeros(1,3,param.n);
  for i=1:param.n
        Ks(:,:,i) = k(i,3*(i-1)+1:3*i);
  end
param.ctrl.Ks = Ks;
param.ctrl.K = k;

if(info)
    disp("Eig of Robust performance")
    eig(value(RoubstPerformance))
    disp("Eig of Q_hat")
    eig(value(Q_hat))
end
end

