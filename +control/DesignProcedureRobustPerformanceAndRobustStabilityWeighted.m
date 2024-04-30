function param = DesignProcedureRobustPerformanceAndRobustStability(param,info)
%DESIGNPROCEDUREROBUSTPERFORMANCEANDROBUSTSTABILITY Summary of this function goes here
yalmip('clear')

q = cell(param.n);
y = cell(param.n);
%fill in p^-1
for i = 1:param.n
    q{i} = sdpvar(3);
    y{i} = sdpvar(1,3,'full');
end
% fill in the rest of p_bar
for i = 1:4
    q{param.n+i} = sdpvar(4);
    y{param.n+i} = sdpvar(4);
end
Q_hat = blkdiag(q{:})
Y_hat = blkdiag(y{:})
gammaPerformance = sdpvar(1);
gammaDesiredPerformance =sdpvar(1);

%constraints
constriants = [];

%Define the lyapunaov constraint for the decoupled system 
    Lyap = param.model.A * Q_hat(1:param.n*3,1:param.n*3) + Q_hat(1:param.n*3,1:param.n*3) * param.model.A.' + param.model.B_Bar * Y_hat(1:param.n,1:param.n*3) + Y_hat(1:param.n,1:param.n*3).' * param.model.B_Bar.';
    constriants = [constriants, Lyap <= 0 Q_hat(1:param.n*3,1:param.n*3) >= 0];

%Define the parameters for desired response
Cz = blkdiag([0 1 0],[0 1 0],[0 1 0],[0 1 0]);
Bw = -blkdiag([0;0;1],[0;0;1],[0;0;1],[0;0;1]);
Dz = eye(4);

A_hat = [param.model.A, Bw*param.model.Cww, zeros(12,4),zeros(12,4),zeros(12,4);
         zeros(4,12),param.model.Aww, zeros(4),zeros(4),zeros(4);
         param.model.Bwz*Cz, -param.model.Bwz*param.model.DWref*param.model.Cww, param.model.Awz, -param.model.Bwz*param.model.CWref,zeros(4);
         zeros(4,12),param.model.BWref*param.model.Cww, zeros(4),param.model.AWref,zeros(4);
         zeros(4,12), zeros(4),zeros(4),zeros(4),param.model.AWu]; 

B_hat = [param.model.B, zeros(12,4),zeros(12,4),zeros(12,4),zeros(12,4);
         zeros(4,20);
         zeros(4,20);
         zeros(4,20);
         param.model.BWu,zeros(4,16)];

B_bar = [Bw*param.model.Dww, param.model.B;
         param.model.Bww, zeros(4);
         -param.model.Bwz*param.model.DWref*param.model.Dww, zeros(4);
         param.model.BWref*param.model.Dww,zeros(4);
         zeros(4,8)];

C1 = [param.model.Dwz*Cz, -param.model.Dwz*param.model.DWref*param.model.Cww, param.model.Cwz,-param.model.Dwz*param.model.CWref,zeros(4);
      zeros(4,12),zeros(4),zeros(4),zeros(4),param.model.CWu;
      zeros(4,28)];


C2 = [zeros(4,20);
      param.model.DWu,zeros(4,16);
      param.model.Dwd, zeros(4,16)];

D_Bar = [-param.model.Dwz*param.model.DWref*param.model.Dww,zeros(4);
        zeros(4,8);
        zeros(4,8)];


DesiredPerformance = [Q_hat*A_hat.'+A_hat*Q_hat+Y_hat.'*B_hat.'+B_hat*Y_hat, B_bar, Q_hat*C1.'+Y_hat.'*C2.';
                     B_bar.', -gammaDesiredPerformance*eye(8),D_Bar.';
                     C1*Q_hat+C2*Y_hat, D_Bar, -gammaDesiredPerformance*eye(12)];

constriants = [constriants, DesiredPerformance <= 0,gammaDesiredPerformance >= 0, Q_hat >=0];

% if info 
%         sol = optimize(constraints, [gammap], options)
%         check(constraints)
%     else
options = sdpsettings('verbose',0,'solver','mosek');    
if(info)
    sol = optimize(constriants,[gammaDesiredPerformance],options)
    check(constriants)
else
    sol = optimize(constriants,[gammaDesiredPerformance],options);

end

disp("gamma per")
value(gammaPerformance)
disp("gamma desired per")
value(gammaDesiredPerformance)

Kfull = value(Y_hat)/value(Q_hat)    
k = Kfull(1:4,1:12)
Ks = zeros(1,3,param.n);
  for i=1:param.n
        Ks(:,:,i) = k(i,3*(i-1)+1:3*i);
  end
param.ctrl.Ks = Ks;
param.ctrl.K = k;

if(info)
    disp("Eig of Q_hat")
    eig(value(Q_hat))
end
end

