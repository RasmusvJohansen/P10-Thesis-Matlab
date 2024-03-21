function param = calculateModel(param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    disp("Calculate bilinear model matrices")
    F = cell(1,param.n);
    G = cell(1,param.n);
    M = cell(1,param.n);
    E =[]; % defines E to be empty
    for i = 1:param.n
        F{i} = [-(param.thermo.B(i))/(param.thermo.C_w*param.thermo.rho_w*param.thermo.V_w(i)), (param.thermo.B(i))/(param.thermo.C_w*param.thermo.rho_w*param.thermo.V_w(i)), 0; 
                 (param.thermo.B(i))/(param.thermo.C_a*param.thermo.rho_a*param.thermo.V_a(i)), -((param.ctrl.Q(i))/(param.thermo.V_a(i)) + (param.thermo.B(i))/(param.thermo.C_a*param.thermo.rho_a*param.thermo.V_a(i))), 0; 
                 0, 1, 0,];
        G{i} = [param.thermo.T_c/param.thermo.V_w(i); 0; 0];
        M{i} = [-1/param.thermo.V_w(i), 0, 0; 0, 0, 0; 0, 0,0];
        E =[E; [0; (param.ctrl.Q(i)*param.thermo.T_A)/(param.thermo.V_a(i)); -param.ctrl.T_ref]];
    end
    
    param.model.F_cali = blkdiag(F{:});
    param.model.G_cali = blkdiag(G{:});
    param.model.M = M;
    param.model.Epsilon = E;
    %H(xi) is found as a function under +model as the function "H_xi.m"
    disp("Calculate linearisation decoupled")
    
    param.model.A = param.model.F_cali + model.H_q(param.ctrl.q_OP,param);
    param.model.B_Bar = (param.model.G_cali + model.H_xi(param.ctrl.xi_OP,param)) / param.ctrl.Lambda_Bar;
    
    disp("Calculate linearisation coupled")
    % Lånt fra John
    n = param.n;
    Lambda = zeros(n,n,n);
    Psi    = zeros(n,n,n);
    Gamma  = zeros(n,n,n);
    S      = zeros(n,n,n);
    
    for i = 1:n
        % Define matrices to use in Gamma
	    [ G2, G3, G4 ] = deal( zeros(n,n) );
	    G1          = ones(n,n);
	    G2(2:4,2:4) = ones(3,3);
	    G3(3:4,3:4) = ones(2,2);
	    G4(4,4)     = 1;
        
        % diag([(i==1) (i==2) (i==3) (i==4)]) create a matrix with 1 in the
        % i'th diagonal entry
	    Lambda(:,:,i) = ( param.pipe.r(i) + param.pump.a(i) )/param.pump.b(i)*diag([(i==1) (i==2) (i==3) (i==4)]);
        % Psi is a scalar multiplied onto a 1 matrix
	    Psi(:,:,i)    = param.pipe.R_c/param.pump.b(i)*ones(n,n);
        % Gamma is a sum of matrices reducing in size, meaning the first 
        % has an entry in ever index, while the second has in every except
        % the first row and column, and so on...
	    Gamma(:,:,i) =  2*param.pipe.R(1)/param.pump.b(i)*G1*(i>0) + ...                    % if i>0 apply this matrix
									    2*param.pipe.R(2)/param.pump.b(i)*G2*(i>1) + ...    % if i>1 apply this matrix
									    2*param.pipe.R(3)/param.pump.b(i)*G3*(i>2) + ...    % if i>2 apply this matrix
									    2*param.pipe.R(4)/param.pump.b(i)*G4*(i>3);         % if i>3 apply this matrix
	    
	    S(:,:,i) = Lambda(:,:,i) + Psi(:,:,i) + Gamma(:,:,i);
    end
    
    % Compute the partial derivative of f
    for ii=1:4
      % Compute rho
      rho(ii) = 1/(2*sqrt(param.ctrl.q_OP(:)'*S(:,:,ii)*param.ctrl.q_OP(:)));
      % Compute Df_q
      Df_q(:,ii) = rho(ii)*S(:,:,ii)*param.ctrl.q_OP(:);
    end
    Dg = inv(Df_q);
    param.model.Dg = Dg;
    param.model.B = (model.H_xi(param.ctrl.xi_OP,param)+param.model.G_cali)*Dg;
    
  
    disp("Testing controlability")
    Ta_ref_singularity= ((param.thermo.T_c*param.thermo.B+param.thermo.C_a*param.thermo.rho_a*param.ctrl.Q*param.thermo.T_A)./(param.thermo.B+param.thermo.C_a*param.thermo.rho_a*param.ctrl.Q));
    if(isequal(round((param.ctrl.T_ref*ones(1,param.n)-Ta_ref_singularity),4),zeros(1,param.n)))
       error("The current selected refference causes a singularity choose another one")
    end
    
    disp("Testing T_ref")
    if(param.ctrl.T_wOP < ones(1,param.n)*param.thermo.T_c)
        error("The choosen T_ref results in operationg point of the water outlet which is smaller then T_c, chose a higher T_ref")
    end


end