% Defines the parameters of the system into a struct called param

    


% Thermodynamics parameters
    % Heat transfer constant between water and air in AHU [J/K]
    param.thermo.B = [24, 14, 12, 20] * 10^3;
    % Specific heat coefficient for water [J/kg K]
    param.thermo.C_w = 4183;
    % Specific heat coefficient for air [J/kg K]
    param.thermo.C_a = 728;
    % mass density of water [kg/m3]
    param.thermo.rho_w = 997;                        
    % mass density of air [kg/m3]
    param.thermo.rho_a = 1.225;
    % Water volume in AHU [m3]
    param.thermo.V_w = [230.00, 134.00, 115.00, 192.00] * 10^-3;
    % Air volume in AHU [m3]
    param.thermo.V_a = [9.00 5.20 4.50 7.50];
    % Nominal ambient temperature [K]
    param.thermo.T_A = 308.15;
    % Nominal supply water temperature [K]
    param.thermo.T_c = 283.15;

% Control parameters
    % Nominal air temperature reference [K]
    param.ctrl.T_ref = 293.15;
    % Nominal water flow through AHU [m3/s]
    param.ctrl.q = [5.75, 3.36, 2.89, 4.81] * 10^-3;
    % Nominal air flow through AHU [m3/s]
    % param.ctrl.Q  = [5.98 3.49 2.99 4.98];
    % QNorm from their matlab but not the paper
    param.ctrl.Q  = [8.9706    5.2329    4.4853    7.4755];
   


% Hydraulic parameters
    % Hydraulic resistance of the chilled water source [pa/(m3 /s)]
    param.pipe.R_c = 8.85 * 10^3;
    % Hydraulic resistance of supply/return pipe to branch [pa/(m3/s)]
    param.pipe.R = [8.85, 20.45, 42.23, 108.26] * 10^3;
    % Hydraulic resistance of branch [pa/(m3/s)]
    param.pipe.r = [151.23, 442.59, 599.11, 216.51] * 10^3;
    % Pump constant [.]
    param.pump.a = [453.69, 2212.96, 4193.79, 1948.61] * 10^3;
    % Pump constant [.]
    param.pump.b = [30.00, 50.00, 70.00, 90.00];

% Model parameters
    
    % Bilinear model 
      param.model.F_cali=0; 
      param.model.G_cali=0;  
      param.model.M=0; %contains all the individual matrices M_i in a cell array used in the functions H_xi and H_q
      param.model.Epsilon=0;
      % Linearisation 
    % Linearised state matrix A of decoupled and coupled hydronics
      param.model.A =0; % This will be set in the code as it depend on the operation point

    % Linearised input matrix B for decoupled hydronics
      param.model.B_Bar =0; % This will be set in the code as it depend on the operation point

    % Linearised input matrix B for decoupled hydronics
      param.model.B= 0; % This will be set in the code as it depend on the operation point
      
    % Output matrix 
      param.model.Cy = eye(12); 

% Weight Parameters
    
    % ww weight (defined as unity gain)
      param.model.Aww = zeros(4);
      param.model.Bww = zeros(4);
      param.model.Cww = zeros(4);
      param.model.Dww = diag([1,1,1,1]);
    % wz weight
      param.model.Awz = zeros(4);
      param.model.Bwz = zeros(4);
      param.model.Cwz = zeros(4);
      param.model.Dwz = diag([1,1,1,1]);

% Dimmensions
    %Number of AHUs obtained from pump parameter
    param.n = length(param.pump.a);  
    

