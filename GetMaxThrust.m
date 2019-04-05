function Tmax = GetMaxThrust(V)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

UEFC = GetUEFC;

ct0 =  0.2093;
ct1 = -0.2484;
ct2 = -0.1386;

Tmax_static = 2; % Maximum thrust desired at static conditions
Rprop = 0.1016; % Propellor radius (m)
Aprop = pi*Rprop^2; % Propellor disk area
Omega = sqrt(Tmax_static/(0.5*UEFC.rho*Rprop^2*Aprop*ct0));

Lambda = V/(Omega*Rprop); % Advance ratio

CT = ct0 + ct1*Lambda + ct2*Lambda.^2; % Thrust coefficient

Tmax = 0.5*UEFC.rho*(Omega*Rprop)^2*Aprop*CT;

