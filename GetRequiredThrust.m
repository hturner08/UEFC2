function T = GetRequiredThrust(N, AR, S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Calculate the required thrust from UEFC parameters and N, AR, S
UEFC = GetUEFC;

rho = UEFC.rho;
g = UEFC.g;
R = UEFC.R;

W = GetWeight(AR, S);

% Calculate CL (in order to find CD)
CL = GetCL(N, AR, S);

% Calculate CD
CD = GetCD(N, AR, S);

% Calculate q (use GetV)
V = GetV(N, AR, S);
q = 0.5*rho*V^2;

% Calculate required thrust from CD, q, S
T = CD*q*S;



