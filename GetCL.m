function CL = GetCL(N, AR, S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Calculate the lift coefficient from UEFC parameters and N, AR, S
UEFC = GetUEFC;

rho = UEFC.rho;
g = UEFC.g;
R = UEFC.R;

V = GetV(N, AR, S);
W = GetWeight(AR, S);
CL = N*W/(0.5*rho*S*V^2);




