function spaneff = Getspaneff(N,AR,S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Set span efficiency including effects of sideslip due to flying
% in a circular path
UEFC = GetUEFC;

R = UEFC.R;
b = sqrt(S*AR);

rbar = 0.5*b/R/N;

dihedral = pi/180*UEFC.dihedral;

CL = GetCL(N, AR, S);

beta = CL/dihedral*(1+4/AR)/(2*pi)*rbar;

spaneff = UEFC.e0*(1-0.5*rbar^2)*(cos(beta))^2;
