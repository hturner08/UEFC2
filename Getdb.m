function db = Getdb(N,AR,S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Calculate the tip d/b from UEFC parameters and S, AR, N
UEFC = GetUEFC;
lambda = UEFC.lambda;
tau = UEFC.tau;
epsilon = Getepsilon(tau);
Wfuse = GetWfuse(AR,S);
Wpay = GetWpay(AR,S);

db = 0.018*N*(Wfuse+Wpay)/(UEFC.Efoam*tau*(tau^2+0.7*epsilon^2))* ...
         (1+lambda)^3*(1+2*lambda)*AR^3/S;

