function [V, N, exitflag] = opt_V(AR, S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Get various fixed parameters
UEFC = GetUEFC;

% Initial guess for optimizer
N0  = 1.1;

% Lower and upper bounds
lbN = 1.0001;
ubN = 5.0;

lb = [lbN];
ub = [ubN];

% Optimization variables are stored in x
x0 = [N0];

% Set-up options for fmincon
options = optimoptions(@fmincon);
options.Algorithm = 'sqp';
options.Display = 'off';

[x,fmin,exitflag] = fmincon(@(x)Calc_nV(x,AR,S),x0,[],[],[],[],lb,ub, ...
                   @(x)Calc_constraints(x,AR,S),options);
if (exitflag > 0),
    V = -fmin;
    N = x;
else
    V = 0;
    N = 0;
end









