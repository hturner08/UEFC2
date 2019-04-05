function [nV] = Calc_nV(x, AR, S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Return the negative of V because fmincon minimizes a function,
% and we desire to maximize speed.  Minimizing the negative of
% V will maximize V.

nV = -GetV(x, AR, S);














