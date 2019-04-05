function [V] = GetV(N, AR, S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Calculate speed from N, g, R

UEFC = GetUEFC;
g = UEFC.g;
R = UEFC.R;

V = sqrt(R*g*sqrt(N^2-1));













