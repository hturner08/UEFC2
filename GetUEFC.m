function [UEFC] = GetUEFC()

% EXCEPT FOR CHANGING THE CONSTANTS AS YOU LOOK AT DIFFERENT
% TAPER, DBMAX, ETC YOU SHOULD NOT NEED TO CHANGE THIS FILE.

UEFC.CLmax = 2.0; % maximum CL wing will be allowed to fly at
UEFC.dbmax = 0.2; % tip displacement bending constraint
UEFC.lambda = 0.5; % taper ratio
UEFC.dihedral = 10.0; % Wing dihedral (degrees)
UEFC.tau = 0.08; % thickness-to-chord ratio
UEFC.e0 = 1.0; % Span efficiency for straight level flight
UEFC.rho = 1.225; % air density kg/m^3
UEFC.mu = 1.789E-5; % dynamic viscosity of air N s/m^2
UEFC.g = 9.81; % gravity, m/s^2
UEFC.R = 12.5; % radius of loop in meters
UEFC.rhofoam = 32.0; % kg/m^3
UEFC.Efoam = 19.3E6; % Pa

