function Wwing = GetWingWeight(AR, S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Calculate the wing weight from UEFC parameters and S and AR
UEFC = GetUEFC;

tau = UEFC.tau;
lambda = UEFC.lambda;
g = UEFC.g;
rhofoam = UEFC.rhofoam;
dihedral = pi/180*UEFC.dihedral;

Afac = 0.66; % Area of airfoil assumed to be Afac*tau*c^2

Wwing = (4/3)*Afac*rhofoam*g*tau/cos(dihedral)*S.^(1.5).*AR.^(-0.5)*(lambda^2+lambda+1)/(lambda+1)^2;

