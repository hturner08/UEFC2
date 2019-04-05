function Wpay = GetWpay(AR,S)

% Return payload weight in N

UEFC = GetUEFC;

% WARNING: the following payload mass is just the burrito.  You
% will need to add any additional mass due to any packaging, etc
% that you decide to use.

mpay_g = 500; % Mass in grams


Wpay = UEFC.g*mpay_g/1000; 
