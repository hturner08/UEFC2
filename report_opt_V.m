function [V, N, exitflag] = report_opt_V(AR, S)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Get various fixed parameters
UEFC = GetUEFC;

[V, N, exitflag] = opt_V(AR, S);

if (exitflag > 0),
    fprintf('AR      = %5.3f\n',AR);
    fprintf('S       = %5.3f sq. m\n',S);
    fprintf('b       = %5.3f m\n',sqrt(AR*S));
    fprintf('cbar    = %5.3f m\n',sqrt(S/AR));
    fprintf('ct      = %5.3f m\n',2*UEFC.lambda/(1+UEFC.lambda)*sqrt(S/AR));
    fprintf('cr      = %5.3f m\n',2/(1+UEFC.lambda)*sqrt(S/AR));
    fprintf('lambda  = %5.3f\n',UEFC.lambda);
    fprintf('tau     = %5.3f\n',UEFC.tau);
    fprintf('eps     = %5.3f\n',Getepsilon(UEFC.tau));
    fprintf('V       = %5.3f m/s\n',V);
    fprintf('N       = %5.3f\n',N);
    [W, Wfuse, Wwing, Wpay] = GetWeight(AR, S);
    fprintf('W/g     = %4.0f g\n',W/UEFC.g*1000);
    fprintf('Wfuse/g = %4.0f g\n',Wfuse/UEFC.g*1000);
    fprintf('Wwing/g = %4.0f g\n',Wwing/UEFC.g*1000);
    fprintf('Wpay/g  = %4.0f g\n',Wpay/UEFC.g*1000);
    CL = GetCL(N, AR, S);
    fprintf('CL      = %5.3f\n',CL);
    fprintf('CLmax   = %5.3f\n',UEFC.CLmax);
    [CD, CDfuse, CDp, CDi, CDpay] = GetCD(N, AR, S);
    fprintf('CD      = %5.3f\n',CD);
    fprintf('CDfuse  = %5.3f\n',CDfuse);
    fprintf('CDp     = %5.3f\n',CDp);
    fprintf('CDi     = %5.3f\n',CDi);
    fprintf('CDpay   = %5.3f\n',CDpay);
    fprintf('e0      = %5.3f\n',UEFC.e0);
    fprintf('e       = %5.3f\n',Getspaneff(N,AR,S));
    T = GetRequiredThrust(N, AR, S);    
    fprintf('T       = %5.3f N\n',T);
    Tmax = GetMaxThrust(V);    
    fprintf('Tmax    = %5.3f N\n',Tmax);
    db = Getdb(N, AR, S);
    fprintf('d/b     = %5.3f\n',db);
    fprintf('d/bmax  = %5.3f\n',UEFC.dbmax);
else
    fprintf('Error in opt_V: exitflag = %i\n',exitflag);
end











