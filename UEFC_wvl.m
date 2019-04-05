function [e, maxcl] = UEFC_wvl(AR,S,agr,agt,CL,iplot)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Inputs to UEFC_wvl:
%
% S: wing area
% AR: wing aspect ratio
% agr: geometric twist angle at root (deg)
% agt: geometric twist angle at tip (deg)
% CL: desired wing CL
% iplot: plots lift and cl distribution if true

% Outputs from UEFC_wvl:
%
% e: span efficiency (then CDi = CL^2/(pi*AR*e)
% maxcl: maximum sectional lift coefficient

%=====================================================================
%
% The inputs to UEFC_wvl are first converted into:
%
%    x   location of leading edge point
%    y   spanwise location  
%    z
%    c     chord length (along x)
%    a0    airfoil zero-lift angle (deg), negative for positive camber
%    tw    section twist angle (deg), positive nose up
%
%  Note: Only  tw-a0  (the net aerodynamic twist) is significant
%

UEFC = GetUEFC;
lambda = UEFC.lambda;
dihedral = UEFC.dihedral;

b = sqrt(S*AR);
cbar = S/b;
cr = 2*cbar/(1+lambda);
ct = cr*lambda;

% Note: we are setting a0r and a0t to -6 degrees.
% These are reasonable values for the angle of attack at which the lift
% will be zero for all of the airfoils (determined by fitting the data
% at cl of about 1 at representative Reynolds number using XFOIL). 
a0r = -6; a0t = -6;

% Geometry is set for x/c=0.25 to be at x=0
%         x         y     z                           c     tw    a0 
geom = [ -0.25*cr   0     0                           cr    agr   a0r   ; 
         -0.25*ct   b/2   b/2*tan(dihedral*pi/180)    ct    agt   a0t  ];

%-------------------------------------------------------------------------

% Set reference values of area, span, and chord
Sref = S;
bref = b;
cref = cbar;

%=================================================================================
%  FLOW CONDITION SETUP 
%
%  Length of alspec vector indicates number of flow conditions
%
%  A specified quantity is imposed only if its flag = 1, 
%     otherwise it is calculated as part of solution

alspec = 0.0;    % alpha
CLspec = CL;     % CL   
bespec = 0.0;    % sideslip angle
pbspec = 0.0;    % normalized roll rate = roll_rate*bref/(2V)
rbspec = 0.0;    % normalized yaw  rate =  yaw_rate*bref/(2V) = cos(bank)*bref/(2R_turn)
Crspec = 0.0;    % roll moment coefficient
Cnspec = 0.0;    % yaw  moment coefficient  

% Flags indicating which quantity above is to be specified at each flow condition
% Must have exactly 4 flags set for each flow condition
ialspec = 0;
iCLspec = 1;
ibespec = 1;
ipbspec = 1;
irbspec = 1;
iCrspec = 0;      
iCnspec = 0;      

%=================================================================================
% NUMERICAL PROBLEM SETUP

% number of horseshoe vortices on half of wing
N = 25;

%ispace = 1;   % uniform spacing
ispace = 2;   % cosine spacing

% max number of iterations which wvl is allowed to perform
itmax = 20;  

% convergence tolerance
toler = 1.0e-6;

%=================================================================================

% Number of geometry sections
ns = size(geom,1); 

% Number of flow conditions
nflow = 1;

% Convert flow angles to radians
alspec = alspec*(pi/180);
bespec = bespec*(pi/180);

geom(:,1:4) = geom(:,1:4)*2/bref;    % Normalize x,y,z,c to unit reference semispan
geom(:,5:6) = geom(:,5:6)*(pi/180);  % Convert twist,a0L angles to radians

%===================================================================
% perform Weissinger VL calculation

  [yv,zv,cl,ccl,vi,wi,alpha,beta,pbar,rbar,CL,CDi,Cr,Cn,Cb] = wvl(geom,N,ispace,Sref/(bref/2)^2,2,cref/(bref/2),itmax,toler, ...
 alspec, bespec, pbspec, rbspec, CLspec, Crspec, Cnspec, ...
ialspec,ibespec,ipbspec,irbspec,iCLspec,iCrspec,iCnspec);

% span efficiency
e = CL^2 / (pi*AR*CDi);

% max cl along span for flow condition
maxcl = max( cl );

if (iplot),
    
    %===================================================================
    % plot spanwise cl, loading, downwash
    
    subplot(4,1,[1 2],'align');
    plot(yv,cl,'r',yv,ccl,'k','Linewidth',1.5); 
    grid on; 
    legend('$c_l$','$L`/q_{\infty}c_{avg}$','Location','south');
    
    titlea = sprintf( ...
        'CL=%3.2f, alfa=%4.1f (deg), maxcl=%5.3f, CDi=%6.4f, e=%3.2f', ...
        CL,alpha*180/pi,maxcl,CDi,e);
    
    title(titlea,'Fontsize',12);
    
    %===================================================================
    % plot wing
    
    xp = [geom(:,1); geom(end:-1:1,1)+geom(end:-1:1,4)];
    yp = [geom(:,2); geom(end:-1:1,2)];
    zp = [geom(:,3); geom(end:-1:1,3)];
    
    % plot top view
    subplot(4,1,3,'align');
    plot( yp,xp,'b','Linewidth',2); hold on;
    plot(-yp,xp,'b','Linewidth',2); hold off;
    axis equal;
    ARS = sprintf('AR = %4.2f,   b = %5.3f,   S = %6.4f,  lambda = %3.2f', ...
                  AR, b, S, lambda);
    title(ARS,'Fontsize',12);
    
    % plot front view
    subplot(4,1,4,'align');
    plot( yp,zp,'b','Linewidth',4);hold on;
    plot(-yp,zp,'b','Linewidth',4);hold off;
    xlabel('$2y/b_{ref}$','Fontsize',12);
    dih = sprintf('dihedral = %4.2f (deg)', dihedral);
    title(dih,'Fontsize',12);
    axis equal;
    
end
