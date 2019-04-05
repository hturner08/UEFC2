function[yv,zv,cl,ccl,vi,wi,alpha,beta,pbar,rbar,CL,CDi,Cr,Cn,Cb] = wvl(geom,N,ispace,Sref,bref,cref,itmax,toler, ...
 alspec, bespec, pbspec, rbspec, CLspec, Crspec, Cnspec, ...
ialspec,ibespec,ipbspec,irbspec,iCLspec,iCrspec,iCnspec)
%========================================================================
%  Weissinger Vortex Lattice program for 3D wing aerodynamic analysis,
%     with dihedral, sideslip, roll-rate, and yaw-rate effects.
%========================================================================

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Number of geometry sections
ns = size(geom,1); 

% Calculate running yz arc length along span
ssec = zeros(ns,1);
for i = 1:ns-1
  dy = geom(i+1,2)-geom(i,2);
  dz = geom(i+1,3)-geom(i,3);
  ssec(i+1) = ssec(i) + sqrt(dy*dy + dz*dz);
end

% Number of flow conditions
nflow = size(alspec,2);

% Number of flow parameters  (alpha,beta,pbar,rbar)
npar = 4;

if (ispace == 1)
% uniform spacing
  dt = ssec(ns)/N;
  t = 0:dt:ssec(ns);
  s = 0:dt:ssec(ns);
else
% cosine spacing
  dt = 0.5*pi/N;
  t = 0:dt:0.5*pi;
  s = ssec(ns)*cos(0.5*pi-t);
end

% set quantities at nodes by interpolation
c = interp1(ssec,geom(:,4),s);
x = interp1(ssec,geom(:,1),s) + 0.25*c;
y = interp1(ssec,geom(:,2),s);
z = interp1(ssec,geom(:,3),s);

% set quantities at vortex midpoints
if (ispace == 1)
 tv = t(2:end) - 0.5*dt;
 sv = s(2:end) - 0.5*dt;
else
 tv = t(2:end) - 0.5*dt;
 sv = ssec(ns)*cos(0.5*pi-tv);
end

cv = interp1(ssec,geom(:,4),sv);
xv = interp1(ssec,geom(:,1),sv) + 0.25*cv;
yv = interp1(ssec,geom(:,2),sv);
zv = interp1(ssec,geom(:,3),sv);

% control point locations
xc = xv + 0.5*cv;
yc = yv;
zc = zv;

% sin, cos of local dihedral angle
sind = zeros(1,N);
cosd = zeros(1,N);
for i=1:N
 sind(1,i) = (z(i+1)-z(i))/(s(i+1)-s(i));
 cosd(1,i) = (y(i+1)-y(i))/(s(i+1)-s(i));
end

twv = interp1(ssec,geom(:,5),sv);
a0v = interp1(ssec,geom(:,6),sv);
twa = twv-a0v;

% Reflect geometry for left wing
c  = [ fliplr(c), c(2:end)];
x  = [ fliplr(x), x(2:end)];
y  = [-fliplr(y), y(2:end)];
z  = [ fliplr(z), z(2:end)];
s  = [-fliplr(s), s(2:end)];

cv  = [ fliplr(cv), cv ];
xv  = [ fliplr(xv), xv ];
yv  = [-fliplr(yv), yv ];
zv  = [ fliplr(zv), zv ];

xc  = [ fliplr(xc), xc ];
yc  = [-fliplr(yc), yc ];
zc  = [ fliplr(zc), zc ];

twa = [ fliplr(twa) , twa ];
cost = cos(twa);
sint = sin(twa);

cosd = [ fliplr(cosd), cosd ];
sind = [-fliplr(sind), sind ];


%---------------------------------------------------------------------
% check specified-quantity vector size consistency

nCL = size(CLspec,2);
if (nCL == 1)
  CLspec = CLspec * ones(nflow,1);
elseif (nCL ~= nflow)
  error('CLspec must be a scalar, or a vector of same length as alspec');
end

nbe = size(bespec,2);
if (nbe == 1)
  bespec = bespec * ones(nflow,1);
elseif (nbe ~= nflow)
  error('bespec must be a scalar, or a vector of same length as alspec');
end

npb = size(pbspec,2);
if (npb == 1)
  pbspec = pbspec * ones(nflow,1);
elseif (npb ~= nflow)
  error('pbspec must be a scalar, or a vector of same length as alspec');
end

nrb = size(rbspec,2);
if (nrb == 1)
  rbspec = rbspec * ones(nflow,1);
elseif (nrb ~= nflow)
  error('rbspec must be a scalar, or a vector of same length as alspec');
end

nCr = size(Crspec,2);
if (nCr == 1)
  Crspec = Crspec * ones(nflow,1);
elseif (nCr ~= nflow)
  error('Crspec must be a scalar, or a vector of same length as alspec');
end

nCn = size(Cnspec,2);
if (nCn == 1)
  Cnspec = Cnspec * ones(nflow,1);
elseif (nCn ~= nflow)
  error('Cnspec must be a scalar, or a vector of same length as alspec');
end

%---------------------------------------------------------------------


% set up AIC matrix
A = zeros(2*N,2*N);
for i=1:2*N
 vvor = vorvel([x(i);y(i);z(i)],[x(i+1);y(i+1);z(i+1)],[xc;yc;zc],1);
 A(i,:) = vvor(1,:).*sint + vvor(3,:).*cost;
end


% initialize parameters to be computed by iteration loop below
alpha= zeros(nflow,1);
beta = zeros(nflow,1);
pbar = zeros(nflow,1);
rbar = zeros(nflow,1);

% set nonzero changes to force at least one iteration loop pass
dalpha= ones(nflow,1);
dbeta = ones(nflow,1);
dpbar = ones(nflow,1);
drbar = ones(nflow,1);


iter = 0;

istop = 0;  % assume convergence will be achieved

%==============================================================
% Nweton iteration loop to converge to specified parameters

while max( [ max(abs(dalpha)) max(abs(dbeta)) max(abs(dpbar)) max(abs(drbar)) ] ) > toler

    %fprintf('\n Iteration %5i ...\n',iter);
    %fprintf(' alpha=%7.3f  %7.3f\n', alpha*180/pi);
    %fprintf(' beta =%7.3f  %7.3f\n', beta*180/pi);
%fprintf(' pbar =%7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f\n', pbar);
%fprintf(' rbar =%7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f\n', rbar);

if (iter >= itmax) 
  fprintf('\n Iteration limit exceeded\n');
  istop = 1;
  break;
end

if (max( abs(alpha) )*180/pi > 45.0)
  fprintf('Alpha > 45 deg reached. CLspec is probably excessive.\n')
  istop = 2;
  break;
end

if (max( abs(beta) )*180/pi > 45.0)
  fprintf('Beta > 45 deg reached, likely due to insufficient dihedral.\n')
  istop = 2;
  break;
end


sina = sin(alpha);
cosa = cos(alpha);
sinb = sin(beta);
cosb = cos(beta);

% set up rhs vectors for each flow
R = zeros(nflow,2*N);
R_a = zeros(nflow,2*N);
R_b = zeros(nflow,2*N);
R_p = zeros(nflow,2*N);
R_r = zeros(nflow,2*N);
for iflow=1:nflow
  Vx =  cosa(iflow)*cosb(iflow) - 2*yc*rbar(iflow);
  Vy = -sinb(iflow)             - 2*zc*pbar(iflow);
  Vz =  sina(iflow)*cosb(iflow) + 2*yc*pbar(iflow);

  Vx_a = -sina(iflow)*cosb(iflow);
  Vy_a = 0;
  Vz_a =  cosa(iflow)*cosb(iflow);

  Vx_b = -cosa(iflow)*sinb(iflow);
  Vy_b = -cosb(iflow);
  Vz_b = -sina(iflow)*sinb(iflow);

  Vx_p =  0;
  Vy_p = -zc;
  Vz_p =  yc;

  Vx_r = -yc;
  Vy_r = 0;
  Vz_r = 0;

  R(iflow,:)   = -(Vx.*sint + Vz.*cost.*cosd - Vy.*sind);
  R_a(iflow,:) = -(Vx_a.*sint + Vz_a.*cost.*cosd - Vy_a.*sind);
  R_b(iflow,:) = -(Vx_b.*sint + Vz_b.*cost.*cosd - Vy_b.*sind);
  R_p(iflow,:) = -(Vx_p.*sint + Vz_p.*cost.*cosd - Vy_p.*sind);
  R_r(iflow,:) = -(Vx_r.*sint + Vz_r.*cost.*cosd - Vy_r.*sind);
end

% circulation vector G = Gamma/Vinf for each flow, and sensitivities wrt parameters
G = R/A;
G_a = R_a/A;
G_b = R_b/A;
G_p = R_p/A;
G_r = R_r/A;

% wake induced velocities
vi = zeros(nflow,2*N);
wi = zeros(nflow,2*N);

vi_a = zeros(nflow,2*N);
vi_b = zeros(nflow,2*N);
vi_p = zeros(nflow,2*N);
vi_r = zeros(nflow,2*N);

wi_a = zeros(nflow,2*N);
wi_b = zeros(nflow,2*N);
wi_p = zeros(nflow,2*N);
wi_r = zeros(nflow,2*N);

for i=1:2*N
  rsqi = (yv-y(i)).^2 + (zv-z(i)).^2;
  rsqp = (yv-y(i+1)).^2 + (zv-z(i+1)).^2;
  wyG =  ((zv-z(i))./rsqi - (zv-z(i+1))./rsqp)/(4*pi);
  wzG = -((yv-y(i))./rsqi - (yv-y(i+1))./rsqp)/(4*pi);
  vi(1:nflow,:) = vi(1:nflow,:) + G(1:nflow,i)*wyG;
  wi(1:nflow,:) = wi(1:nflow,:) + G(1:nflow,i)*wzG;

  vi_a(1:nflow,:) = vi_a(1:nflow,:) + G_a(1:nflow,i)*wyG;
  vi_b(1:nflow,:) = vi_b(1:nflow,:) + G_b(1:nflow,i)*wyG;
  vi_p(1:nflow,:) = vi_p(1:nflow,:) + G_p(1:nflow,i)*wyG;
  vi_r(1:nflow,:) = vi_r(1:nflow,:) + G_r(1:nflow,i)*wyG;

  wi_a(1:nflow,:) = wi_a(1:nflow,:) + G_a(1:nflow,i)*wzG;
  wi_b(1:nflow,:) = wi_b(1:nflow,:) + G_b(1:nflow,i)*wzG;
  wi_p(1:nflow,:) = wi_p(1:nflow,:) + G_p(1:nflow,i)*wzG;
  wi_r(1:nflow,:) = wi_r(1:nflow,:) + G_r(1:nflow,i)*wzG;
end

cl  = zeros(nflow,2*N);  % local L'/(0.5 rho V^2 c)        =  cl
ccl = zeros(nflow,2*N);  % local L'/(0.5 rho Vinf^2 cref)  =  cl*(c/cref)*(V/Vinf)^2   (local loading)
CL = zeros(nflow,1);  % overall CL
CDi= zeros(nflow,1);  % overall CDi
Cr = zeros(nflow,1);  % overall Cr  (roll moment coefficient)
Cn = zeros(nflow,1);  % overall Cn  (yaw  moment coefficient)
Cb = zeros(nflow,1);  % root bending moment coefficient

CL_a = zeros(nflow,1);
CL_b = zeros(nflow,1);
CL_p = zeros(nflow,1);
CL_r = zeros(nflow,1);

Cr_a = zeros(nflow,1);
Cr_b = zeros(nflow,1);
Cr_p = zeros(nflow,1);
Cr_r = zeros(nflow,1);

Cn_a = zeros(nflow,1);
Cn_b = zeros(nflow,1);
Cn_p = zeros(nflow,1);
Cn_r = zeros(nflow,1);

for i=1:2*N
  dx = x(i+1)-x(i);
  dy = y(i+1)-y(i);
  dz = z(i+1)-z(i);

% normalized local velocity relative to wing station,  (Vx,Vy,Vz) = (u,v,w)/Vinf
  Vx =  cosa.*cosb - yc(i)*rbar;
  Vy = -sinb       - zc(i)*pbar;
  Vz =  sina.*cosb + yc(i)*pbar;

  Vx_a = -sina.*cosb;
  Vy_a = 0;
  Vz_a =  cosa.*cosb;

  Vx_b = -cosa.*sinb;
  Vy_b = -cosb;
  Vz_b = -sina.*sinb;

  Vx_p =  0;
  Vy_p = -zc(i);
  Vz_p =  yc(i);

  Vx_r = -yc(i);
  Vy_r = 0;
  Vz_r = 0;
  
  Vperp = sqrt(Vx.^2 + Vz.^2);
  Vperp_Vx = Vx/Vperp;
  Vperp_Vz = Vz/Vperp;
  
  Vperp_a = Vperp_Vx*Vx_a + Vperp_Vz*Vz_a;
  Vperp_b = Vperp_Vx*Vx_b + Vperp_Vz*Vz_b;
  Vperp_p = Vperp_Vx*Vx_p + Vperp_Vz*Vz_p;
  Vperp_r = Vperp_Vx*Vx_r + Vperp_Vz*Vz_r;
%
  cl(1:nflow,i)  = 2*G(1:nflow,i)./(Vperp*cv(i));
  ccl(1:nflow,i) = 2*G(1:nflow,i).*Vperp/cref;

  CL(1:nflow) = CL(1:nflow) + 2*G(1:nflow,i).* Vperp*dy/Sref;
  Cr(1:nflow) = Cr(1:nflow) - 2*G(1:nflow,i).* Vx*(yv(i)*dy+zv(i)*dz)/(bref*Sref);
  Cn(1:nflow) = Cn(1:nflow) - 2*G(1:nflow,i).*(Vz+wi(1:nflow,i))* yv(i)*dy/(bref*Sref);

  Cb(1:nflow) = Cb(1:nflow) + 2*G(1:nflow,i).* Vx*abs(yv(i)*dy+zv(i)*dz)/(bref*Sref);
  CDi(1:nflow)= CDi(1:nflow)- 2*G(1:nflow,i).*(wi(1:nflow,i)*dy-vi(1:nflow,i)*dz)/Sref;

  CL_a(1:nflow) = CL_a(1:nflow) + 2*(G_a(1:nflow,i).*Vperp ...
                                     + G(1:nflow,i).*Vperp_a ) * dy/Sref;
  Cr_a(1:nflow) = Cr_a(1:nflow) - 2*(G_a(1:nflow,i).*Vx ...
                                     + G(1:nflow,i).*Vx_a ) * (yv(i)*dy+zv(i)*dz)/(bref*Sref);
  Cn_a(1:nflow) = Cn_a(1:nflow) - 2*(G_a(1:nflow,i).*(Vz+wi(1:nflow,i)) ...
                                     + G(1:nflow,i).*(Vz_a+wi_a(1:nflow,i) )) * yv(i)*dy/(bref*Sref);

  CL_b(1:nflow) = CL_b(1:nflow) + 2*(G_b(1:nflow,i).*Vperp ...
                                     + G(1:nflow,i).*Vperp_b ) * dy/Sref;
  Cr_b(1:nflow) = Cr_b(1:nflow) - 2*(G_b(1:nflow,i).*Vx ...
                                     + G(1:nflow,i).*Vx_b ) * (yv(i)*dy+zv(i)*dz)/(bref*Sref);
  Cn_b(1:nflow) = Cn_b(1:nflow) - 2*(G_b(1:nflow,i).*(Vz+wi(1:nflow,i)) ...
                                     + G(1:nflow,i).*(Vz_b+wi_b(1:nflow,i) )) * yv(i)*dy/(bref*Sref);

  CL_p(1:nflow) = CL_p(1:nflow) + 2*(G_p(1:nflow,i).*Vperp ...
                                     + G(1:nflow,i).*Vperp_p ) * dy/Sref;
  Cr_p(1:nflow) = Cr_p(1:nflow) - 2*(G_p(1:nflow,i).*Vx ...
                                     + G(1:nflow,i).*Vx_p ) * (yv(i)*dy+zv(i)*dz)/(bref*Sref);
  Cn_p(1:nflow) = Cn_p(1:nflow) - 2*(G_p(1:nflow,i).*(Vz+wi(1:nflow,i)) ...
                                     + G(1:nflow,i).*(Vz_p+wi_p(1:nflow,i) )) * yv(i)*dy/(bref*Sref);

  CL_r(1:nflow) = CL_r(1:nflow) + 2*(G_r(1:nflow,i).*Vperp ...
                                     + G(1:nflow,i).*Vperp_r ) * dy/Sref;
  Cr_r(1:nflow) = Cr_r(1:nflow) - 2*(G_r(1:nflow,i).*Vx ...
                                     + G(1:nflow,i).*Vx_r ) * (yv(i)*dy+zv(i)*dz)/(bref*Sref);
  Cn_r(1:nflow) = Cn_r(1:nflow) - 2*(G_r(1:nflow,i).*(Vz+wi(1:nflow,i)) ...
                                     + G(1:nflow,i).*(Vz_r+wi_r(1:nflow,i) )) * yv(i)*dy/(bref*Sref);
end

% Induced drag scales inversely with square of projected span which is reduced by cos(sideslip)
CDi = CDi ./ cosb.^2;

%-------------------------------------------------------------------
% Update parameters for each flow condition
for iflow = 1:nflow

% Set up Newton system for parameter changes
  nsys = npar;
  Asys = zeros(nsys,nsys);
  Rsys = zeros(nsys,1);

  ksys = 0;
  if (ialspec(iflow) == 1)
    ksys = ksys+1;
    Rsys(ksys) = alpha(iflow) - alspec(iflow);
    Asys(ksys,1) = 1;
  end
  if (ibespec(iflow) == 1)
    ksys = ksys+1;
    Rsys(ksys) = beta(iflow) - bespec(iflow);
    Asys(ksys,2) = 1;
  end
  if (ipbspec(iflow) == 1)
    ksys = ksys+1;
    Rsys(ksys) = pbar(iflow) - pbspec(iflow);
    Asys(ksys,3) = 1;
  end
  if (irbspec(iflow) == 1)
    ksys = ksys+1;
    Rsys(ksys) = rbar(iflow) - rbspec(iflow);
    Asys(ksys,4) = 1;
  end
  if (iCLspec(iflow) == 1)
    ksys = ksys+1;
    Rsys(ksys) = CL(iflow) - CLspec(iflow);
    Asys(ksys,1) = CL_a(iflow);
    Asys(ksys,2) = CL_b(iflow);
    Asys(ksys,3) = CL_p(iflow);
    Asys(ksys,4) = CL_r(iflow);
  end
  if (iCrspec(iflow) == 1)
    ksys = ksys+1;
    Rsys(ksys) = Cr(iflow) - Crspec(iflow);
    Asys(ksys,1) = Cr_a(iflow);
    Asys(ksys,2) = Cr_b(iflow);
    Asys(ksys,3) = Cr_p(iflow);
    Asys(ksys,4) = Cr_r(iflow);
  end
  if (iCnspec(iflow) == 1)
    ksys = ksys+1;
    Rsys(ksys) = Cn(iflow) - Cnspec(iflow);
    Asys(ksys,1) = Cn_a(iflow);
    Asys(ksys,2) = Cn_b(iflow);
    Asys(ksys,3) = Cn_p(iflow);
    Asys(ksys,4) = Cn_r(iflow);
  end

  if (ksys ~= nsys)
    fprintf('Error: %3i quantities are specified for flow condition %3i .  Must have 4. \n',ksys, iflow);
    error('Stopping execution');
  end

% Asys
  dpar = -Asys\Rsys;

  dalpha(iflow)= dpar(1);
  dbeta(iflow) = dpar(2);
  dpbar(iflow) = dpar(3);
  drbar(iflow) = dpar(4);
%-------------------------------------------------------------------

end

% Update parameters
alpha= alpha+ dalpha;
beta = beta + dbeta;
pbar = pbar + dpbar;
rbar = rbar + drbar;

iter = iter + 1;

end  % of parameter iteration loop
%================================================

if (istop ~= 0)
  error('SOLUTION NOT CONVERGED')
end

