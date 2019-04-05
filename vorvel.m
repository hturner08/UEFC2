function vel = vorvel(ra,rb,r,ibound)

% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

N  = size(r,2);
a  = r - repmat(ra,1,N);
b  = r - repmat(rb,1,N);
ih = repmat([1;0;0],1,N);

am  = sqrt(dot(a,a));
bm  = sqrt(dot(b,b));
adb = dot(a,b);
axb = cross(a,b);
axi = cross(a,ih);
bxi = cross(b,ih);

den  = am.*bm + adb;
dena = am - a(1,:);
denb = bm - b(1,:);

v = repmat((1./am)./dena,3,1).*axi ...
  - repmat((1./bm)./denb,3,1).*bxi;

if(ibound == 1) 
  v = v + repmat((1./am + 1./bm)./den,3,1).*axb;
end

vel = v/(4*pi);
