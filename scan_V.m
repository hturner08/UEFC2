% YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

% Get various fixed parameters
UEFC = GetUEFC;

% Loop over AR and S and determine max speed

nAR = 41;
ARvals = linspace(  7,  14, nAR);

nS = 41;
Svals  = linspace(0.2, 0.5, nS);

N = zeros(nAR,nS);
V = zeros(nAR,nS);

Vmin = 1e6;
Vmax = 0;
CL = zeros(nAR,nS);
T  = zeros(nAR,nS);
db = zeros(nAR,nS);
CLmin = 1e6;
CLmax = 0;
Tmin = 1e6;
Tmax = 0;
dbmin = 1e6;
dbmax = 0;
for iAR = 1:nAR,
    for iS = 1:nS,
        ARtmp = ARvals(iAR);
        Stmp = Svals(iS);
        [Vtmp, Ntmp] = opt_V(ARtmp, Stmp);  
        V(iAR,iS) = Vtmp;
        N(iAR,iS) = Ntmp;
        if ((Vtmp > 0) & (Vtmp < Vmin)),
            Vmin = Vtmp;
        end
        if (Vtmp > Vmax),
           Vmax = Vtmp;
        end
        if (Vtmp > 0),
            CLtmp = GetCL(Ntmp, ARtmp, Stmp);
            CL(iAR,iS) = CLtmp;
            if (CLtmp < CLmin),
                CLmin = CLtmp;
            end
            if (CLtmp > CLmax),
                CLmax = CLtmp;
            end
            Ttmp = GetRequiredThrust(Ntmp, ARtmp, Stmp);
            T(iAR,iS) = Ttmp;
            if (Ttmp < Tmin),
                Tmin = Ttmp;
            end
            if (Ttmp > Tmax),
                Tmax = Ttmp;
            end            
            dbtmp = Getdb(Ntmp, ARtmp, Stmp);
            db(iAR,iS) = dbtmp;
            if (dbtmp < dbmin),
                dbmin = dbtmp;
            end
            if (dbtmp > dbmax),
                dbmax = dbtmp;
            end            
        end
    end
    fprintf('Completed %3.1f%% of AR,S scan\n',iAR/nAR*100);
end

figure(1);
[C_V,h_V] = contour(ARvals,Svals,V',linspace(Vmin,Vmax,21));
xlabel('AR'); ylabel('S (sq. m)');title('Flight speed (m/s)');
grid on;
clabel(C_V,h_V);

figure(2);
Nmax = max(max(N));
[C_N,h_N] = contour(ARvals,Svals,N',linspace(1,Nmax,11));
xlabel('AR'); ylabel('S (sq. m)');title('N');
grid on;
clabel(C_N,h_N);

figure(3);
[C_CL,h_CL] = contour(ARvals,Svals,CL',linspace(CLmin,CLmax,11));
xlabel('AR'); ylabel('S (sq. m)');title('CL');
grid on;
clabel(C_CL,h_CL);

figure(4);
[C_T,h_T] = contour(ARvals,Svals,T',linspace(Tmin,Tmax,11));
xlabel('AR'); ylabel('S (sq. m)');title('T (N)');
grid on;
clabel(C_T,h_T);

figure(5);
[C_db,h_db] = contour(ARvals,Svals,db',linspace(dbmin,dbmax,11));
xlabel('AR'); ylabel('S (sq. m)');title('d/b');
grid on;
clabel(C_db,h_db);

[Vopt,Iopt] = max(V(:));
[jAR,jS] = ind2sub(size(V),Iopt);
ARopt = ARvals(jAR);
Sopt = Svals(jS);

fprintf('\n\n');
fprintf('Maximum flight speed aircraft characteristics:\n');
fprintf('----------------------------------------------\n');
report_opt_V(ARopt,Sopt);

for ii = 1:5,
    figure(ii);
    hold on;
    plot(ARopt,Sopt,'r*','MarkerSize',14);
    hold off;
end
