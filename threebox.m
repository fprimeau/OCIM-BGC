% indices
n  = 3       ; % number of boxes
nx = 2       ; % number of surface boxes
k  = 1:n     ; % index of all ocean grid boxes
s_ndx = 1:nx ; % index of surface grid boxes
i_ndx = setdiff(k,s_ndx); % index of interior grid boxes

% volumes and areas from T99
Na  = 1.773e20     ; % molar volume of atmosphere
SA  = 3.49e14      ; % m^2 area of ocean
VT  = 1.292e18     ; % m^3 volume of ocean
dz  = [250; 100]   ; % m thickness of boxes P thru A
fa1 = 0.15         ;
fa  = [fa1; 1-fa1] ; % fractional surface areas
tmp = [0.0131e18; 0.0297e18] ;
VOL = [tmp; VT-sum(tmp)]     ;

% conversion factors
rho  = 1024.5   ; % kg m^-3
spyr = 3.1558e7 ; % seconds per year

% volume fluxes from T99 (Sv)
Th  = 20  ; % overturning circulation
fhd = 100 ;
fld = 0   ;
fhl = 0   ;

% diffusion (Sv)
A1 = [0, fhl, fhd; fhl, 0, fld; fhd, fld, 0];

% circulation (Sv)
A2 = [0, 0, Th; Th, 0, 0; 0, Th, 0];

% transport operator
A   = (A1 + A2) * 1e6 * spyr       ; % m^3/yr
TR1 = -diag(sum(A(k,:))) + A(:,k)' ;
TR  = inv(diag(VOL)) * TR1         ; % yr^-1

% temperature (K)
T = [2.0; 21.5] + 273.15 ;

% salinity (psu)
S = [1; 1] * 34.7 ;

% alkalinity (mol/kg)
ALK = zeros(n,1)+(587.05+50.560*34.7)./1e6 ;

% equilibrium constants
[k0,k1,k2,KB,BT,KW] = equic(T,S)  ;

% piston velocity (m/yr)
pv = zeros(nx, 1) + 5.5e-5 * spyr ;
  
% gas exchange coefficient for CO2 (mol/kg/atm/yr)
kappa = pv.*k0.*fa*SA./VOL(s_ndx) ;

% transport operator for carbon (4 x 4) (yr^-1)
AC = [zeros(1,4);[zeros(3,1),TR]] ;
  
%%% Solve for steady-state carbon distribution
% number of time-steps
l  = 1e5 ;
% record C distribution at each time step
C  = zeros(n+1,l) ;
% step forward in time with step dt (yr)
dt = .5 ;
% implicit time-stepping scheme C(t+1) = M\N, M holds stiff matrix
M  = eye(n+1)-(dt/2)*AC ;
% total carbon inventory
CT = 2.87805e18 ;
% a vector to ensure conservation of carbon
VC = [Na; VOL*rho]   ;
% initial condition (mix the carbon uniformally)
C(:,1) = CT./sum(VC) ;
for t  = 1 : l
    % compute surface ocean pco2 at each time step
    pco2 = cpco2(C(s_ndx+1,t),ALK(s_ndx),k1,k2,KB,BT,KW,k0);
    % compute gas exchange vector GC at each time step (atm/yr, mol/kg/yr)
    GC   = [[(-sum(kappa.*C(1,t).*VOL(s_ndx)) + ...
              sum(kappa.*pco2.*VOL(s_ndx)))*rho/Na];...
            [kappa.*(C(1,t)-pco2)]; ...
            [zeros(1,1)]] ;
    % step forward
    C(:,t+1) = M\(C(:,t) + dt*((AC/2)*C(:,t) + GC)) ;
    %C(1,t+1) = 0;
    % compute residual vector to decide if the model is in equilibrium
    R  = AC * C(:,t) + GC ;
    er = max(abs(R))      ;
    if(er < 1.5e-11)
        break
    end
end
C = C(:,1:t) ;

figure(1)
pco2atm = C(1,:)*1e6;
time = (1:t)*dt;
plot(time,pco2atm); 
title('Atmospheric pCO2 (ppm)')