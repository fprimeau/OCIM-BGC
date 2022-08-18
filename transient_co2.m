% transient anthropogenic co2 simulation
% path(path,'/Volumes/Ocean HD/RESEARCH/MOCM/BOMBC14/CO2SYS');
% Kw = a*(1-fice)*(u10^2)*(Sc/660)^-.5, u10 = sqrt(u2+v), u2 is
% the squared wind speed [m^2/s^2], v is the wind variance [m^2/s^2],
% and Sc is the Schmidt number for cfc [unitless]

secperyear = 365*24*60^2;

load co2data.mat co2data
co2syspar = co2data.co2syspar; % co2 system parameters
fice = co2data.fice;           % fraction of sea ice cover
u10 = co2data.u10;             % wind speed at 10 m above sea level
Sc660 = (co2data.Sc_co2./660); % Schmitd number

a = 20.5837;   % scales the wind-speed dependence of the air-sea
               % gas exchange operator
               
% load the atmospheric CO2 time history
atm = importdata('atmhist.txt');
pco2atm = atm.data(:,5);
year = atm.data(:,1);

% load the transport operator, wet-dry mask, and grid info
load  omega2_ad1e-05_ai1000.mat output
grd = output.grid;             % grid info
M3d = output.M3d;              % mask: wet=1 dry=0
TRdiv = -output.TR*secperyear; % transport operator [s^-1]
dic0 = output.dic; % initial guess for equilibrium
pic0 = 0*dic0;
%
iwet = find(M3d(:));           % indices to wet grid boxes
nh = length(output.msk.hkeep); % number of surface points



% piston velocity (m/yr) 
Kw = a*(1-fice).*(u10.^2).*Sc660.^-.5;


% assume the system is in equilibrium at time t = year(1)
t = year(1);
fprintf('solving for the pre-industrial equilibrium state...');
tic
% solve for equilibrium DIC with Newton's method
options.atol = 5e-7; options.rtol = 1e-16; options.iprint = 1;
load ../PO4/po4_sol.mat remin prod
rcp = 106;
remin = rcp*remin;
prod = rcp*prod;

p.RR = 0.081;  % pic:poc ratio
p.d = 1600; % pic remin e-folding length scale (m)
p.taup = 720*60^2; % (s) pic dissolution time-scale (Has not been
                   % calibrated yet!))

[dicpic,ierr,FDIC] = ...
    nsnew([dic0;pic0],@(x) eqdic(x,p,grd,M3d,TRdiv,Kw,pco2atm(1),co2syspar,remin,prod),options);
% need air-sea gas-exchange operator (Lgas)
[F,FD,Lgas] = eqdic(dicpic,p,grd,M3d,TRdiv,Kw,pco2atm(1),co2syspar,remin,prod);
toc

% now time-step from year(1) to length(year) using a semi-implicit scheme
% 
% advection-diffusion:   trapezoidal rule (a.k.a. Crank-Nicholson)
% air-sea gas exchange:  Euler forward 
% bio pump:              Euler forward                               
%      dic(n+1)-dic(n) +dt*TRdiv*((dic(n+1)+dic(n)/2) = dt*flux(dic(n)))
% ==> (I+(dt/2)*TRdiv)*dic(n+1) = (I-(dt/2)*TRdiv)*dic(n) + dt*flux(dic(n))
% ==>                A*dic(n+1) =                B*dic(n) + dt*flux(dic(n))
nstep_per_year = 5; % should be 5 or bigger.
dt = 1/nstep_per_year;
t = year(1):dt:year(end);

m = length(iwet);
dic = dicpic(1:m); 
pic = dicpic(m+1:end); 
% save the initial equilibrium state so that we can compute the
% anthropogenic inventory from the difference dic(t)-dic0
dic0 = dic;
dVt = grd.DXT3d.*grd.DYT3d.*grd.DZT3d;
dvol = dVt(iwet);
%
I = speye(length(dic));
A = I+(dt/2)*TRdiv;
B = I-(dt/2)*TRdiv;
%
fprintf('factoring the differential system propagator...');
tic
  FA = mfactor(A);
toc
% interpolate the annually sampled atmospheric pCO2 to frequency required
% by dt
pco2atm = interp1(year,pco2atm,t);
% initialize a transient DIC variable
DIC = zeros(length(dic0),length(year));
%
kk = 1;    
DIC(:,kk) = dic;
tt(kk) = t(1);
tic
fg1 = 1;

figure(fg1);
dV = M3d;
dV(iwet) = dvol;
tmp = M3d;
tmp(iwet) = dic0;
horiz_avg = @(v,dV) squeeze(sum(sum(v.*dV,1),2))./squeeze(sum(sum(dV,1),2));
havg0 = horiz_avg(tmp,dV);

for i = 2:length(t);
    [co2surf,R,k0] = eqco2(dic(1:nh),co2syspar);
    dic = mfactor(FA,(B*dic+dt*Lgas*(k0*pco2atm(1)-co2surf)+ ...
                      dt*(remin(iwet)-prod(iwet)) ...
                      -dt*p.RR*prod(iwet)+(dt/p.taup)*pic));
  
  if(mod(i,nstep_per_year)==1) % save solution once per year
    kk = kk+1;
    DIC(:,kk) = dic;
    inventory = 12*sum((dic-dic0).*dvol)/1e18; % PgC
    tmp(iwet) = dic;
    havg = horiz_avg(tmp,dV);
    set(0,'CurrentFigure',fg1)
    clf
    plot(havg0,-grd.zt,'-sr'); hold on
    plot(havg,-grd.zt,'-sb');  
    ylabel('depth (m)');
    xlabel('dic (mmol/m3)');
    s = sprintf('      [C_{ant} inventory = %3.1f PgC]',inventory);
    title([num2str(t(i)),s]);
    grid on
    drawnow

    % compute the inventory of 
    tt(kk) = t(i);
    %
    fprintf('t = %4.1f pco2atm = %3.1f max(pco2ocn) = %3.1f min(pco2ocn) = %3.1f inventory = %3.1f PgC\n ',...
            tt(kk),pco2atm(i),max(co2surf./k0),min(co2surf./k0),inventory)
  end
end

