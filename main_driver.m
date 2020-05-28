clc; clear all; close all
addpath('/DFS-L/DATA/primeau/weilewang/my_func/')
addpath('/DFS-L/DATA/primeau/weilewang/DATA/')
addpath('/DFS-L/DATA/primeau/weilewang/GREG/Couple_CP/')
on = true; off = false;
global GC
% load constraint data
load transport_v4.mat
load GLODAP_grid_dic
load GLODAP_grid_Alk
load human_co2
load splco2_mod_monthly % monthly CO2 data
load co2syspar90.mat co2syspar

load Sobs_90x180x24.mat
load tempobs_90x180x24.mat
load po4obs_90x180x24.mat % WOA PO4 observation
load npp_90x180.mat % Satellite NPP in C unit;
load kw660.mat Kw660 p4
load tmpC.mat 
format long
grd  = grid         ;
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
% define some constants;
spd  = 24*60^2;
spa  = 365*spd;
zc   = sum(grd.zt(1:2));

format long
MSK_now = M3d;
parm.human_co2 = human_co2;
parm.Salt  = Sobs    ;
parm.Temp  = tempobs ;
parm.ALK   = TAstar  ;
parm.dVt   = dVt     ;
parm.Kw660 = Kw660   ;
parm.p4    = p4      ;
parm.c2p   = 110     ;
parm.M3d   = M3d     ;
parm.iwet  = iwet    ;
parm.nwet  = nwet    ;
parm.TRdiv = -TR     ;
parm.grd   = grd     ;
parm.I     = speye(nwet);
parm.aveT  = nanmean(tempobs(:,:,1:8),3);
parm.Tobs  = tempobs;
parm.rho   = 1024.5         ; % seawater density;
permil     = parm.rho*1e-3  ; % from umol/kg to mmol/m3;

parm.DICobs = DICstar*permil; % GLODAP dic obs [mmol/m3];
parm.TAobs  = TAstar*permil ; % GLODAP TA obs [mmol/m3];
parm.po4obs = po4obs        ;
GC = [real(DIC);zeros(3*nwet,1)];

% transiant CO2 concentraion;
parm.year      = splco2_mod(:,1) ;
parm.pco2_air  = splco2_mod(:,2) ;
parm.co2syspar = co2syspar       ;

iwet_msk = find(MSK_now(:))      ;

parm.kappa_g = 1/(1e6*spa); % geological restoring time [1/s];
parm.DIPbar  = nansum(po4obs(iwet).*dVt(iwet))/nansum(dVt(iwet)); % volume
parm.taup    = 720*60^2; % (s) pic dissolution time-scale
parm.tau_TA  = 1./parm.taup;

par.biogeochem.sigma    = 0.30;
par.biogeochem.slopep   = 0.00      ; % Martin curve exponent of POP
par.biogeochem.interpp  = 9.750e-01 ;
par.biogeochem.kappa_dp = 7.824e-08 ;
par.biogeochem.alpha    = 9.151e-03 ;
par.biogeochem.beta     = 4.807e-01 ;
par.biogeochem.d        = 4048      ;   % pic remin e-folding length scale (m)
                         % 

par.biogeochem.slopec   = 0         ; % Martin curve exponent of OC
par.biogeochem.interpc  = 1.015e+00 ;
par.biogeochem.kappa_dc = 9.569e-09 ;
par.biogeochem.kappa_da = 0.5e-7    ;
par.biogeochem.RR       = 5.294e-02 ;  % pic:poc ratio

par.biogeochem.opt_sigma = off;
par.biogeochem.opt_kappa_dp = on;
par.biogeochem.opt_slopep = off;
par.biogeochem.opt_interpp = on;
par.biogeochem.opt_alpha = on;
par.biogeochem.opt_beta = on;
par.biogeochem.opt_d = on;
par.biogeochem.opt_slopec = off;
par.biogeochem.opt_interpc = on;
par.biogeochem.opt_RR = on;
par.biogeochem.opt_kappa_dc = on;

p0 = [];
if (par.biogeochem.opt_sigma == on)
    sigma = par.biogeochem.sigma; lsigma = log(sigma);
    strt = length(p0) + 1;
    p0 = [p0; lsigma];
    par.pindx.lsigma = strt : length(p0);
end 

if (par.biogeochem.opt_slopep == on)
    slopep = par.biogeochem.slopep; lslopep = log(slopep);
    strt = length(p0) + 1;
    p0 = [p0; lslopep];
    par.pindx.lslopep = strt : length(p0);
end 

if (par.biogeochem.opt_interpp == on)
    interpp = par.biogeochem.interpp; linterpp = log(interpp);
    strt = length(p0) + 1;
    p0 = [p0; linterpp];
    par.pindx.linterpp = strt : length(p0);
end 

if (par.biogeochem.opt_alpha == on)
    alpha = par.biogeochem.alpha; lalpha = log(alpha);
    strt = length(p0) + 1;
    p0 = [p0; lalpha];
    par.pindx.lalpha = strt : length(p0);
end 

if (par.biogeochem.opt_beta == on)
    beta = par.biogeochem.beta; lbeta = log(beta);
    strt = length(p0) + 1;
    p0 = [p0; lbeta];
    par.pindx.lbeta = strt : length(p0);
end 

if (par.biogeochem.opt_kappa_dp == on)
    kappa_dp = par.biogeochem.kappa_dp; lkappa_dp = log(kappa_dp);
    strt = length(p0) + 1;
    p0 = [p0; lkappa_dp];
    par.pindx.lkappa_dp = strt : length(p0);
end 

if (par.biogeochem.opt_slopec == on)
    slopec = par.biogeochem.slopec; lslopec = log(slopec);
    strt = length(p0) + 1;
    p0 = [p0; lslopec];
    par.pindx.lslopec = strt : length(p0);
end 

if (par.biogeochem.opt_interpc == on)
    interpc = par.biogeochem.interpc; linterpc = log(interpc);
    strt = length(p0) + 1;
    p0 = [p0; linterpc];
    par.pindx.linterpc = strt : length(p0);
end 

if (par.biogeochem.opt_d == on)
    d = par.biogeochem.d; ld = log(d);
    strt = length(p0) + 1;
    p0 = [p0; ld];
    par.pindx.ld = strt : length(p0);
end 

if (par.biogeochem.opt_kappa_dc == on)
    kappa_dc = par.biogeochem.kappa_dc; lkappa_dc = log(kappa_dc);
    strt = length(p0) + 1;
    p0 = [p0; lkappa_dc];
    par.pindx.lkappa_dc = strt : length(p0);
end 

if (par.biogeochem.opt_RR == on)
    RR = par.biogeochem.RR; lRR = log(RR);
    strt = length(p0) + 1;
    p0 = [p0; lRR];
    par.pindx.lRR = strt : length(p0);
end 

parm.kappa_p = 1/(720*60^2) ;


parm.p2c = 0.006+0.0069*po4obs;
parm.nzo = 2;
%%%%%%% prepare NPP for the model %%%%%%%%
inan        = find(isnan(npp(:)));
npp(inan)   = 0;
parm.npp    = npp/(12*spd);
parm.Lambda = M3d*0;
parm.Lambda(:,:,1) = 0.5*(1/grd.dzt(1))*parm.p2c(:,:,1)./(1e-9+po4obs(:,:,1));
parm.Lambda(:,:,2) = 0.5*(1/grd.dzt(2))*parm.p2c(:,:,2)./(1e-9+po4obs(:,:,2));
parm.Lambda(:,:,3:end) = 0;
%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%
parm.p0 = p0;
x0    = p0;
myfun = @(x) neglogpost(x, parm, par);

options = optimoptions(@fminunc                  , ...
                       'Algorithm','trust-region', ...
                       'GradObj','on'            , ...
                       'Hessian','off'           , ...
                       'Display','iter'          , ...
                       'MaxFunEvals',2000        , ...
                       'MaxIter',2000            , ...
                       'TolX',1e-9               , ...
                       'TolFun',1e-9             , ...
                       'DerivativeCheck','off'   , ...
                       'FinDiffType','central'   , ...
                       'PrecondBandWidth',Inf)   ;
%
G_test = on        ;
nip    = length(x0) ;
if(G_test);
    dx = sqrt(-1)*eps.^3*eye(nip);
    for ii = 1:8 %nip
        x  = real(x0)+dx(:,ii)    ;
        parm.interpp  = exp(x(1)) ;
        % parm.slopep   = x(2)      ;
        parm.interpc  = exp(x(2)) ;
        % parm.slopec   = x(4)      ;
        parm.kappa_dp = exp(x(3)) ;
        parm.alpha    = exp(x(4)) ;
        parm.beta     = exp(x(5)) ;
        parm.d        = exp(x(6)) ;
        parm.kappa_c  = exp(x(7)) ;
        parm.RR       = exp(x(8)) ;
        [f,fx] = neglogpost(x, parm, par)       ;
        diff = real(fx(ii))-imag(f)/eps^3 ;
        fprintf('%i %e  \n',ii,diff)      ;
    end
    keyboard
else
    [xhat,fval,exitflag] = fminunc(myfun,x0,options);
    [f,fx] = neglogpost(xhat,parm,par);
    save(fname,'xhat')
end

fprintf('------------ END! ---------------');