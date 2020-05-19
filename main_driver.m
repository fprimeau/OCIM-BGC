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
load GC.mat GC
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
GC = real(GC);

% transiant CO2 concentraion;
parm.year      = splco2_mod(:,1) ;
parm.pco2_air  = splco2_mod(:,2) ;
parm.co2syspar = co2syspar       ;

iwet_msk = find(MSK_now(:))      ;

parm.kappa_g = 1/(1e6*spa); % geological restoring time [1/s];
parm.DIPbar  = nansum(po4obs(iwet).*dVt(iwet))/nansum(dVt(iwet)); % volume
parm.taup    = 720*60^2; % (s) pic dissolution time-scale
parm.tau_TA  = 1./parm.taup;
sigma    = 0.30      ;
slopep   = 0.00      ; % Martin curve exponent of POP
interpp  = 9.750e-01 ;
kappa_dp =  7.824e-08;
alpha    =  9.151e-03;
beta     = 4.807e-01 ;
d        = 4048      ;   % pic remin e-folding length scale (m)
% 

slopec   = 0         ; % Martin curve exponent of OC
interpc  = 1.015e+00 ;
kappa_dc = 9.569e-09 ;
kappa_da = 0.5e-7    ;
RR       = 5.294e-02 ;  % pic:poc ratio

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

x       = [interpp;slopep;interpc;slopec;sigma;kappa_dp;alpha;beta;d;kappa_dc;RR];
% x0      = [x0(1:4); log(x0(5:end))];
x0      = [log(x([1,3,6:end]))];
parm.x  = x;   % all interested parameters;
parm.x0 = x0;  % all parameters to be trained;
myfun   = @(x) neglogpost(x,parm);

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
G_test = off        ;
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
        parm.RR       = exp(x(8));
        [f,fx] = neglogpost(x,parm)       ;
        diff = real(fx(ii))-imag(f)/eps^3 ;
        fprintf('%i %e  \n',ii,diff)      ;
    end
    keyboard
else
    [xhat,fval,exitflag] = fminunc(myfun,x0,options);
    [f,fx] = neglogpost(xhat,parm);
    save(fname,'xhat')
end

fprintf('------------ END! ---------------');