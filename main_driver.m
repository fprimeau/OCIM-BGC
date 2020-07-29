clc; clear all; close all
addpath('/DFS-L/DATA/primeau/weilewang/my_func/')
addpath('/DFS-L/DATA/primeau/weilewang/DATA/')
addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
addpath('/DFS-L/DATA/primeau/weilewang/GREG/Couple_CP/')
format long
on = true; off = false;
global GC GO iter
iter = 0;
spd  = 24*60^2;
spa  = 365*spd;
%
version = 90;
par.optim  = on;
par.Cmodel = on;
par.Omodel = on;
par.Simodel = off;
%
if version == 90
    % ATTENTION: please change this directory to where you wanna
    % save your output files
    parm.VER = '/DFS-L/DATA/primeau/weilewang/COP4WWF/MSK90/fixedPO';
    
    load transport_v4.mat grid M3d TR
    load Sobs_90x180x24.mat
    load tempobs_90x180x24.mat
    load po4obs_90x180x24.mat % WOA PO4 observation
    load Siobs_91x180x24.mat Siobs
    %
    load DICant_90x180x24.mat
    load GLODAPv2_90x180x24raw.mat
    load splco2_mod_monthly % monthly CO2 data
    load co2syspar90.mat co2syspar
    load cbpm_npp_annual_90x180.mat
    load kw660_90x180.mat
    load /DFS-L/DATA/primeau/weilewang/OutputCoupledCPO/MSK90/temp_dep_b_C.mat
    load /DFS-L/DATA/primeau/weilewang/OutputCoupledCPO/MSK90/temp_dep_b_O2.mat
    grd  = grid ;
    
elseif version == 91
    % ATTENTION: please change this directory to where you wanna
    % save your output files
    parm.VER = '/DFS-L/DATA/primeau/weilewang/COP4WWF/MSK91/CTL_He_fixedPO';

    load OCIM2_CTL_He.mat output 
    % load OCIM2_KiLOW_He.mat output 
    % load OCIM2_KiHIGH_He.mat output 
    % load OCIM2_KvHIGH_He.mat output 
    % load OCIM2_KvHIGH_KiLOW_He.mat output 
    
    % load OCIM2_CTL_noHe.mat output 
    % load OCIM2_KiLOW_noHe.mat output 
    % load OCIM2_KiHIGH_noHe.mat output 
    % load OCIM2_KvHIGH_noHe.mat output 
    % load OCIM2_KvHIGH_KiLOW_noHe.mat output 
    % load OCIM2_KvHIGH_KiHIGH_noHe.mat output 
    load Sobs_91x180x24.mat
    load po4obs_91x180x24.mat % WOA PO4 observation
    load tempobs_91x180x24.mat
    load Siobs_91x180x24.mat Siobs
    %
    load DICant_91x180x24.mat
    load GLODAPv2_91x180x24raw.mat
    load splco2_mod_monthly % monthly CO2 data
    load co2syspar91.mat co2syspar
    load cbpm_npp_annual_91x180.mat
    load kw660_91x180.mat
    load /DFS-L/DATA/primeau/weilewang/COP4WWF/MSK91/CTL_He_fixedPO_C.mat
    load /DFS-L/DATA/primeau/weilewang/COP4WWF/MSK91/CTL_He_fixedPO_O2.mat
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
end

iwet = find(M3d(:)) ;
nwet = length(iwet) ;
dVt = grd.DXT3d.*grd.DYT3d.*grd.DZT3d;
%
parm.Salt  = Sobs    ;
parm.Temp  = tempobs ;
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
                              %
parm.SIL     = Siobs;
parm.po4obs  = po4obs;
parm.o2raw   = o2raw;
parm.po4raw  = po4raw;
parm.sio4raw = sio4raw;
parm.dicraw  = dicraw*permil; % GLODAP dic obs [mmol/m3];
parm.human_co2 = DICant*permil;

GC = real([DIC(iwet); POC(iwet); DOC(iwet); CaC(iwet)]);
GO = real(O2(iwet)) + 1e-5*randn(parm.nwet,1);

% transiant CO2 concentraion;
parm.year      = splco2_mod(:,1) ;
parm.pco2_air  = splco2_mod(:,2) ;
parm.co2syspar = co2syspar       ;

parm.kappa_g = 1/(1e6*spa); % geological restoring time [1/s];
parm.SILbar  = nansum(Siobs(iwet).*dVt(iwet))/nansum(dVt(iwet));
parm.DIPbar  = nansum(po4obs(iwet).*dVt(iwet))/nansum(dVt(iwet)); % volume
parm.taup    = 720*60^2; % (s) pic dissolution time-scale
parm.tau_TA  = 1./parm.taup;
par.kappa_da = 0.5e-7    ;
par.kappa_gs = 1/(1e6*spa); % geological restoring time [1/s];
% PME part;
[sst,ss] = PME(parm) ;
parm.ss  = ss        ;
parm.sst = sst       ;
% load optimal parameters if they exist
fname = strcat(parm.VER,'_xhat.mat');
if isfile(fname)
    load(fname)
end

% P model parameters;
if exist('xhat') & isfield(xhat,'sigma')
    par.sigma = xhat.sigma;
else 
    par.sigma = 1/3  ;
end
if exist('xhat') & isfield(xhat,'kappa_dp')
    par.kappa_dp = xhat.kappa_dp;
else 
    par.kappa_dp = 4.44e-08 ;
end 
if exist('xhat') & isfield(xhat,'slopep')
    par.slopep = xhat.slopep;
else 
    par.slopep = 0 ;
end 
if exist('xhat') & isfield(xhat,'interpp')
    par.interpp  = xhat.interpp ;
else 
    par.interpp  = 0.89 ;
end 
if exist('xhat') & isfield(xhat,'alpha')
    par.alpha = xhat.alpha;
else 
    par.alpha = 9.33e-04 ;
end 
if exist('xhat') & isfield(xhat,'beta')
    par.beta = xhat.beta;
else
    par.beta = 1.16e-01 ;
end 

% C model parameters                                      
if exist('xhat') & isfield(xhat,'slopec')
    par.slopec = xhat.slopec ;
else
    par.slopec = 0 ;
end 
if exist('xhat') & isfield(xhat,'interpc')
    par.interpc = xhat.interpc ;
else
    par.interpc = 1.06e+00 ;
end 
if exist('xhat') & isfield(xhat,'d')
    par.d = xhat.d   ;
else
    par.d = 2.25e+03 ;
end 
if exist('xhat') & isfield(xhat,'kappa_dc')
    par.kappa_dc = xhat.kappa_dc ;
else 
    par.kappa_dc = 3.06e-08 ;
end 
if exist('xhat') & isfield(xhat,'RR')
    par.RR = xhat.RR  ;
else
    par.RR = 6.37e-02 ;
end
if exist('xhat') & isfield(xhat,'cc')
    par.cc = xhat.cc  ;
else
    par.cc = 5.77e-03 ;
end 
if exist('xhat') & isfield(xhat,'dd')
    par.dd = xhat.dd  ;
else 
    par.dd = 3.39e-03 ;
end 
%
% O model parameters
if exist('xhat') & isfield(xhat,'slopeo')
    par.slopeo = xhat.slopeo ;
else 
    par.slopeo = 0.0e+00 ;
end 
if exist('xhat') & isfield(xhat,'interpo')
    par.interpo = xhat.interpo ;
else 
    par.interpo = 1.70e+02 ;
end 
%
% Si model parameters
if exist('xhat') & isfield(xhat,'bsi')
    par.bsi = xhat.bsi;
else
    par.bsi = 0.33;
end 
if exist('xhat') & isfield(xhat,'at')
    par.at = xhat.at;
else
    par.at = 1.32e16/spd;
end 
if exist('xhat') & isfield(xhat,'bt')
    par.bt = xhat.bt ;
else 
    par.bt = 11481;
end 
if exist('xhat') & isfield(xhat,'aa')
    par.aa = xhat.aa ;
else
    par.aa = 1;
end 
if exist('xhat') & isfield(xhat,'bb')
    par.bb = xhat.bb;
else 
    par.bb = 0.968;
end 

%% -------------------------------------------------------------
%
% P model parameters
par.opt_beta = on;
par.opt_alpha = on;
par.opt_sigma = off; 
par.opt_slopep = on; 
par.opt_interpp = on;
par.opt_kappa_dp = on;
% C model parameters
par.opt_d  = on; % 
par.opt_RR = on; % 
par.opt_cc = on;
par.opt_dd = on;
par.opt_slopec   = on;
par.opt_interpc  = on; %
par.opt_kappa_dc = on; %
% O model parameters
par.opt_slopeo  = off; 
par.opt_interpo = off; 
% Si model parameters
par.opt_bsi = on;
par.opt_at = on;
par.opt_bt = off;
par.opt_aa = on;
par.opt_bb = on;
%% -------------------------------------------------------------
p0 = [];
% sigma 
if (par.opt_sigma == on)
    sigma = par.sigma; lsigma = log(sigma);
    strt = length(p0) + 1;
    p0 = [p0; lsigma];
    par.pindx.lsigma = strt : length(p0);
end 
% kappa_dp 
if (par.opt_kappa_dp == on)
    kappa_dp = par.kappa_dp; lkappa_dp = log(kappa_dp);
    strt = length(p0) + 1;
    p0 = [p0; lkappa_dp];
    par.pindx.lkappa_dp = strt : length(p0);
end
% slopep
if (par.opt_slopep == on)
    slopep = par.slopep;
    strt = length(p0) + 1;
    p0 = [p0; slopep];
    par.pindx.slopep = strt : length(p0);
end 
% interpp
if (par.opt_interpp == on)
    interpp = par.interpp; linterpp = log(interpp);
    strt = length(p0) + 1;
    p0 = [p0; linterpp];
    par.pindx.linterpp = strt : length(p0);
end 
% alpha 
if (par.opt_alpha == on)
    alpha = par.alpha; lalpha = log(alpha);
    strt = length(p0) + 1;
    p0 = [p0; lalpha];
    par.pindx.lalpha = strt : length(p0);
end 
% beta
if (par.opt_beta == on)
    beta = par.beta; lbeta = log(beta);
    strt = length(p0) + 1;
    p0 = [p0; lbeta];
    par.pindx.lbeta = strt : length(p0);
end
%
if (par.Cmodel == on)
    % slopec
    if (par.opt_slopec == on)
        slopec = par.slopec;
        strt = length(p0) + 1;
        p0 = [p0; slopec];
        par.pindx.slopec = strt : length(p0);
    end 
    % interpc
    if (par.opt_interpc == on)
        interpc = par.interpc; linterpc = log(interpc);
        strt = length(p0) + 1;
        p0 = [p0; linterpc];
        par.pindx.linterpc = strt : length(p0);
    end 
    % d
    if (par.opt_d == on)
        d = par.d; ld = log(d);
        strt = length(p0) + 1;
        p0 = [p0; ld];
        par.pindx.ld = strt : length(p0);
    end 
    % kappa_dc
    if (par.opt_kappa_dc == on)
        kappa_dc = par.kappa_dc; lkappa_dc = log(kappa_dc);
        strt = length(p0) + 1;
        p0 = [p0; lkappa_dc];
        par.pindx.lkappa_dc = strt : length(p0);
    end 
    % RR
    if (par.opt_RR == on)
        RR = par.RR; lRR = log(RR);
        strt = length(p0) + 1;
        p0 = [p0; lRR];
        par.pindx.lRR = strt : length(p0);
    end
    % cc
    if (par.opt_cc == on)
        cc = par.cc; lcc = log(cc);
        strt = length(p0) + 1;
        p0 = [p0; lcc];
        par.pindx.lcc = strt : length(p0);
    end
    % dd
    if (par.opt_dd == on)
        dd = par.dd; ldd = log(dd);
        strt = length(p0) + 1;
        p0 = [p0; ldd];
        par.pindx.ldd = strt : length(p0);
    end
end
if par.Omodel == on
    % slopeo
    if (par.opt_slopeo == on)
        slopeo = par.slopeo; 
        strt = length(p0) + 1;
        p0 = [p0; slopeo];
        par.pindx.slopeo = strt : length(p0);
    end 
    % interpo
    if (par.opt_interpo == on)
        interpo = par.interpo; linterpo = log(interpo);
        strt = length(p0) + 1;
        p0 = [p0; linterpo];
        par.pindx.linterpo = strt : length(p0);
    end
end 
if par.Simodel == on 
    % bsi
    if (par.opt_bsi == on)
        bsi = par.bsi; lbsi = log(bsi);
        strt = length(p0) + 1;
        p0 = [p0; lbsi];
        par.pindx.lbsi = strt : length(p0);
    end
    % at 
    if (par.opt_at == on)
        at = par.at; lat = log(at);
        strt = length(p0) + 1;
        p0 = [p0; lat];
        par.pindx.lat = strt : length(p0);
    end
    % bt
    if (par.opt_bt == on)
        bt = par.bt; lbt = log(bt);
        strt = length(p0) + 1;
        p0 = [p0; lbt];
        par.pindx.lbt = strt : length(p0);
    end
    % aa
    if (par.opt_aa == on)
        aa = par.aa;
        strt = length(p0) + 1;
        p0 = [p0; aa];
        par.pindx.aa = strt : length(p0);
    end
    % bb
    if (par.opt_bb == on)
        bb = par.bb; lbb = log(bb);
        strt = length(p0) + 1;
        p0 = [p0; lbb];
        par.pindx.lbb = strt : length(p0);
    end
end
%
%% -------------------------------------------------------------
%
parm.kappa_p = 1/(720*60^2) ;
parm.p2c = 0.006+0.0069*po4obs;
parm.nzo = 2;
%%%%%%% prepare NPP for the model %%%%%%%%
inan = find(isnan(npp(:)) | npp(:) < 0);
npp(inan) = 0;

% parm.VER = '/DFS-L/DATA/primeau/weilewang/COP4WWF/MSK91/CTL_He_fixedPO_halfNP';
% tmp = squeeze(M3d(:,:,1));
% tmp(1:15,:) = nan; % SO
% tmp(65:78,55:125) = nan; % NP
% tmp(35:55,90:145) = nan; % EP
% iso = find(isnan(tmp));
% npp(iso) = npp(iso)*0.5;
%
parm.npp    = npp/(12*spd);
parm.Lambda = M3d*0;
parm.Lambda(:,:,1) = 0.5*(1/grd.dzt(1))*parm.p2c(:,:,1)./(1e-9+po4obs(:,:,1));
parm.Lambda(:,:,2) = 0.5*(1/grd.dzt(2))*parm.p2c(:,:,2)./(1e-9+po4obs(:,:,2));
parm.Lambda(:,:,3:end) = 0;

%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%
parm.p0 = p0;
x0 = p0;
%
myfun = @(x) neglogpost(x, parm, par);
%
options = optimoptions(@fminunc                  , ...
                       'Algorithm','trust-region', ...
                       'GradObj','on'            , ...
                       'Hessian','on'            , ...
                       'Display','iter'          , ...
                       'MaxFunEvals',2000        , ...
                       'MaxIter',2000            , ...
                       'TolX',5e-7               , ...
                       'TolFun',5e-7             , ...
                       'DerivativeCheck','off'   , ...
                       'FinDiffType','central'   , ...
                       'PrecondBandWidth',Inf)   ;
%
G_test = off;
nip    = length(x0);
if(G_test);
    dx = sqrt(-1)*eps.^3*eye(nip);
    for ii = 1:nip
        x  = real(x0)+dx(:,ii);
        [f,fx,fxx] = neglogpost(x, parm, par) ;
        diff = real(fx(ii)) - imag(f)/eps.^3 ;
        fprintf('%i %e  \n',ii,diff);
        diffx = real(fxx(:,ii)) - imag(fx)/eps.^3;
        for jj = 1:length(fx)
            fprintf('%e  ', diffx(jj));
        end
        fprintf('\n');
        exit 
    end
else
    [xhat,fval,exitflag] = fminunc(myfun,x0,options);
    [f,fx,fxx] = neglogpost(xhat,parm,par);
    fname = strcat(parm.VER,'_xhat.mat');
    load(fname)
    xhat.f   = f   ;
    xhat.fx  = fx  ;
    xhat.fxx = fxx ;
    save(fname, 'xhat')
end

fprintf('-------------- end! ---------------\n');