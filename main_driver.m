clc; clear all; close all
on = true; off = false;
% addpath according to opterating system
if ismac 
    addpath('~/Dropbox/myfunc')
    addpath('~/Documents/DATA/')
    addpath('~/Documents/DATA/OCIM')
elseif isunix
    addpath('/DFS-L/DATA/primeau/weilewang/my_func/')
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/')
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
end

format long
global GC GO iter
iter = 0;
spd  = 24*60^2;
spa  = 365*spd;
%
Gtest = off ;
Htest = off ;
%
TR_ver  = 91 ;
mod_ver = 'CTL_He_varP2O_noArc';
par.optim   = on ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ; 
% save results 
% ATTENTION: please change this directory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',TR_ver); 
    % load optimal parameters if they exist
    fxhat = append(output_dir,mod_ver,'_xhat.mat');
elseif isunix
    output_dir = sprintf('/DFS-L/DATA/primeau/weilewang/COP4WWF/MSK%2d/',TR_ver);
    fxhat = append(output_dir,mod_ver,'_xhat.mat');
end
par.fxhat = fxhat;
VER = strcat(output_dir,mod_ver);
%
% load optimal parameters if they exist
% if isfile(fxhat)
% load(fxhat)
% end
%
% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_P');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_PC');
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
        fname = strcat(VER,'_PCO');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
        fname = strcat(VER,'_PCSi');
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
        fname = strcat(VER,'_PCOSi');
    end
end
par.fname = fname;
%
if TR_ver == 90
    load transport_v4.mat grid M3d TR
    load M3d90x180x24v2.mat MSKS 
    load Sobs_90x180x24.mat
    load tempobs_90x180x24.mat
    load po4obs_90x180x24.mat       % WOA PO4 observation
    load Siobs_90x180x24.mat Siobs
    %
    load DICant_90x180x24.mat
    load GLODAPv2_90x180x24raw.mat
    load splco2_mod_monthly.mat     % monthly CO2 data
    load co2syspar90.mat co2syspar
    load cbpm_npp_annual_90x180.mat
    load kw660_90x180.mat
    if ismac
        load MSK90/fixedPO_C.mat
        load MSK90/fixedPO_O2.mat
    elseif isunix
        load /DFS-L/DATA/primeau/weilewang/COP4WWF/MSK90/fixedPO_C.mat
        load /DFS-L/DATA/primeau/weilewang/COP4WWF/MSK90/fixedPO_O2.mat
    end 
    grd = grid ;
    
elseif TR_ver == 91
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
    load M3d91x180x24.mat MSKS 
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
    if ismac
        load MSK91/CTL_He_fixedPO_C.mat
        load MSK91/CTL_He_fixedPO_O2.mat
    elseif isunix
        load /DFS-L/DATA/primeau/weilewang/COP4WWF/MSK91/CTL_He_fixedPO_C.mat
        load /DFS-L/DATA/primeau/weilewang/COP4WWF/MSK91/CTL_He_fixedPO_O2.mat
    end
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
end
% get rid of arctice o2 observations
ARC = MSKS.ARC;
iarc = find(ARC(:));
o2raw(iarc)   = nan ;
dicraw(iarc)  = nan ;
po4raw(iarc)  = nan ;
sio4raw(iarc) = nan ; 
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
dVt  = grd.DXT3d.*grd.DYT3d.*grd.DZT3d;
%
[par.kw,par.P] = kw(M3d,grd);
par.Salt  = Sobs    ;
par.Temp  = tempobs ;
par.dVt   = dVt     ;
par.Kw660 = Kw660   ;
par.p4    = p4      ;
par.c2p   = 110     ;
par.M3d   = M3d     ;
par.iwet  = iwet    ;
par.nwet  = nwet    ;
par.TRdiv = -TR     ;
par.grd   = grd     ;
par.I     = speye(nwet);
par.rho   = 1024.5         ; % seawater density;
permil    = par.rho*1e-3  ; % from umol/kg to mmol/m3;
                            %
par.SIL     = Siobs;
par.po4obs  = po4obs;
par.o2raw   = o2raw;
par.po4raw  = po4raw;
par.sio4raw = sio4raw;
par.dicraw  = dicraw*permil; % GLODAP dic obs [mmol/m3];
par.human_co2 = DICant*permil;

GC = real([DIC(iwet); POC(iwet); DOC(iwet); CaC(iwet)]);
GO = real(O2(iwet)) + 1e-7*randn(par.nwet,1);

% transiant CO2 concentraion;
par.year      = splco2_mod(:,1) ;
par.pco2_air  = splco2_mod(:,2) ;
par.co2syspar = co2syspar       ;

par.kappa_g  = 1/(1e6*spa); % geological restoring time [1/s];
par.SILbar   = nansum(Siobs(iwet).*dVt(iwet))/nansum(dVt(iwet));
par.DIPbar   = nansum(po4obs(iwet).*dVt(iwet))/nansum(dVt(iwet)); % volume
par.taup     = 720*60^2; % (s) pic dissolution time-scale
par.tau_TA   = 1./par.taup;
par.kappa_gs = 1/(1e6*spa); % geological restoring time [1/s];
par.kappa_p  = 1/(720*60^2) ;
% PME part;
[modT,modS] = PME(par) ;
par.modS = modS        ;
par.modT = modT        ;
par.aveT = nanmean(modT(:,:,1:8),3);
%
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
if exist('xhat') & isfield(xhat,'bP_T')
    par.bP_T = xhat.bP_T;
else 
    par.bP_T = 0 ;
end 
if exist('xhat') & isfield(xhat,'bP')
    par.bP  = xhat.bP ;
else 
    par.bP  = 0.89 ;
end 
if exist('xhat') & isfield(xhat,'alpha')
    par.alpha = xhat.alpha;
else 
    par.alpha = 9.33e-08 ;
end 
if exist('xhat') & isfield(xhat,'beta')
    par.beta = xhat.beta;
else
    par.beta = 1.16e-01 ;
end 

% C model parameters                                      
if exist('xhat') & isfield(xhat,'bC_T')
    par.bC_T = xhat.bC_T ;
else
    par.bC_T = 0 ;
end 
if exist('xhat') & isfield(xhat,'bC')
    par.bC = xhat.bC ;
else
    par.bC = 1.06e+00 ;
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
if exist('xhat') & isfield(xhat,'dsi')
    par.dsi = xhat.dsi;
else
    par.dsi = 3300;
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
par.opt_beta  = on;
par.opt_alpha = on;
par.opt_sigma = on; 
par.opt_bP_T  = on; 
par.opt_bP    = on;
par.opt_kappa_dp = on;
% C model parameters
par.opt_d  = on; % 
par.opt_RR = on; % 
par.opt_cc = on;
par.opt_dd = on;
par.opt_bC_T = on;
par.opt_bC   = on; %
par.opt_kappa_dc = on; %
% O model parameters
par.opt_slopeo  = on; 
par.opt_interpo = on; 
% Si model parameters
par.opt_dsi = on;
par.opt_at = off;
par.opt_bt = on;
par.opt_aa = on;
par.opt_bb = on;
%% -------------------------------------------------------------
npx = 0; ncx = 0;
nox = 0; nsx = 0;
p0 = [];
% sigma
if (par.opt_sigma == on)
    npx = npx+1;
    sigma = par.sigma; lsigma = log(sigma);
    strt = length(p0) + 1;
    p0 = [p0; lsigma];
    par.pindx.lsigma = strt : length(p0);
end 
% kappa_dp 
if (par.opt_kappa_dp == on)
    npx = npx+1;
    kappa_dp = par.kappa_dp; lkappa_dp = log(kappa_dp);
    strt = length(p0) + 1;
    p0 = [p0; lkappa_dp];
    par.pindx.lkappa_dp = strt : length(p0);
end
% bP_T
if (par.opt_bP_T == on)
    npx = npx+1;
    bP_T = par.bP_T;
    strt = length(p0) + 1;
    p0 = [p0; bP_T];
    par.pindx.bP_T = strt : length(p0);
end 
% bP
if (par.opt_bP == on)
    npx = npx+1;
    bP = par.bP; lbP = log(bP);
    strt = length(p0) + 1;
    p0 = [p0; lbP];
    par.pindx.lbP = strt : length(p0);
end 
% alpha 
if (par.opt_alpha == on)
    npx = npx+1;
    alpha = par.alpha; lalpha = log(alpha);
    strt = length(p0) + 1;
    p0 = [p0; lalpha];
    par.pindx.lalpha = strt : length(p0);
end 
% beta
if (par.opt_beta == on)
    npx = npx+1;
    beta = par.beta; lbeta = log(beta);
    strt = length(p0) + 1;
    p0 = [p0; lbeta];
    par.pindx.lbeta = strt : length(p0);
end
%
if (par.Cmodel == on)
    % bC_T
    if (par.opt_bC_T == on)
        ncx = ncx + 1;
        bC_T = par.bC_T;
        strt = length(p0) + 1;
        p0 = [p0; bC_T];
        par.pindx.bC_T = strt : length(p0);
    end 
    % bC
    if (par.opt_bC == on)
        ncx = ncx + 1;
        bC = par.bC; lbC = log(bC);
        strt = length(p0) + 1;
        p0 = [p0; lbC];
        par.pindx.lbC = strt : length(p0);
    end 
    % d
    if (par.opt_d == on)
        ncx = ncx + 1;
        d = par.d; ld = log(d);
        strt = length(p0) + 1;
        p0 = [p0; ld];
        par.pindx.ld = strt : length(p0);
    end 
    % kappa_dc
    if (par.opt_kappa_dc == on)
        ncx = ncx + 1;
        kappa_dc = par.kappa_dc; lkappa_dc = log(kappa_dc);
        strt = length(p0) + 1;
        p0 = [p0; lkappa_dc];
        par.pindx.lkappa_dc = strt : length(p0);
    end 
    % RR
    if (par.opt_RR == on)
        ncx = ncx + 1;
        RR = par.RR; lRR = log(RR);
        strt = length(p0) + 1;
        p0 = [p0; lRR];
        par.pindx.lRR = strt : length(p0);
    end
    % cc
    if (par.opt_cc == on)
        ncx = ncx + 1;
        cc = par.cc; lcc = log(cc);
        strt = length(p0) + 1;
        p0 = [p0; lcc];
        par.pindx.lcc = strt : length(p0);
    end
    % dd
    if (par.opt_dd == on)
        ncx = ncx + 1;
        dd = par.dd; ldd = log(dd);
        strt = length(p0) + 1;
        p0 = [p0; ldd];
        par.pindx.ldd = strt : length(p0);
    end
    par.ncx = ncx;
end
%
if par.Omodel == on
    % slopeo
    if (par.opt_slopeo == on)
        nox = nox + 1;
        slopeo = par.slopeo; 
        strt = length(p0) + 1;
        p0 = [p0; slopeo];
        par.pindx.slopeo = strt : length(p0);
    end 
    % interpo
    if (par.opt_interpo == on)
        nox = nox + 1;
        interpo = par.interpo; linterpo = log(interpo);
        strt = length(p0) + 1;
        p0 = [p0; linterpo];
        par.pindx.linterpo = strt : length(p0);
    end
end 
if par.Simodel == on
    % dsi
    if (par.opt_dsi == on)
        nsx = nsx + 1;
        dsi = par.dsi; ldsi = log(dsi);
        strt = length(p0) + 1;
        p0 = [p0; ldsi];
        par.pindx.ldsi = strt : length(p0);
    end
    % at 
    if (par.opt_at == on)
        nsx = nsx + 1;
        at = par.at; lat = log(at);
        strt = length(p0) + 1;
        p0 = [p0; lat];
        par.pindx.lat = strt : length(p0);
    end
    % bt
    if (par.opt_bt == on)
        nsx = nsx + 1;
        bt = par.bt; lbt = log(bt);
        strt = length(p0) + 1;
        p0 = [p0; lbt];
        par.pindx.lbt = strt : length(p0);
    end
    % aa
    if (par.opt_aa == on)
        nsx = nsx + 1;
        aa = par.aa;
        strt = length(p0) + 1;
        p0 = [p0; aa];
        par.pindx.aa = strt : length(p0);
    end
    % bb
    if (par.opt_bb == on)
        nsx = nsx + 1;
        bb = par.bb; lbb = log(bb);
        strt = length(p0) + 1;
        p0 = [p0; lbb];
        par.pindx.lbb = strt : length(p0);
    end
end
par.npx = npx; par.ncx = ncx;
par.nox = nox; par.nsx = nsx;
%
%% -------------------------------------------------------------
%
par.p2c = 0.006+0.0069*po4obs;
par.nzo = 2;
%%%%%%% prepare NPP for the model %%%%%%%%
inan = find(isnan(npp(:)) | npp(:) < 0);
npp(inan) = 0;

% tmp = squeeze(M3d(:,:,1)) ;
% tmp(1:15,:) = nan         ; % SO
% tmp(65:78,55:125) = nan   ; % NP
% tmp(35:55,90:145) = nan   ; % EP
% itarg = find(isnan(tmp(:)))  ;
% npp(itarg) = npp(itarg)*0.5  ;
%
par.npp    = npp/(12*spd);
par.npp1   = (0.5*par.npp./grd.dzt(1)).*par.p2c(:,:,1) ; 
par.npp2   = (0.5*par.npp./grd.dzt(2)).*par.p2c(:,:,2) ; 
par.Lambda = M3d*0;
par.Lambda(:,:,1) = 1./(1e-6+po4obs(:,:,1));
par.Lambda(:,:,2) = 1./(1e-6+po4obs(:,:,2));
% par.Lambda(:,:,1) = 0.5*(1/grd.dzt(1))*par.p2c(:,:,1)./(0.1+po4obs(:,:,1));
% par.Lambda(:,:,2) = 0.5*(1/grd.dzt(2))*par.p2c(:,:,2)./(0.1+po4obs(:,:,2));
par.Lambda(:,:,3:end) = 0;
%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%
par.p0 = p0;
x0 = p0;
%
myfun = @(x) neglogpost(x, par);
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
nip    = length(x0);
if(Gtest);
    dx = sqrt(-1)*eps.^3*eye(nip);
    for ii = 1 : nip
        x  = real(x0)+dx(:,ii);
        if Htest == on
            [f,fx,fxx] = neglogpost(x, par) ;
            % print relative errors
            diff = (real(fx(ii)) - imag(f)/eps.^3)/(imag(f)/eps.^3);
            fprintf('%i % .3e  \n',ii,diff);
            diffx = (real(fxx(:,ii))-imag(fx)/eps.^3)./(imag(fx)/eps.^3+eps.^3);
            for jj = 1:length(fx)
                fprintf('% .3e  ', diffx(jj));
            end
        else
            [f,fx] = neglogpost(x, par) ;
            diff = (real(fx(ii)) - imag(f)/eps.^3)/(imag(f)/eps.^3) ;
            fprintf('%i % .3e  \n',ii,diff);
            fprintf('\n');
        end 
        fprintf('\n');
    end
    exit 
else
    [xhat,fval,exitflag] = fminunc(myfun,x0,options);
    [f,fx,fxx] = neglogpost(xhat,par);
    load(fxhat)
    xhat.f   = f   ;
    xhat.fx  = fx  ;
    xhat.fxx = fxx ;
    save(fxhat, 'xhat')
end

fprintf('-------------- end! ---------------\n');