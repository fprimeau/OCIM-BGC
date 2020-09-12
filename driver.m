clc; clear all; close all
global GC GO iter
iter = 0       ;
on   = true    ; off = false    ;
spd  = 24*60^2 ; spa  = 365*spd ;
% addpath according to opterating system
if ismac 
    addpath('~/Dropbox/myfunc'     )
    addpath('~/Documents/DATA/'    )
    addpath('~/Documents/DATA/OCIM')
elseif isunix
    addpath('/DFS-L/DATA/primeau/weilewang/my_func/'  )
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/'     )
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
end
format long
%
Gtest = on ;
Htest = on ;
par.optim   = on ; 
par.Cmodel  = on ; 
par.Omodel  = off ; 
par.Simodel = off ;
%
GridVer   = 90 ;
operator = 'A' ;
if GridVer == 90
    TRdivVer = 'Tv4' ;
elseif GridVer == 91 
    switch(operator)
      case 'A'
        TRdivVer = 'CTL_He'   ;
      case 'B'
        TRdivVer = 'CTL_noHe' ;
      case 'C'
        TRdivVer = 'KiHIGH_He'   ;
      case 'D'
        TRdivVer = 'KiHIGH_noHe' ;
      case 'E'
        TRdivVer = 'KvHIGH_KiLOW_He'  ;
      case 'F'
        TRdivVer = 'KvHIGH_KiLOW_noHe';
      case 'G'
        TRdivVer = 'KiLOW_He'   ;
      case 'H'
        TRdivVer = 'KiLOW_noHe' ;
      case 'I'
        TRdivVer = 'KvHIGH_He'  ;
      case 'J'
        TRdivVer = 'KvHIGH_noHe';
      case 'K'
        TRdivVer = 'KvHIGH_KiHIGH_noHe';
    end 
end 

fprintf('Transport version: % s \n', TRdivVer)
if par.Cmodel == on
    fprintf('---- C model is on ---- \n')
end
if par.Omodel == on
    fprintf('---- O model is on ---- \n')
end 
if par.Simodel == on
    fprintf('---- Si model is on ---- \n')
end 
fprintf('\n')
% P model parameters
par.opt_sigma = off ; 
par.opt_kP_T  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ; 
par.opt_bP    = on ;
par.opt_beta  = on ;
par.opt_alpha = on ;
% C model parameters
par.opt_bC_T  = on ;
par.opt_bC    = on ; 
par.opt_d     = on ;
par.opt_kC_T  = on ;
par.opt_kdC   = on ; 
par.opt_RR    = on ; 
par.opt_cc    = on ;
par.opt_dd    = on ;
% O model parameters
par.opt_O2C_T = on ;
par.opt_rO2C  = on ;
par.opt_O2P_T = on ; 
par.opt_rO2P  = on ; 
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;
%
% save results 
% ATTENTION: please change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer); 
elseif isunix
    output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
                        'MSK%2d/'],GridVer);
end
VER = strcat(output_dir,TRdivVer);
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
par.fname = fname ; 
% load optimal parameters if they exist
fxhat = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ; 
%
if GridVer == 90
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
        if isfile(fname)
            load(fname)
        else
            load initCO_90x180x24.mat
        end
    end 
    grd = grid ;
    
elseif GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    %
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
        load initCO_91x180x24.mat
    end
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
end

% get rid of arctice o2 observations
ARC  = MSKS.ARC;
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
par.I     = speye(nwet)  ;
par.rho   = 1024.5       ; % seawater density;
permil    = par.rho*1e-3 ; % from umol/kg to mmol/m3;

par.SIL     = Siobs   ;
par.po4obs  = po4obs  ;
par.o2raw   = o2raw   ;
par.po4raw  = po4raw  ;
par.sio4raw = sio4raw ;
par.dicraw  = dicraw*permil; % GLODAP dic obs [mmol/m3];
par.human_co2 = DICant*permil;

GC = real([DIC(iwet); POC(iwet); DOC(iwet); CaC(iwet)]);
GO = real(O2(iwet)) + 1e-9*randn(par.nwet,1);

% transiant CO2 concentraion;
par.year      = splco2_mod(:,1) ;
par.pco2_air  = splco2_mod(:,2) ;
par.co2syspar = co2syspar       ;

par.SILbar = nansum(Siobs(iwet).*dVt(iwet))/nansum(dVt(iwet))  ;
par.DIPbar = nansum(po4obs(iwet).*dVt(iwet))/nansum(dVt(iwet)) ;

% load optimal parameters from a file
% or set them to default values 
par = SetPara(par) ;
%
% pack adjustable parameters in an array and
% assign them corresponding indices.
[p0, par] = PackPar(par) ;
%
% PME part;
[modT,modS] = PME(par)   ;
par.modS    = modS       ;
par.modT    = modT       ;
Tz  = zscore(modT(iwet)) ;
Tz3d        = M3d + nan  ;
Tz3d(iwet)  = Tz         ;
par.aveT    = nanmean(Tz3d(:,:,1:3),3) ;
par.Tz      = Tz*1e-8    ;
%
%%%%%%% prepare NPP for the model %%%%%%%%
par.nzo = 2 ;
par.p2c = 0.006+0.0069*po4obs ;
inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan) = 0 ;

% tmp = squeeze(M3d(:,:,1)) ;
% tmp(1:15,:) = nan         ; % SO
% tmp(65:78,55:125) = nan   ; % NP
% tmp(35:55,90:145) = nan   ; % EP
% itarg = find(isnan(tmp(:)))  ;
% npp(itarg) = npp(itarg)*0.5  ;
%
par.npp    = npp/(12*spd) ;
par.npp1   = (0.5*par.npp./grd.dzt(1)).*par.p2c(:,:,1) ; 
par.npp2   = (0.5*par.npp./grd.dzt(2)).*par.p2c(:,:,2) ; 
par.Lambda = M3d*0 ;
par.Lambda(:,:,1) = 1./(1e-6+po4obs(:,:,1)) ;
par.Lambda(:,:,2) = 1./(1e-6+po4obs(:,:,2)) ;

par.Lambda(:,:,3:end) = 0 ;
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
nip = length(x0);
if(Gtest);
    dx = sqrt(-1)*eps.^3*eye(nip);
    for ii = 6 : nip
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
    keyboard 
else
    [xhat,fval,exitflag] = fminunc(myfun,x0,options);
    [f,fx,fxx,data] = neglogpost(xhat,par);
    load(fxhat)
    xhat.f   = f   ;
    xhat.fx  = fx  ;
    xhat.fxx = fxx ;
    % save results 
    save(fxhat, 'xhat')
    save(par.fname, 'data')
end

fprintf('-------------- end! ---------------\n');