clc; clear alul; close all
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
Gtest = off ;
Htest = off ;
par.optim   = on ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.LoadOpt = off ; % if load optimial par. 
par.pscale  = 0.0 ;
par.cscale  = 1.0 ; % factor to weigh DOC in the objective function
                    %
GridVer  = 91     ;
operator = 'A'    ;
if GridVer == 90
    TRdivVer = 'Tv4' ;
elseif GridVer == 91 
    switch(operator)
      case 'A'
        TRdivVer = 'CTL_He'             ;
      case 'B'
        TRdivVer = 'CTL_noHe'           ;
      case 'C'
        TRdivVer = 'KiHIGH_He'          ;
      case 'D'
        TRdivVer = 'KiHIGH_noHe'        ;
      case 'E'
        TRdivVer = 'KvHIGH_KiLOW_He'    ;
      case 'F'
        TRdivVer = 'KvHIGH_KiLOW_noHe'  ;
      case 'G'
        TRdivVer = 'KiLOW_He'           ;
      case 'H'
        TRdivVer = 'KiLOW_noHe'         ;
      case 'I'
        TRdivVer = 'KvHIGH_He'          ;
      case 'J'
        TRdivVer = 'KvHIGH_noHe'        ;
      case 'K'
        TRdivVer = 'KvHIGH_KiHIGH_noHe' ;
    end 
end 

fprintf('Transport version: % s \n', TRdivVer)
if par.Cmodel == on
    fprintf('---- C model is on ---- \n')
    fprintf('DOP scaling factor is %2.2e \n', par.pscale)
    fprintf('DOC scaling factor is %2.2e \n', par.cscale)
end
if par.Omodel == on
    fprintf('---- O model is on ---- \n')
end 
if par.Simodel == on
    fprintf('---- Si model is on ---- \n')
end 
fprintf('\n')

% P model parameters
par.opt_sigma = on ; 
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
par.opt_R_Si  = on ; 
par.opt_rR    = on ; 
par.opt_cc    = on ;
par.opt_dd    = on ;
% O model parameters
par.opt_O2C_T = on ;
par.opt_rO2C  = on ;
par.opt_O2P_T = on ; 
par.opt_rO2P  = on ; 
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = on ;
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
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
                        % 'MSK%2d/'],GridVer);
    output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
                        'MSK%2d/PME4DICALK/'],GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
                        % 'MSK%2d/'],GridVer);
end
VER = strcat(output_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_P');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
        base_name = strcat(VER,'_PCv1');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
        base_name = strcat(VER,'_PCOv1');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
        base_name = strcat(VER,'_PCSi');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
        base_name = strcat(VER,'_PCOSi');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    end
end
par.fname = strcat(fname,'.mat') ; 
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
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
    load GLODAPv2_talk.mat 
    load PME_TS_90x180x24.mat modT modS pme  
    load DICant_90x180x24.mat
    load GLODAPv2_90x180x24raw.mat
    load splco2_mod_monthly.mat     % monthly CO2 data
    load co2syspar90.mat co2syspar
    load cbpm_npp_annual_90x180.mat
    load DOMobs_90x180x24.mat
    load kw660_90x180.mat
    if ismac
        load MSK90/fixedPO_C.mat
        load MSK90/fixedPO_O2.mat
    elseif isunix
        if isfile(par.fname)
            load(par.fname)
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
    load GLODAPv2_talk_91x180x24.mat
    load PME_TS_91x180x24.mat modT modS pme 
    load DICant_91x180x24.mat
    load GLODAPv2_91x180x24raw.mat
    load splco2_mod_monthly % monthly CO2 data
    load co2syspar91.mat co2syspar
    load cbpm_npp_annual_91x180.mat
    load DOMobs_91x180x24.mat
    load kw660_91x180.mat
    if ismac
        load MSK91/CTL_He_PCO.mat
        load MSK91/CTL_He_PCO.mat
    elseif isunix
        if isfile(par.fname)
            load(par.fname)
        else 
            load initCO_91x180x24.mat
        end 
    end
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
end

%---------------------- constants -------------------
[par.kw,par.P] = kw(M3d,grd);
par.Kw660 = Kw660   ;
par.p4    = p4      ;
par.c2p   = 110     ;
par.rho   = 1024.5       ; % seawater density;
permil    = par.rho*1e-3 ; % from umol/kg to mmol/m3;
par.permil = permil ;
% transiant CO2 concentraion;
par.year      = splco2_mod(:,1) ;
par.pco2_air  = splco2_mod(:,2) ;
par.pco2atm   = splco2_mod(1,2) ;
par.co2syspar = co2syspar       ;

%------------------- model grid info ----------------
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
dVt  = grd.DXT3d.*grd.DYT3d.*grd.DZT3d;
ARC  = MSKS.ARC     ;
MED  = MSKS.MED     ;
PAC  = MSKS.PAC     ;
ATL  = MSKS.ATL     ;
IND  = MSKS.IND     ;
iarc = find(ARC(:)) ;
imed = find(MED(:)) ;
ipac = find(PAC(:)) ;
iatl = find(ATL(:)) ;
iind = find(IND(:)) ;
%
par.dVt   = dVt     ;
par.M3d   = M3d     ;
par.iwet  = iwet    ;
par.nwet  = nwet    ;
par.TRdiv = -TR     ;
par.grd   = grd     ;
par.MSKS  = MSKS    ;
par.I = speye(nwet) ;

%------------------------  data  -----------------------
% get rid of arctice o2 observations
% DOCobs(iarc)  = nan ;
% DOCobs(imed)  = nan ;
alkraw(iarc)  = nan ;
alkraw(imed)  = nan ; 
dicraw(iarc)  = nan ;
dicraw(imed)  = nan ;
% po4raw(iarc)  = nan ;
% po4raw(imed)  = nan ; 
% sio4raw(iarc) = nan ;
% sio4raw(imed) = nan ;
% o2raw(iarc)   = nan ;
% o2raw(imed)   = nan ;

par.Salt    = Sobs    ;
par.Temp    = tempobs ;
par.DSi     = Siobs   ;
par.po4obs  = po4obs  ;
par.o2raw   = o2raw   ;
par.po4raw  = po4raw  ;
par.sio4raw = sio4raw ;
par.DOCobs  = DOCobs  ;
par.alkraw  = alkraw*permil ;
par.dicraw  = dicraw*permil ; 
par.dicant  = DICant*permil ;
par.dopraw  = DOPobs - 0.03 ; % less refractory DOP 
DOCclean   = RemoveRef(par) ;
ibad = find( DOCclean(iarc) > 50 ) ;
DOCclean(iarc(ibad)) = nan ;
par.docraw = DOCclean ;

%---------------- inital guesses on C and O ---------------
DIC = data.DIC - par.dicant ;
GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); ...
       data.PIC(iwet); data.ALK(iwet)];
GC  = GC + 1e-6*randn(5*nwet,1) ;

%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPara(par) ;
% pack parameters into an array, assign them corresponding indices.
[p0, par] = PackPar(par) ;

%-------------------- prepare virtual flux -------------------
% PME part;
% [modT,modS,pme] = PME(par) ; 
par.modS    = modS ;
par.modT    = modT ;
par.pme     = pme  ;   
junk = M3d ;
junk(:,:,2:end) = 0 ;
isrf = find(junk(iwet)) ;
sdic = find(par.dicraw(iwet(isrf)) > 0);
salk = find(par.alkraw(iwet(isrf)) > 0);
% surface mean concentraions
par.sDICbar = sum(par.dicraw((iwet(isrf(sdic)))).* ...
                  dVt(iwet(isrf(sdic))))./sum(dVt(iwet(isrf(sdic))));
par.sALKbar = sum(par.alkraw((iwet(isrf(salk)))).* ...
                  dVt(iwet(isrf(salk))))./sum(dVt(iwet(isrf(salk))));

%-------------------- normalize temperature --------------------
% Zscore is tried but it generate inf values(span from neg to pos).
% Tz   = (vT-mean(vT))./std(vT) ; 
vT = par.modT(iwet) ;
Tz = (vT - min(vT))./(max(vT) - min(vT)) ;
Tz( Tz < 0.05) = 0.05 ;

Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par.Tz     = Tz*1e-8 ;
par.aveT   = nanmean(Tz3d(:,:,1:3),3) ;

%------------------- prepare for restoring ---------------------
% calculating global mean DIP, ALK, and DSi concentraions for
% restoring 
idip = find(par.po4raw(iwet)>0)  ;
ialk = find(par.alkraw(iwet)>0)  ;
isil = find(par.sio4raw(iwet)>0) ;

par.DIPbar = sum(par.po4raw(iwet(idip)).*dVt(iwet(idip)))/sum(dVt(iwet(idip))) ;
par.ALKbar = sum(par.alkraw(iwet(ialk)).*dVt(iwet(ialk)))/sum(dVt(iwet(ialk))) ;
par.DSibar = sum(par.sio4raw(iwet(isil)).*dVt(iwet(isil)))/sum(dVt(iwet(isil)));

%------------------ prepare NPP for the model --------------------
par.nzo = 2 ;
par.p2c = 0.006+0.0069*po4obs ;
inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan) = 0 ;

par.npp   = npp/(12*spd) ;
par.npp1  = (0.5*par.npp./grd.dzt(1)).*par.p2c(:,:,1) ; 
par.npp2  = (0.5*par.npp./grd.dzt(2)).*par.p2c(:,:,2) ; 
par.Lambda = M3d*0 ;
par.Lambda(:,:,1) = 1./(1e-6+po4obs(:,:,1)) ;
par.Lambda(:,:,2) = 1./(1e-6+po4obs(:,:,2)) ;
par.Lambda(:,:,3:end) = 0 ;

%-------------------set up fminunc -------------------------------
par.p0 = p0 ;
x0     = p0 ;
myfun = @(x) neglogpost(x, par);
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
    for ii = nip : -1 : 13
    % for ii = 8 : nip
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