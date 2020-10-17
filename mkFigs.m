clc; clear all; close all
on = true;      off = false;
spd  = 24*60^2; spa  = 365*spd;
% addpath according to opterating system
if ismac 
    addpath('~/Dropbox/myfunc'     )
    addpath('~/Documents/DATA/'    )
    addpath('~/Documents/DATA/OCIM')
else 
    addpath('/DFS-L/DATA/primeau/weilewang/DATA')
    addpath('/DFS-L/DATA/primeau/weilewang/my_func')
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
end 
format long
%
Cmodel  = on ;
Omodel  = off ;
Simodel = off ;
fscale  = 1.0 ;
%
GridVer   = 91 ;
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

% save results
% ATTENTION: please change this directory to where you wanna
if ismac
    input_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer); 
elseif isunix 
    input_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
                        'MSK%2d/'],GridVer);
end
VER   = strcat(input_dir,TRdivVer);
if (Cmodel == off & Omodel == off & Simodel == off)
    fname = strcat(VER,'_P');
elseif (Cmodel == on & Omodel == off & Simodel == off)
    base_name = strcat(VER,'_PCv1Zscore');
    catDOC = sprintf('_DOC%2.0e',fscale);
    fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == on & Simodel == off)
    base_name = strcat(VER,'_PCOv1');
    catDOC = sprintf('_DOC%2.0e',fscale);
    fname = strcat(base_name,catDOC);
elseif (Cmodel == on & par.Omodel == off & Simodel == on)
    base_name = strcat(VER,'_PCSi');
    catDOC = sprintf('_DOC%2.0e',fscale);
    fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == on & Simodel == on)
    base_name = strcat(VER,'_PCOSi');
    catDOC = sprintf('_DOC%2.0e',fscale);
    fname = strcat(base_name,catDOC);
end

% load optimal parameters if they exist
fxhat = strcat(fname,'_xhat.mat');
if GridVer == 90
    load transport_v4.mat
    load Sobs_90x180x24.mat     % woa2013 salinity data.
    load tempobs_90x180x24.mat
    load Siobs_90x180x24.mat Siobs
    load po4obs_90x180x24.mat
    load PME_TS_90x180x24.mat modT modS 
    grd  = grid; 
elseif GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    load Sobs_91x180x24.mat     % woa2013 salinity data.
    load tempobs_91x180x24.mat
    load po4obs_91x180x24.mat
    load PME_TS_91x180x24.mat modT modS
    M3d = output.M3d    ;
    grd = output.grid   ;
    TR  = output.TR/spa ;
end
%
load(fxhat)
iwet = find(M3d(:))   ;
nwet = length(iwet)   ;
I    = speye(nwet)    ;  
PO4  = po4obs(iwet)   ;  
dAt  = grd.DXT3d.*grd.DYT3d ;
dVt  = dAt.*grd.DZT3d ;
dzt  = grd.dzt;

kappa_p    = 1/(720*60^2)   ;   
par.taup   = 720*60^2    ; % (s) pic dissolution time-scale
par.tau_TA = 1./par.taup ;
par.tauPIC = (1/0.38)*spd;
par.M3d    = M3d     ;
par.grd    = grd     ;
par.iwet   = iwet    ;
par.nwet   = nwet    ;
par.TRdiv  = -TR     ;
par.dVt    = dVt     ;
par.Temp   = tempobs ;
par.Salt   = Sobs    ;
par.modT   = modT    ;
vT = modT(iwet) ;
% Tz = (vT - min(vT))./(max(vT) - min(vT)) ;
% Tz( Tz < 0.05 ) = 0.05 ;
Tz   = (vT-mean(vT))/std(vT) ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
aveT   = nanmean(Tz3d(:,:,1:3),3) ;

nfig = 0;
if isfield(xhat,'bP_T')
    nfig = nfig + 1;
    figure(nfig)
    bP   = xhat.bP   ;    
    bP_T = xhat.bP_T ;
    bP2D = bP_T*aveT + bP   ;
    pcolor(bP2D) ; colorbar ; shading flat 
    title('b4P')
    % saveas(gcf,'Figs91/b4P.png')
end
    
if isfield(xhat,'kP_T')
    nfig = nfig + 1  ;
    figure(nfig) 
    kP_T  = xhat.kP_T ;
    kdP   = xhat.kdP  ;
    
    kP3d = M3d + nan ;
    kP3d(iwet) = (1./(kP_T * Tz * 1e-8 + kdP))/spd ; 
    pcolor(kP3d(:,:,2));colorbar;shading flat
    % caxis([80 140])
    title('ka4P')
    % saveas(gcf,'Figs91/kappa4P.png')
end

if Cmodel == on
    % DOP remineralization rate constant. 
    if isfield(xhat,'bC_T')
        nfig = nfig + 1        ;
        figure(nfig)
        bC   = xhat.bC   ;
        bC_T = xhat.bC_T ;
        
        bC2D = bC_T*aveT + bC  ; 
        pcolor(bC2D); colorbar ; shading flat
        title('b4C')
        % saveas(gcf,'Figs91/b4C.png')
    end 

    if isfield(xhat,'kC_T') 
        nfig = nfig + 1  ;
        figure(nfig)
        kC_T = xhat.kC_T ;
        kdC  = xhat.kdC  ;
        kC3d = M3d + nan ;
        kC3d(iwet) = (1./(kC_T * Tz * 1e-8 + kdC))/spd; 
        pcolor(kC3d(:,:,2));colorbar;shading flat
        % caxis([100 600])
        title('kd4C')
        % saveas(gcf,'Figs91/kappa4C.png')
    end

    if isfield(xhat,'cc')
        nfig = nfig + 1  ;
        figure(nfig)
        cc   = xhat.cc   ;
        dd   = xhat.dd   ;
        
        C2P = M3d + nan  ;
        C2P(iwet)  = 1./(cc*PO4 + dd) ;
        pcolor(C2P(:,:,1)); colorbar;shading flat
        title('C:P uptake ratio')
        % saveas(gcf,'Figs91/CP ratio.png')
    end 
end

if Omodel == on
    if isfield(xhat,'O2C_T')
        nfig = nfig + 1    ;
        figure(nfig)
        O2C_T = xhat.O2C_T ;
        rO2C  = xhat.rO2C  ;
        O2C   = M3d + nan  ;
        O2C(iwet) = (O2C_T * Tz + rO2C) ; 
        pcolor(O2C(:,:,10)); colorbar; shading flat
        title('O2C consumption ratio')
    end 

    if isfield(xhat,'O2P_T')
        par.opt_O2P_T = on      ;
        nfig = nfig + 1         ;
        figure(nfig)
        par.O2P_T  = xhat.O2P_T ;
        par.rO2P   = xhat.rO2P  ;
        vout = mkO2P(par)       ;
        O2P  = M3d + nan        ;
        O2P(iwet) = vout.O2P    ;
        pcolor(nanmean(O2P(:,:,1:2),3)); colorbar; shading flat
        title('O2 production to P ratio')
    end 
end

if Simodel == on
    nfig  = nfig + 1 ;
    aa    = xhat.aa  ;
    bb    = xhat.bb  ;
    Z     = Siobs(iwet)     ;
    mu    = surface_mean(Z) ;
    Delta = sqrt(surface_mean((Z-mu).^2)) ;
    
    % standardize the regressor variables
    ZR = 0.5+0.5*tanh((Z-mu)/Delta) ;
    %
    Si2C = (aa*ZR + bb)./Z ;
    figure(nfig)
    pcolor(Si2C);colorbar; shading flat;
end 