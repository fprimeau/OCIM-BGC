clc; clear all; close all
addpath('/DFS-L/DATA/primeau/weilewang/DATA')
addpath('/DFS-L/DATA/primeau/weilewang/my_func')
addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
on = true;      off = false;
spd  = 24*60^2; spa  = 365*spd;
%
TR_ver  = 90 ;
% mod_ver = 'CTL_He_noArc';
mod_ver = 'noArc';
%
Cmodel  = on ;
Omodel  = off ;
Simodel = off ;
% save results
% ATTENTION: please change this directory to where you wanna
if ismac
    input_dir = sprintf('~/Documents/CP-model/MSK%2d/',TR_ver); 
    % load optimal parameters if they exist
    fxhat = append(output_dir,mod_ver,'_xhat.mat');
elseif isunix 
    input_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
                        'MSK%2d/'],TR_ver);
    fxhat = append(input_dir,mod_ver,'_xhat.mat');
end
VER   = strcat(input_dir,mod_ver);
if (Cmodel == off & Omodel == off & Simodel == off)
    fname = strcat(VER,'_P');
elseif (Cmodel == on & Omodel == off & Simodel == off)
    fname = strcat(VER,'_PC');
elseif (Cmodel == on & Omodel == on & Simodel == off)
    fname = strcat(VER,'_PCO');
elseif (Cmodel == on & Omodel == off & Simodel == on)
    fname = strcat(VER,'_PCSi');
elseif (Cmodel == on & Omodel == on & Simodel == on)
    fname = strcat(VER,'_PCOSi');
end

if TR_ver == 90
    load transport_v4.mat
    load Sobs_90x180x24.mat     % woa2013 salinity data.
    load tempobs_90x180x24.mat
    load Siobs_90x180x24.mat Siobs
    load po4obs_90x180x24.mat
    grd  = grid; 
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
    load Sobs_91x180x24.mat     % woa2013 salinity data.
    load tempobs_91x180x24.mat
    load Siobs_91x180x24.mat Siobs
    load po4obs_91x180x24.mat
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
end
load(fname)
load(fxhat)
iwet = find(M3d(:));
nwet = length(iwet);
I    = speye(nwet);  
DIP  = DIP(iwet) ;
DIC  = DIC(iwet) ;
POC  = POC(iwet) ;
DOC  = DOC(iwet) ;
PO4  = po4obs(iwet) ;  
dAt  = grd.DXT3d.*grd.DYT3d;
dVt  = dAt.*grd.DZT3d;

dzt   = grd.dzt;
sigma = xhat.sigma ; 
% DOP remineralization rate constant. 
kappa4p = xhat.kappa_dp ; 
% linear parameter of npp to DIP assimilation function. 
alpha   = xhat.alpha ; 
% exponential parameter of npp to DIN assimilation function.
beta    = xhat.beta ; 
% POP disolution constant [s^-1];
kappa_p = 1/(720*60^2) ;   
% DOP remineralization rate constant. 
kappa4c = xhat.kappa_dc ;
cc      = xhat.cc ;
dd      = xhat.dd ;
bP      = xhat.bP ;
bP_T    = xhat.bP_T ;
bC      = xhat.bC   ;
bC_T    = xhat.bC_T ;
par.taup   = 720*60^2; % (s) pic dissolution time-scale
par.tau_TA = 1./par.taup;
par.M3d    = M3d   ;
par.grd    = grd   ;
par.iwet   = iwet  ;
par.nwet   = nwet  ;
par.TRdiv  = -TR   ;
par.dVt    = dVt   ;
par.Temp   = tempobs ;
par.Salt   = Sobs  ;

% PME part;
[modT,modS] = PME(par) ;
par.modT    = modT ;
par.modS    = modS ;
aveT = nanmean(modT(:,:,1:8),3);

nfig = 1;
figure(nfig)
bP2D = bP_T*aveT + bP ;
pcolor(bP2D); colorbar;shading flat 

if Cmodel == on
    nfig = nfig + 1;
    figure(nfig)
    bC2D = bC_T*aveT + bC ; 
    pcolor(bC2D); colorbar;shading flat 
    %
    nfig = nfig + 1;
    figure(nfig)
    C2P = M3d + nan ;
    C2P(iwet)  = 1./(cc*PO4 + dd);
    pcolor(C2P(:,:,1)); colorbar;shading flat
end

if Omodel == on
    nfig = nfig + 1;
    figure(nfig)
    par.slopeo  = xhat.slopeo  ;
    par.interpo = xhat.interpo ;
    vout  = mkO2P(par);
    O2P   = M3d + nan;
    O2P(iwet) = vout.O2P ;
    pcolor(nanmean(O2P(:,:,1:2),3)); colorbar;shading flat 
end

if Simodel == on
    nfig  = nfig + 1;
    aa    = xhat.aa ;
    bb    = xhat.bb ;
    Z     = Siobs(iwet);
    mu    = surface_mean(Z);
    Delta = sqrt(surface_mean((Z-mu).^2));
    
    % standardize the regressor variables
    ZR = 0.5+0.5*tanh((Z-mu)/Delta);
    %
    Si2C = (aa*ZR + bb)./Z;
    figure(nfig)
    pcolor(Si2C);colorbar;shading flat;
end 