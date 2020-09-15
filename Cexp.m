clc; clear all; close all
addpath('/DFS-L/DATA/primeau/weilewang/DATA')
addpath('/DFS-L/DATA/primeau/weilewang/my_func')
addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
spd  = 24*60^2; spa  = 365*spd;
%
TR_ver  = 91 ;
mod_ver = 'CTL_He_varP2O_noArc';
%
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
fname = strcat(VER,'_PCO');

if TR_ver == 90
    load transport_v4.mat
    load Sobs_90x180x24.mat     % woa2013 salinity data.
    load tempobs_90x180x24.mat
    load Mouw_POC_90x180x24.mat % sediment trap data MOUW
    load po4obs_90x180x24.mat
    load cbpm_npp_annual_90x180.mat
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
    load Mouw_POC_91x180x24.mat % sediment trap data MOUW
    load cbpm_npp_annual_91x180.mat
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
TRdiv = -TR;

DIP  = DIP(iwet) ;
DIC  = DIC(iwet) ;
POC  = POC(iwet) ;
DOC  = DOC(iwet) ;
PO4  = po4obs(iwet) ;  
dAt = grd.DXT3d.*grd.DYT3d;
dVt = dAt.*grd.DZT3d;

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
bC      = xhat.bC ;
bC_T    = xhat.bC_T ;
% (s) pic dissolution time-scale
par.taup   = 720*60^2 ; 
par.tau_TA = 1./par.taup ;
par.M3d    = M3d   ;
par.grd    = grd   ;
par.iwet   = iwet  ;
par.nwet   = nwet  ;
par.TRdiv  = TRdiv ;
par.dVt    = dVt   ;
par.Temp   = tempobs ;
par.Salt   = Sobs  ;
[modT,modS] = PME(par) ;
aveT   = nanmean(modT(:,:,1:8),3) ;
b = bC_T*aveT + bC ;
par.bC   = bC ; 
par.bC_T = bC_T ; 
par.aveT = aveT ;
par.kappa_p = kappa_p;

%%%%%%% prepare NPP for the model %%%%%%%%
nzo = 2;
p2c = 0.006 + 0.0069*po4obs;

inan = find(isnan(npp(:)) | npp(:) < 0);
npp(inan) = 0;
npp    = npp/(12*spd);
ismall = find(po4obs(iwet) < 1e-3);
% po4obs(iwet(ismall)) = 1e-3;
Lambda = M3d*0;
Lambda(:,:,1) = 0.5*(1/grd.dzt(1))*p2c(:,:,1)./(1e-9+po4obs(:,:,1));
Lambda(:,:,2) = 0.5*(1/grd.dzt(2))*p2c(:,:,2)./(1e-9+po4obs(:,:,2));
Lambda(:,:,3:end) = 0;

% DIP assimilation
Lambda(:,:,1) = (npp.^beta).*Lambda(:,:,1) ;
Lambda(:,:,2) = (npp.^beta).*Lambda(:,:,2) ;
L             = d0(Lambda(iwet))           ; % per second

% preparation for adjoint method.
W = d0(dVt(iwet)) ;
C2P3D = M3d + nan ;
C2P3D(iwet) = 1./(cc*PO4 + dd) ;
nn = 3 ;

%%%%%%%%%----------------%%%%%%%%%%%%%%
% calculate model primary production.
G        = M3d*0        ;
G(iwet)  = alpha*L*DIP  ; % primary production in P unit.
Int_CNPP = 0*M3d(:,:,1) ;
Int_PNPP = 0*M3d(:,:,1) ;

for ij = 1:nn
    Int_CNPP = Int_CNPP + G(:,:,ij).*grd.dzt(ij).*C2P3D(:,:,ij)*12;
    Int_PNPP = Int_PNPP + G(:,:,ij).*grd.dzt(ij); 
end
PNPP = Int_PNPP*spa*1e-3 ;
CNPP = Int_CNPP*spa*1e-3 ; % convert production from mg C/m^3/s to g
                           % C/m^2/year;
tem_CNPP = CNPP.*dAt(:,:,1)*1e-15 ;
Sum_CNPP = nansum(tem_CNPP(:))    ;
fprintf('Model NPP is %3.3e \n',Sum_CNPP) ;

%%%%%%%%% -------------- %%%%%%%%%%%%%%%
PFD = buildPFD(par, 'POC'); 
% calculate total export.
F_diag_p = inv(W)*PFD'*W;
T_diag   = inv(W)*TRdiv'*W;

junk = M3d;
junk(:,:,1:nn) = 0;
Omega = junk(iwet);
Prod  = G(iwet).*C2P3D(iwet);
% adjoint method.
Jex_C = kappa4c*d0(Prod)*(sigma*I+kappa_p*(1-sigma) * ...
        inv(F_diag_p+kappa_p*I))*((T_diag+kappa4c*I)\Omega); 

C3d = M3d+nan;
C3d(iwet) = Jex_C; 

Int_c = 0*M3d(:,:,1);
for ij = 1:nn
    Int_c = Int_c+C3d(:,:,ij).*grd.dzt(ij)*12;
end
% convert P export from mmol P/m^3/s to mg C/m^2/day;
TOCexp   = Int_c*spd; 
tem_Cexp = TOCexp.*dAt(:,:,3);
Sum_Cexp = nansum(tem_Cexp(:))*365*1e-18;
fprintf('Model C export is %3.3e \n\n',Sum_Cexp);
keyboard
%%%%%%%%%%%%%%%%%%% compare to ANCP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOCexp = smoothit(grd,M3d,TOCexp,3,1e5);
POCexp = smoothit(grd,M3d,POCexp,3,1e5);
DOCexp = TOCexp-POCexp;

Lat_HOTS = 22+45/60; Lon_HOTS = mod(-158,360);
Lat_BATS = 31+40/60; Lon_BATS = mod((-64-10/60),360);
Lat_OSP  = 50+1/60;  Lon_OSP  = mod((-144-9/60),360);

indx_hots = length(find(grd.xt<Lon_HOTS));
indy_hots = length(find(grd.yt<Lat_HOTS));

indx_bats = length(find(grd.xt<Lon_BATS));
indy_bats = length(find(grd.yt<Lat_BATS));

indx_osp = length(find(grd.xt<Lon_OSP));
indy_osp = length(find(grd.yt<Lat_OSP));

% find ANCP at specific location and convert unit from
% mg/m2/day to mol/m2/year;
TOCexp_HOTS = TOCexp(indy_hots,indx_hots)/12/1000*365;
TOCexp_BATS = TOCexp(indy_bats,indx_bats)/12/1000*365;
TOCexp_OSP  = TOCexp(indy_osp,indx_osp)/12/1000*365;
fprintf('TOC export at HOT is %2.2f mol/m2/year\n', TOCexp_HOTS)
fprintf('TOC export at BATS is %2.2f mol/m2/year \n',TOCexp_BATS)
fprintf('TOC export at OSP is %2.2f mol/m2/year \n\n', TOCexp_OSP)

msk_tropical = M3d(:,:,1)*0;
msk_tropical(length(find(grd.yt<-15)):length(find(grd.yt<15)),:) = 1;

junk1 = M3d(:,:,1)*0;
junk1(length(find(grd.yt<-30)):length(find(grd.yt<30)),:) = 1;
msk_subtro = junk1-msk_tropical;

junk2  = M3d(:,:,1)*0;
junk2(length(find(grd.yt<-45)):length(find(grd.yt<45)),:) = 1;
msk_subtro_subpo = junk2-junk1;

junk3 =  M3d(:,:,1)*0;
junk3(length(find(grd.yt<-60)):length(find(grd.yt<60)),:) = 1;
msk_subpolar = junk3-junk2;

% units mg/m2/day;
TOCexp_tropical = msk_tropical.*TOCexp;
TOCexp_subtro = msk_subtro.*TOCexp;
TOCexp_subtro_subpo = msk_subtro_subpo.*TOCexp;
TOCexp_subpolar = msk_subpolar.*TOCexp;

TOCexp_tropical(TOCexp_tropical(:)==0) = nan;
TOCexp_subtro(TOCexp_subtro(:)==0) = nan;
TOCexp_subtro_subpo(TOCexp_subtro_subpo(:)==0) = nan;
TOCexp_subpolar(TOCexp_subpolar(:)==0) = nan;

mean_TOC_tropical = nanmean(TOCexp_tropical(:))/12/1000*365;
mean_TOC_subtro = nanmean(TOCexp_subtro(:))/12/1000*365;
mean_TOC_subtro_subpo = nanmean(TOCexp_subtro_subpo(:))/12/1000*365;
mean_TOC_subpolar = nanmean(TOCexp_subpolar(:))/12/1000*365;
fprintf('mean TOC export at tropical is %2.2f mol/m2/yr\n', mean_TOC_tropical)
fprintf('TOC export at subtropical is %2.2f mol/m2/yr \n',mean_TOC_subtro)
fprintf('TOC export at subtropical-subpolar is %2.2f mol/m2/yr \n', mean_TOC_subtro_subpo)
fprintf('TOC export at subpolar is %2.2f mol/m2/year \n\n', mean_TOC_subpolar)

% calculate DOC to TOC export ratio for the four biomes.
D2T = DOCexp./TOCexp;
D2T_tropical = (msk_tropical.*D2T);
D2T_subtro = (msk_subtro.*D2T);
D2T_subtro_subpo = (msk_subtro_subpo.*D2T);
D2T_subpolar = (msk_subpolar.*D2T);

D2T_tropical(D2T_tropical(:)==0) = nan;
D2T_subtro(D2T_subtro(:)==0) = nan;
D2T_subtro_subpo(D2T_subtro_subpo(:)==0) = nan;
D2T_subpolar(D2T_subpolar(:)==0) = nan;

mean_D2T_tropical = nanmean(D2T_tropical(:));
mean_D2T_subtro = nanmean(D2T_subtro(:));
mean_D2T_subtro_subpo = nanmean(D2T_subtro_subpo(:));
mean_D2T_subpolar = nanmean(D2T_subpolar(:));

fprintf('tropial zonal mean DOC to TOC export ratio is %2.2f percent\n', ...
        mean_D2T_tropical*100)
fprintf('subtropical zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subtro*100)
fprintf('subtropical subpolar zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subtro_subpo*100)
fprintf('subpolar zonal mean DOC to TOC export ratio is %2.2f percent\n\n', ...
        mean_D2T_subpolar*100)

