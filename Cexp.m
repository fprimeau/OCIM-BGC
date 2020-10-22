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
par.optim   = on ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.LoadOpt = on ; % if load optimial par. 
par.pscale  = 0.0 ;
par.cscale  = 0.0 ; % factor to weigh DOC in the objective function
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
if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
    fname = strcat(VER,'_P');
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
    base_name = strcat(VER,'_PCv2');
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

load(par.fname)
load(par.fxhat)
%------------------- model grid info ----------------
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
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
I    = speye(nwet)  ;
TRdiv = -TR         ;
dAt   = grd.DXT3d.*grd.DYT3d ;
dVt   = dAt.*grd.DZT3d ;
par.M3d  = M3d  ;
par.grd  = grd  ;
par.iwet = iwet ;
par.nwet = nwet ;
%
%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPara(par) ;

%------------------ extract parameters --------------------
% POP disolution constant [s^-1];
sigma = par.sigma ; 
% linear parameter of npp to DIP assimilation function. 
alpha = par.alpha ; 
% exponential parameter of npp to DIN assimilation function.
beta  = par.beta ; 

%-------------------- normalize temperature --------------------
% Zscore is tried but it generate inf values(span from neg to pos).
% Tz   = (vT-mean(vT))./std(vT) ; 
vT = modT(iwet) ;
Tz = (vT - min(vT))./(max(vT) - min(vT)) ;
Tz( Tz < 0.05) = 0.05 ;

Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
Tz     = Tz*1e-8 ;
aveT   = nanmean(Tz3d(:,:,1:3),3) ;
par.aveT = aveT  ;
%------------------ prepare NPP for the model --------------------
nzo = 2 ;
p2c = 0.006+0.0069*po4obs ;
inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan) = 0 ;

npp   = npp/(12*spd) ;
npp1  = (0.5*npp./grd.dzt(1)).*p2c(:,:,1) ; 
npp2  = (0.5*npp./grd.dzt(2)).*p2c(:,:,2) ; 
Lambda = M3d*0 ;
Lambda(:,:,1) = 1./(1e-6+po4obs(:,:,1)) ;
Lambda(:,:,2) = 1./(1e-6+po4obs(:,:,2)) ;

% DIP assimilation
LAM        = 0*M3d;
LAM(:,:,1) = (npp1.^beta).*Lambda(:,:,1);
LAM(:,:,2) = (npp2.^beta).*Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];

%------------------ extract data ---------------------------
DIP  = data.DIP(iwet) ;
DIC  = data.DIC(iwet) ;
POC  = data.POC(iwet) ;
DOC  = data.DOC(iwet) ;
PO4  = po4obs(iwet)   ;  

% -------------- C:P uptake ratio -------------------------
W = d0(dVt(iwet)) ;
C2P3D = M3d + nan ;
C2P3D(iwet) = 1./(par.cc*PO4 + par.dd) ;
nn = 3 ;

%%%%%%%%%----------------%%%%%%%%%%%%%%
% calculate model primary production.
G        = M3d*0        ;
G(iwet)  = alpha*L*DIP  ; % primary production in P unit.
Int_CNPP = 0*M3d(:,:,1) ;
Int_PNPP = 0*M3d(:,:,1) ;

for ij = 1 : nn
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
kC    = d0(par.kC_T * Tz + par.kdC) ;
Jex_C = kC*d0(Prod)*(sigma*I+par.kappa_p*(1-sigma) * ...
        inv(F_diag_p+par.kappa_p*I))*((T_diag+kC*I)\Omega); 

C3d = M3d+nan;
C3d(iwet) = Jex_C;

% convert export from mmol C/m^3/s to mg C/m^2/day;
TOCexp = C3d(:,:,2).*grd.dzt(2)*12*spd;

tem_Cexp = TOCexp.*dAt(:,:,2);
Sum_Cexp = nansum(tem_Cexp(:))*365*1e-18;
fprintf('Model C export is %3.3e \n\n',Sum_Cexp);

% POC export
[~,Gout] = buildPFD(par,'POC') ;
w = -Gout.w ;
POCexp   = data.POC(:,:,2).*w(:,:,2)*12*spd ;
tem_POCexp = POCexp.*dAt(:,:,2);
Sum_POCexp = nansum(tem_POCexp(:))*365*1e-18;
fprintf('Model C export is %3.3e \n\n',Sum_POCexp);

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

