clc; clear alul; close all
global iter
iter = 0 ;
on   = true  ;
off  = false ;
format long
% 
GridVer  = 91  ;
operator = 'A' ;

Gtest = off ;
Htest = off ;
par.optim   = on ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.LoadOpt = on ; % if load optimial par. 
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigma = on ; 
par.opt_kP_T  = off ;
par.opt_kdP   = on ;
par.opt_bP_T  = off ; 
par.opt_bP    = on ;
par.opt_beta  = on ;
par.opt_alpha = on ;
% C model parameters
par.opt_bC_T  = off ;
par.opt_bC    = on ; 
par.opt_d     = on ;
par.opt_kC_T  = off ;
par.opt_kdC   = on ; 
par.opt_R_Si  = off ; 
par.opt_rR    = on ; 
par.opt_cc    = on ;
par.opt_dd    = on ;
% O model parameters
par.opt_O2C_T = off ;
par.opt_rO2C  = on ;
par.opt_O2P_T = off ; 
par.opt_rO2P  = on ; 
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;
%
%-------------load data and set up parameters---------------------
SetUp ;

% save results 
% ATTENTION: Change this direcrtory to where you wanna save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer); 
elseif isunix
    output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/Cexp/']);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
                        % 'MSK%2d/PME4DICALK/'],GridVer);
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
load(par.fname) ;
load(par.fxhat) ;

%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par) ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

%------------------ extract parameters ---------------------------
% POP disolution constant [s^-1];
sigma = par.sigma ; 
% linear parameter of npp to DIP assimilation function. 
alpha = par.alpha ; 
% exponential parameter of npp to DIN assimilation function.
beta  = par.beta ; 

%------------------ prepare NPP for the model --------------------
% DIP assimilation
LAM        = 0*M3d;
LAM(:,:,1) = (par.npp1.^beta).*par.Lambda(:,:,1);
LAM(:,:,2) = (par.npp2.^beta).*par.Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];

%------------------ extract data ---------------------------------
DIP  = data.DIP(iwet) ;
DIC  = data.DIC(iwet) ;
POC  = data.POC(iwet) ;
DOC  = data.DOC(iwet) ;
PO4  = po4obs(iwet)   ;  
TRdiv= par.TRdiv      ;
I    = par.I          ;
% -------------- C:P uptake ratio --------------------------------
W = d0(dVt(iwet)) ;
C2P3D = M3d + nan ;
C2P3D(iwet) = 1./(par.cc*PO4 + par.dd) ;
nn = 2 ;

%--------------- calculate primary production --------------------
G        = M3d*0        ;
G(iwet)  = alpha*L*DIP  ; % primary production in P unit.
Int_CNPP = 0*M3d(:,:,1) ;
Int_PNPP = 0*M3d(:,:,1) ;

for ij = 1 : nn
    Int_PNPP = Int_PNPP + G(:,:,ij).*grd.dzt(ij); 
    Int_CNPP = Int_CNPP + G(:,:,ij).*grd.dzt(ij).*C2P3D(:,:,ij)*12;
end
PNPP = Int_PNPP*spa*1e-3 ;
CNPP = Int_CNPP*spa*1e-3 ; % convert production from mg C/m^3/s to g
                           % C/m^2/year;
tem_CNPP = CNPP.*dAt(:,:,1)*1e-15 ;
Sum_CNPP = nansum(tem_CNPP(:))    ;
fprintf('Model NPP is %3.3e \n',Sum_CNPP) ;

%---------------- calculate phosphorus export --------------------
PFD = buildPFD(par, 'POP') ; 

F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;

junk = M3d ;
junk(:,:,1:nn) = 0 ;
Omega = junk(iwet) ;
Prod  = G(iwet)    ;
% adjoint method.
kP    = d0(par.kP_T * par.Tz + par.kdP) ;
Jex_P = d0(kP*Prod)*(sigma*I+par.kappa_p*(1-sigma) * ...
                     inv(F_diag_p+par.kappa_p*I))*((T_diag+kP)\Omega); 

P3d = M3d+nan;
P3d(iwet) = Jex_P;

%---------------- calculate carbon export -------------------------
PFD = buildPFD(par, 'POC') ; 

F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;

junk = M3d ;
junk(:,:,1:nn) = 0 ;
Omega = junk(iwet) ;
Prod  = G(iwet).*C2P3D(iwet) ;
% adjoint method.
kC    = d0(par.kC_T * par.Tz + par.kdC) ;
Jex_C = d0(kC*Prod)*(sigma*I+par.kappa_p*(1-sigma) * ...
                     inv(F_diag_p+par.kappa_p*I))*((T_diag+kC)\Omega); 

C3d = M3d+nan ;
C3d(iwet) = Jex_C ;

% convert export from mmol C/m^3/s to mg C/m^2/day;
TOCexp = C3d(:,:,2).*grd.dzt(2)*12*spd;

tem_Cexp = TOCexp.*dAt(:,:,2);
Sum_Cexp = nansum(tem_Cexp(:))*365*1e-18;
fprintf('Model TOC export is %3.3e \n\n',Sum_Cexp);

% POC export
[~,Gout] = buildPFD(par,'POC') ;
w = -Gout.w ;
POCexp   = data.POC(:,:,2).*w(:,:,2)*12*spd ;
tem_POCexp = POCexp.*dAt(:,:,2);
Sum_POCexp = nansum(tem_POCexp(:))*365*1e-18;
fprintf('Model POC export is %3.3e \n\n',Sum_POCexp);

%--------------- compare to ANCP -----------------------------
DOCexp = TOCexp - POCexp; DOCexp(DOCexp(:)<0) = 0 ;
TOCexp = smoothit(grd,M3d,TOCexp,3,1e5);
POCexp = smoothit(grd,M3d,POCexp,3,1e5);
DOCexp = smoothit(grd,M3d,DOCexp,3,1e5);

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

