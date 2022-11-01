%Cexp_future
%Cexp_Tup5C

% Use adjoint to compute C and P export for optimized model and future
% projections.

% or compute C export for each model run and then save model output structure with
% another field for Cexp and Pexp
%or should i put something like this inside driver_future??
% or inside the normal Cexp

% Plot Optimized zonal average uptake (solid) and export (dashed lines)

clc; clear all; close all
on   = true  ;
off  = false ;
spd  = 24*60^2 ;
spa  = 365*spd ;


GridVer  = 91  ;
operator = 'A' ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function
par.dynamicP = off;

cd ../

%% load T only model output
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/';
figDir = strcat(outputDir,'futureproj/v2_');


% load Cell model output
fname = strcat(outputDir,'optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0.mat');
load(fname);
model_Cell = data;
fxhat = strcat(outputDir,'optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_Cell = xhat.allparams;
clear data xhat

% load P decrease by 20% Cell output fields
fname_future = strcat(outputDir, 'futureproj/Pdn20pct_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_Cell_Pdn20 = d.model;
clear d

% load T increase by 5C Cell output fields
fname_future = strcat(outputDir, 'futureproj/Tup5Cv2_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_Cell_Tup5C = d.model;
clear d

% load CMIP5projection to 2100 Cell output fields
fname_future = strcat(outputDir, 'futureproj/CMIP2100_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_Cell_CMIP = d.model;
clear d


%% ------ set up figures ---------

set(groot,'defaultAxesFontName','Times',...
    'defaultAxesFontSize',14,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','on',...
    'defaultAxesYMinorTick','on');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Times',...
    'defaultTextInterpreter','latex');

% Define Some Colors ROYGBIV

colors.maroon 		= [128/255 0 0];
colors.tomato 		= [255/255 99/255 71/255];	% light red-orange
colors.indianred 	= [205/255 92/255 92/255]; % light red-brown
colors.limegreen 	= [50/255 205/255 50/255];
colors.darkgreen 	= [0 100/255 0];
colors.teal 		= [0 128/255 128/255];
colors.aqua 		= [0.2 0.8 0.8];
colors.lblue 		= [0 191/255 255/255];
colors.navy 		= [ 0 0 128/255];
colors.darkmagenta 	= [139/255 0 139/255];

%% -------------load data and set up parameters---------------------

% load grid
addpath('/DFS-L/DATA/primeau/weilewang/DATA/')
addpath('/DFS-L/DATA/primeau/weilewang/my_func/'  )

load M3d91x180x24.mat MSKS
OperName = sprintf('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2/OCIM2_CTL_He.mat');
    load(OperName,'output') ;
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
    %
iwet = find(M3d(:)) ;
dAt  = grd.DXT3d.*grd.DYT3d;
dVt  = dAt.*grd.DZT3d;
I = speye(length(iwet)) ;

%volume weight diagonal matrix
W = d0(dVt(iwet)) ;

par.M3d = M3d;
par.grd = grd;
par.TRdiv = -TR     ;
TRdiv= par.TRdiv      ;

% number of euphotic zone layers
nn = 2;
par.nzo = nn;

% ---------- load phosphate -----------------
load po4obs_91x180x24.mat % WOA PO4 observation
par.po4obs  = po4obs  ;
%PO4  = par.po4obs(iwet)   ;

%-------------- load and set up NPP -----------------
load cbpm_npp_annual_91x180.mat

inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan) = 0 ;

par.npp   = npp/(12*spd) ;		% units: mmol C/m^2/s /(1 mmol C/m^2/s)
par.npp1  = (0.5*par.npp./grd.dzt(1)) ; % units: mmol C/m^2/s
par.npp2  = (0.5*par.npp./grd.dzt(2)) ;
par.Lambda = M3d*0 ;
par.Lambda(:,:,1) = 1./(1e-6+po4obs(:,:,1)) ;
par.Lambda(:,:,2) = 1./(1e-6+po4obs(:,:,2)) ;
par.Lambda(:,:,3:end) = 0 ;

% ----------- load and set up temperature ----------------
load tempobs_91x180x24.mat
tempobs(tempobs(:)<-2.0) = -2.0 ;
tempobs = tempobs;
par.Temp    = tempobs ;

% --- normalize temperature --
for ji = 1:24
    t2d = par.Temp(:,:,ji);
    par.Temp(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
end
vT = par.Temp(iwet) ;
Tz = (vT - min(vT))./(max(vT) - min(vT)) ;
% Tz = zscore(vT)  ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par.Tz     = Tz*1e-8 ;
par.aveT   = nanmean(Tz3d(:,:,1:2),3) ;
clear t2d Tz Tz3d

%% -------- new par structure for Tup5C ------------
% get original min and Max temperature for Normalization
vT0min = min(vT);
vT0max = max(vT);
clear vT

par_Tup5C = par;
% ----overwrite temperature obs from SetUp -----
	load tempobs_91x180x24.mat
tempobs(tempobs(:)<-2.0) = -2.0 ;
% Uniform  +5 degree C temperature increase
tempobs = tempobs+5;
par_Tup5C.Temp_proj    = tempobs ;

%-- normalize temperature --
for ji = 1:24
    t2d = par_Tup5C.Temp_proj(:,:,ji);
    par_Tup5C.Temp_proj(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
end
vT = par_Tup5C.Temp_proj(iwet) ;
Tz = (vT - vT0min)./(vT0max - vT0min) ;
% Tz = zscore(vT)  ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par_Tup5C.Tz_proj     = Tz*1e-8 ;
par_Tup5C.aveT_proj   = nanmean(Tz3d(:,:,1:2),3) ;
clear t2d vT Tz Tz3d tempobs


%% Cell model


%--------------- calculate primary production --------------------
% DIP assimilation; same for cell optim and Tup5C
LAM        = 0*M3d;
LAM(:,:,1) = (par.npp1.^xhat_Cell.beta).*par.Lambda(:,:,1);
LAM(:,:,2) = (par.npp2.^xhat_Cell.beta).*par.Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
L_Cell = L;
clear LAM

%% Cell model opt compute C export using adjoint
G_Cell        = M3d*0        ;
G_Cell(iwet)  = xhat_Cell.alpha*L_Cell*model_Cell.DIP(iwet)  ; % primary production [unit: mmol P/m^3/s]

%---------------- calculate phosphorus export --------------------
par.bP = xhat_Cell.bP;
par.bP_T = xhat_Cell.bP_T;
par.bC = xhat_Cell.bC;
par.bC_T = xhat_Cell.bC_T;
par.kappa_p = xhat_Cell.kappa_p;

PFD = buildPFD(par, 'POP') ;
F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;

junk = M3d ;
junk(:,:,1:2) = 0 ;
Omega = junk(iwet) ;
% adjoint method.
kP    = d0(xhat_Cell.kP_T * par.Tz + xhat_Cell.kdP) ;
Jex_P = d0(kP*G_Cell(iwet))*(xhat_Cell.sigma*I+xhat_Cell.kappa_p*(1-xhat_Cell.sigma) * ...
                     inv(F_diag_p+xhat_Cell.kappa_p*I))*((T_diag+kP)\Omega);

P3d_Cell = M3d+nan;
P3d_Cell(iwet) = Jex_P;
clear kP Jex_P PFD F_diag_p T_diag

% convert export from mmol P/m^3/s to mg P/m^2/day;
tem_Pexp = (P3d_Cell(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31*spd).*dAt(:,:,1:nn);			% [mg P/day]
fprintf('Cell Model (optimal): TOP export is %3.4f Pg P /yr   (Integrated to %4.1f m) \n',nansum(tem_Pexp(:))*365*1e-18, sum(grd.dzt(1:nn)));

% convert export from mmol P/m^3/s to mol P/m^2/yr;
TOPexp3d_Cell_opt = P3d_Cell(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
TOPexp_Cell_opt = sum(TOPexp3d_Cell_opt,3); %sum(TOPexp3d,3,'omitnan');

%---------------- calculate carbon export -------------------------
	PFD_C = buildPFD(par, 'POC') ;
	F_diag_p = inv(W)*PFD_C'*W   ;
	T_diag   = inv(W)*TRdiv'*W ;

	C2P_Cell = model_Cell.CellOut.C2P(iwet);
	Prod_Cell  = G_Cell(iwet).*C2P_Cell ;
	% adjoint method.
	kC    = d0(xhat_Cell.kC_T * par.Tz + xhat_Cell.kdC) ;
	Jex_C = d0(kC*Prod_Cell)*(xhat_Cell.sigma*I+xhat_Cell.kappa_p*(1-xhat_Cell.sigma) * ...
	                     inv(F_diag_p+xhat_Cell.kappa_p*I))*((T_diag+kC)\Omega);

	C3d_Cell = M3d+nan ;
	C3d_Cell(iwet) = Jex_C ;
    clear kC Jex_C PFD_C F_diag_p T_diag

    % convert export from mmol P/m^3/s to mol P/m^2/yr;
    TOCexp3d_Cell_opt = C3d_Cell(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
    TOCexp_Cell_opt = sum(TOCexp3d_Cell_opt,3); %sum(TOPexp3d,3,'omitnan');

	% convert export from mmol C/m^3/s to mg C/m^2/day;
	tem_Cexp = (C3d_Cell(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*12*spd).*dAt(:,:,1:nn);
	fprintf('Cell Model (optimal): TOC export is %3.4f Pg C /yr   (Integrated to %4.1f m) \n',nansum(tem_Cexp(:))*365*1e-18, sum(grd.dzt(1:nn)));

% save
    EXPORT.TOPexp_Cell_opt = TOPexp_Cell_opt;
    EXPORT.TOCexp_Cell_opt = TOCexp_Cell_opt;
    %save('export_Cell_Tup5C.mat','TOPexp3d_Cell_Tup5C','TOCexp3d_Cell_Tup5C')


%%------------- Cell_Tup5C compute C export using adjoint -----------------------
G_Cell_Tup5C        = M3d*0        ;
G_Cell_Tup5C(iwet)  = xhat_Cell.alpha*L_Cell*model_Cell_Tup5C.DIP(iwet)  ; % primary production [unit: mmol P/m^3/s]

%---------------- calculate phosphorus export --------------------
par_Tup5C.bP = xhat_Cell.bP;
par_Tup5C.bP_T = xhat_Cell.bP_T;
par_Tup5C.bC = xhat_Cell.bC;
par_Tup5C.bC_T = xhat_Cell.bC_T;
par_Tup5C.kappa_p = xhat_Cell.kappa_p;

PFD = buildPFD(par_Tup5C, 'POP') ;

F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;

junk = M3d ;
junk(:,:,1:2) = 0 ;
Omega = junk(iwet) ;
% adjoint method.
kP    = d0(xhat_Cell.kP_T * par_Tup5C.Tz + xhat_Cell.kdP) ;
Jex_P = d0(kP*G_Cell_Tup5C(iwet))*(xhat_Cell.sigma*I+xhat_Cell.kappa_p*(1-xhat_Cell.sigma) * ...
                     inv(F_diag_p+xhat_Cell.kappa_p*I))*((T_diag+kP)\Omega);

P3d_Cell_Tup5C = M3d+nan;
P3d_Cell_Tup5C(iwet) = Jex_P;
clear kP Jex_P

% convert export from mmol P/m^3/s to mol P/m^2/yr;
TOPexp3d_Cell_Tup5C = P3d_Cell_Tup5C(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
TOPexp_Cell_Tup5C = sum(TOPexp3d_Cell_Tup5C,3); %sum(TOPexp3d,3,'omitnan');

% convert export from mmol P/m^3/s to mg P/m^2/day;
tem_Pexp = (P3d_Cell_Tup5C(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31*spd).*dAt(:,:,1:nn);			% [mg P/day]
fprintf('Cell Model (Temp +5C): TOP export is %3.4f Pg P /yr   (Integrated to %4.1f m) \n',nansum(tem_Pexp(:))*365*1e-18, sum(grd.dzt(1:nn)));

%---------------- calculate carbon export -------------------------
	PFD_C = buildPFD(par_Tup5C, 'POC') ;
	F_diag_p = inv(W)*PFD_C'*W   ;
	T_diag   = inv(W)*TRdiv'*W ;

	C2P_Cell_Tup5C = model_Cell_Tup5C.CellOut.C2P(iwet);
	Prod_Cell_Tup5C  = G_Cell_Tup5C(iwet).*C2P_Cell_Tup5C ;
	% adjoint method.
	kC    = d0(xhat_Cell.kC_T * par_Tup5C.Tz + xhat_Cell.kdC) ;
	Jex_C = d0(kC*Prod_Cell_Tup5C)*(xhat_Cell.sigma*I+xhat_Cell.kappa_p*(1-xhat_Cell.sigma) * ...
	                     inv(F_diag_p+xhat_Cell.kappa_p*I))*((T_diag+kC)\Omega);

	C3d_Cell_Tup5C = M3d+nan ;
	C3d_Cell_Tup5C(iwet) = Jex_C ;
    clear kC Jex_C

    % convert export from mmol P/m^3/s to mol P/m^2/yr;
    TOCexp3d_Cell_Tup5C = C3d_Cell_Tup5C(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
    TOCexp_Cell_Tup5C = sum(TOCexp3d_Cell_Tup5C,3); %sum(TOPexp3d,3,'omitnan');

	% convert export from mmol C/m^3/s to mg C/m^2/day;
	tem_Cexp = (C3d_Cell_Tup5C(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*12*spd).*dAt(:,:,1:nn);
	fprintf('Cell Model (Temp +5C): TOC export is %3.4f Pg C /yr   (Integrated to %4.1f m) \n',nansum(tem_Cexp(:))*365*1e-18, sum(grd.dzt(1:nn)));

% save
    EXPORT.TOPexp_Cell_Tup5C = TOPexp_Cell_Tup5C;
    EXPORT.TOCexp_Cell_Tup5C = TOCexp_Cell_Tup5C;
    %save('export_Cell_Tup5C.mat','TOPexp3d_Cell_Tup5C','TOCexp3d_Cell_Tup5C')

clear par_Tup5C

%% ----------- Pdn20pct --------------------------------------
%%------------ new par structure for Pdn20pct --------------

par_Pdn20 = par;
iprod = find(M3d(:,:,1:2));
po4obs(iprod) = po4obs(iprod).*0.8;

% Uniform  20 percent Phosphate decrease in euphotic zone
par_Pdn20.po4obs  = po4obs  ;

clear po4obs

%%------------- Cell_Pdn20pct compute C export using adjoint -----------------------
G_Cell_Pdn20        = M3d*0        ;
G_Cell_Pdn20(iwet)  = xhat_Cell.alpha*L_Cell*model_Cell_Pdn20.DIP(iwet)  ; % primary production [unit: mmol P/m^3/s]

%---------------- calculate phosphorus export --------------------
par_Pdn20.bP = xhat_Cell.bP;
par_Pdn20.bP_T = xhat_Cell.bP_T;
par_Pdn20.bC = xhat_Cell.bC;
par_Pdn20.bC_T = xhat_Cell.bC_T;
par_Pdn20.kappa_p = xhat_Cell.kappa_p;

PFD = buildPFD(par_Pdn20, 'POP') ;

F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;

junk = M3d ;
junk(:,:,1:2) = 0 ;
Omega = junk(iwet) ;
% adjoint method.
kP    = d0(xhat_Cell.kP_T * par_Pdn20.Tz + xhat_Cell.kdP) ;
Jex_P = d0(kP*G_Cell_Pdn20(iwet))*(xhat_Cell.sigma*I+xhat_Cell.kappa_p*(1-xhat_Cell.sigma) * ...
                     inv(F_diag_p+xhat_Cell.kappa_p*I))*((T_diag+kP)\Omega);

P3d_Cell_Pdn20 = M3d+nan;
P3d_Cell_Pdn20(iwet) = Jex_P;
clear kP Jex_P

% convert export from mmol P/m^3/s to mol P/m^2/yr;
TOPexp3d_Cell_Pdn20 = P3d_Cell_Pdn20(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
TOPexp_Cell_Pdn20 = sum(TOPexp3d_Cell_Pdn20,3); %sum(TOPexp3d,3,'omitnan');

% convert export from mmol P/m^3/s to mg P/m^2/day;
tem_Pexp = (P3d_Cell_Pdn20(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31*spd).*dAt(:,:,1:nn);			% [mg P/day]
fprintf('Cell Model (PO4 -20%%): TOP export is %3.4f Pg P /yr   (Integrated to %4.1f m) \n',nansum(tem_Pexp(:))*365*1e-18, sum(grd.dzt(1:nn)));

%---------------- calculate carbon export -------------------------
	PFD_C = buildPFD(par_Pdn20, 'POC') ;
	F_diag_p = inv(W)*PFD_C'*W   ;
	T_diag   = inv(W)*TRdiv'*W ;

	C2P_Cell_Pdn20 = model_Cell_Pdn20.CellOut.C2P(iwet);
	Prod_Cell_Pdn20  = G_Cell_Pdn20(iwet).*C2P_Cell_Pdn20 ;
	% adjoint method.
	kC    = d0(xhat_Cell.kC_T * par_Pdn20.Tz + xhat_Cell.kdC) ;
	Jex_C = d0(kC*Prod_Cell_Pdn20)*(xhat_Cell.sigma*I+xhat_Cell.kappa_p*(1-xhat_Cell.sigma) * ...
	                     inv(F_diag_p+xhat_Cell.kappa_p*I))*((T_diag+kC)\Omega);

	C3d_Cell_Pdn20 = M3d+nan ;
	C3d_Cell_Pdn20(iwet) = Jex_C ;
    clear kC Jex_C

    % convert export from mmol P/m^3/s to mol P/m^2/yr;
    TOCexp3d_Cell_Pdn20 = C3d_Cell_Pdn20(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
    TOCexp_Cell_Pdn20 = sum(TOCexp3d_Cell_Pdn20,3); %sum(TOPexp3d,3,'omitnan');

	% convert export from mmol C/m^3/s to mg C/m^2/day;
	tem_Cexp = (C3d_Cell_Pdn20(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*12*spd).*dAt(:,:,1:nn);
	fprintf('Cell Model (PO4 -20%%): TOC export is %3.4f Pg C /yr   (Integrated to %4.1f m) \n',nansum(tem_Cexp(:))*365*1e-18, sum(grd.dzt(1:nn)));

% save
    EXPORT.TOPexp_Cell_Pdn20 = TOPexp_Cell_Pdn20;
    EXPORT.TOCexp_Cell_Pdn20 = TOCexp_Cell_Pdn20;
    %save('export_Cell_Pdn20.mat','TOPexp3d_Cell_Pdn20','TOCexp3d_Cell_Pdn20')

clear par_Pdn20

%%--- CMIP -------------------------------

%% ------- new par stucture for CMIP2100 --------
par_CMIP = par;
if GridVer ==91
	load('/DFS-L/DATA/primeau/meganrs/DATA/CMIP5/CMIP5mean_no3_thetao_91x180x24.mat')
end

% overwrite po4obs
par_CMIP.po4obs = PO4_CMIP5mean;
par_CMIP.po4proj  = PO4_CMIP5mean  ;

% get original min and Max temperature for Normalization
vT0 = par.Temp(iwet) ;
vT0min = min(vT0);
vT0max = max(vT0);

% ----overwrite temperature obs from SetUp -----
tempobs = T_CMIP5mean;
tempobs(tempobs(:)<-2.0) = -2.0 ;
% Uniform  +5 degree C temperature increase
par_CMIP.Temp_proj    = tempobs ;

vT = par_CMIP.Temp_proj(iwet) ;
Tz = (vT - vT0min)./(vT0max - vT0min) ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par_CMIP.Tz_proj     = Tz*1e-8 ;
par_CMIP.aveT_proj   = nanmean(Tz3d(:,:,1:2),3) ;

clear t2d vT Tz Tz3d tempobs
%%------------- Cell_CMIP compute C export using adjoint -----------------------
G_Cell_CMIP        = M3d*0        ;
G_Cell_CMIP(iwet)  = xhat_Cell.alpha*L_Cell*model_Cell_CMIP.DIP(iwet)  ; % primary production [unit: mmol P/m^3/s]

%---------------- calculate phosphorus export --------------------
par_CMIP.bP = xhat_Cell.bP;
par_CMIP.bP_T = xhat_Cell.bP_T;
par_CMIP.bC = xhat_Cell.bC;
par_CMIP.bC_T = xhat_Cell.bC_T;
par_CMIP.kappa_p = xhat_Cell.kappa_p;

PFD = buildPFD(par_CMIP, 'POP') ;

F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;

junk = M3d ;
junk(:,:,1:2) = 0 ;
Omega = junk(iwet) ;
% adjoint method.
kP    = d0(xhat_Cell.kP_T * par_CMIP.Tz + xhat_Cell.kdP) ;
Jex_P = d0(kP*G_Cell_CMIP(iwet))*(xhat_Cell.sigma*I+xhat_Cell.kappa_p*(1-xhat_Cell.sigma) * ...
                     inv(F_diag_p+xhat_Cell.kappa_p*I))*((T_diag+kP)\Omega);

P3d_Cell_CMIP = M3d+nan;
P3d_Cell_CMIP(iwet) = Jex_P;
clear kP Jex_P

% convert export from mmol P/m^3/s to mol P/m^2/yr;
TOPexp3d_Cell_CMIP = P3d_Cell_CMIP(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
TOPexp_Cell_CMIP = sum(TOPexp3d_Cell_CMIP,3); %sum(TOPexp3d,3,'omitnan');

% convert export from mmol P/m^3/s to mg P/m^2/day;
tem_Pexp = (P3d_Cell_CMIP(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31*spd).*dAt(:,:,1:nn);			% [mg P/day]
fprintf('Cell Model (CMIP2100): TOP export is %3.4f Pg P /yr   (Integrated to %4.1f m) \n',nansum(tem_Pexp(:))*365*1e-18, sum(grd.dzt(1:nn)));

%---------------- calculate carbon export -------------------------
	PFD_C = buildPFD(par_CMIP, 'POC') ;
	F_diag_p = inv(W)*PFD_C'*W   ;
	T_diag   = inv(W)*TRdiv'*W ;

	C2P_Cell_CMIP = model_Cell_CMIP.CellOut.C2P(iwet);
	Prod_Cell_CMIP  = G_Cell_CMIP(iwet).*C2P_Cell_CMIP ;
	% adjoint method.
	kC    = d0(xhat_Cell.kC_T * par_CMIP.Tz + xhat_Cell.kdC) ;
	Jex_C = d0(kC*Prod_Cell_CMIP)*(xhat_Cell.sigma*I+xhat_Cell.kappa_p*(1-xhat_Cell.sigma) * ...
	                     inv(F_diag_p+xhat_Cell.kappa_p*I))*((T_diag+kC)\Omega);

	C3d_Cell_CMIP = M3d+nan ;
	C3d_Cell_CMIP(iwet) = Jex_C ;
    clear kC Jex_C

    % convert export from mmol P/m^3/s to mol P/m^2/yr;
    TOCexp3d_Cell_CMIP = C3d_Cell_CMIP(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
    TOCexp_Cell_CMIP = sum(TOCexp3d_Cell_CMIP,3); %sum(TOPexp3d,3,'omitnan');

	% convert export from mmol C/m^3/s to mg C/m^2/day;
	tem_Cexp = (C3d_Cell_CMIP(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*12*spd).*dAt(:,:,1:nn);
	fprintf('Cell Model (CMIP2100): TOC export is %3.4f Pg C /yr   (Integrated to %4.1f m) \n',nansum(tem_Cexp(:))*365*1e-18, sum(grd.dzt(1:nn)));

% save
    EXPORT.TOPexp_Cell_CMIP = TOPexp_Cell_CMIP;
    EXPORT.TOCexp_Cell_CMIP = TOCexp_Cell_CMIP;


%% --------------------
save([figDir 'EXPORT_Cell_adjoint.mat'],'EXPORT')
