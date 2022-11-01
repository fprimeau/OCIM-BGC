%mkplot_C2Pzonalavg

% Plot Optimized zonal average uptake (solid) and export (dashed lines)
% C:P across all 3 models (P-only, T-only, and Cell model)
% for A) Modern, B) Temperature increase, C) Phosphate decrease

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


%% load model output
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/';
figDir = strcat(outputDir,'FIGS_CellGM15Tzv2/');


% load cell model output fields
fname = strcat(outputDir,'optC_Cell_CTL_He_PCCell_DOC0.25_DOP0.mat');
load(fname);
model_cell = data;
% load optimal parameter values
fxhat = strcat(outputDir,'optC_Cell_CTL_He_PCCell_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_cell = xhat.allparams;
clear data xhat

% load Tup 5C cell model output fields
fname_future = strcat(outputDir, 'futureproj/Tup5Cv2_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_cell_Tup5C = d.model;
clear d

% load P04obs down 20 percent cell model output fields
fname_future = strcat(outputDir, 'futureproj/Pdn20pct_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_cell_Pdn20 = d.model;
clear d

% load CMIP 2100 cell model output fields
fname_future = strcat(outputDir, 'futureproj/CMIP2100_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_cell_CMIP = d.model;
clear d

% % load PO4obs on;y down 20 percent cell model output fields
% fname_future = strcat(outputDir, 'futureproj/PO4dn20pct_Cell_CTL_He_PCCell.mat');
% d = load(fname_future);
% model_cell_PO4dn20 = d.model;
% clear d

% load GM15 model output
fname = strcat(outputDir,'optC_GM15_CTL_He_PC_DOC0.25_DOP0.mat');
load(fname);
model_GM15 = data;
fxhat = strcat(outputDir,'optC_GM15_CTL_He_PC_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_GM15 = xhat.allparams;
clear data xhat

% load Temp increase by 5C GM15 output fields
fname_future = strcat(outputDir, 'futureproj/Tup5Cv2_GM15_CTL_He_PC.mat');
d = load(fname_future);
model_GM15_Tup5C = d.model;
clear d

% load PO4obs decrease by 20% GM15 output fields
fname_future = strcat(outputDir, 'futureproj/Pdn20pct_GM15_CTL_He_PC.mat');
d = load(fname_future);
model_GM15_Pdn20 = d.model;
clear d

% load CMIP 2100 GM15 output fields
fname_future = strcat(outputDir, 'futureproj/CMIP2100_GM15_CTL_He_PC.mat');
d = load(fname_future);
model_GM15_CMIP = d.model;
clear d

% load Tonly model output
% load Tz model output
fname = strcat(outputDir,'optC_Tz_CTL_He_PC_DOC0.25_DOP0.mat');
load(fname);
model_Tz = data;
fxhat = strcat(outputDir,'optC_Tz_CTL_He_PC_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_Tz = xhat.allparams;
clear data xhat

% load Temp increase by 5C Tz output fields
fname_future = strcat(outputDir, 'futureproj/Tup5Cv2_Tz_CTL_He_PC.mat');
d = load(fname_future);
model_Tz_Tup5C = d.model;
clear d

% load PO4obs decrease by 20% Tz output fields
fname_future = strcat(outputDir, 'futureproj/Pdn20pct_Tz_CTL_He_PC.mat');
d = load(fname_future);
model_Tz_Pdn20 = d.model;
clear d

% load CMIP2100 Tz output fields
fname_future = strcat(outputDir, 'futureproj/CMIP2100_Tz_CTL_He_PC.mat');
d = load(fname_future);
model_Tz_CMIP = d.model;
clear d


%% Load Adjoint export
EXP = struct;

X = load([outputDir 'futureproj/v2_EXPORT_Cell_adjoint.mat'])
EXP.TOPexp_Cell_opt = X.EXPORT.TOPexp_Cell_opt;
EXP.TOCexp_Cell_opt = X.EXPORT.TOCexp_Cell_opt;
EXP.TOPexp_Cell_Tup5C = X.EXPORT.TOPexp_Cell_Tup5C;
EXP.TOCexp_Cell_Tup5C = X.EXPORT.TOCexp_Cell_Tup5C;
EXP.TOPexp_Cell_Pdn20 = X.EXPORT.TOPexp_Cell_Pdn20;
EXP.TOCexp_Cell_Pdn20 = X.EXPORT.TOCexp_Cell_Pdn20;
EXP.TOPexp_Cell_CMIP = X.EXPORT.TOPexp_Cell_CMIP;
EXP.TOCexp_Cell_CMIP = X.EXPORT.TOCexp_Cell_CMIP;

X = load([outputDir 'futureproj/v2_EXPORT_GM15_adjoint.mat'])
EXP.TOPexp_GM15_opt = X.EXPORT.TOPexp_GM15_opt;
EXP.TOCexp_GM15_opt = X.EXPORT.TOCexp_GM15_opt;
EXP.TOPexp_GM15_Tup5C = X.EXPORT.TOPexp_GM15_Tup5C;
EXP.TOCexp_GM15_Tup5C = X.EXPORT.TOCexp_GM15_Tup5C;
EXP.TOPexp_GM15_Pdn20 = X.EXPORT.TOPexp_GM15_Pdn20;
EXP.TOCexp_GM15_Pdn20 = X.EXPORT.TOCexp_GM15_Pdn20;
EXP.TOPexp_GM15_CMIP = X.EXPORT.TOPexp_GM15_CMIP;
EXP.TOCexp_GM15_CMIP = X.EXPORT.TOCexp_GM15_CMIP;

X = load([outputDir 'futureproj/v2_EXPORT_Tz_adjoint.mat'])
EXP.TOPexp_Tz_opt = X.EXPORT.TOPexp_Tz_opt;
EXP.TOCexp_Tz_opt = X.EXPORT.TOCexp_Tz_opt;
EXP.TOPexp_Tz_Tup5C = X.EXPORT.TOPexp_Tz_Tup5C;
EXP.TOCexp_Tz_Tup5C = X.EXPORT.TOCexp_Tz_Tup5C;
EXP.TOPexp_Tz_Pdn20 = X.EXPORT.TOPexp_Tz_Pdn20;
EXP.TOCexp_Tz_Pdn20 = X.EXPORT.TOCexp_Tz_Pdn20;
EXP.TOPexp_Tz_CMIP = X.EXPORT.TOPexp_Tz_CMIP;
EXP.TOCexp_Tz_CMIP = X.EXPORT.TOCexp_Tz_CMIP;

%% ------- compute zonal average C:P export fro adjoint------------------
% Cell model
C2Pexp_adj_latavg_Cell = mean(EXP.TOCexp_Cell_opt./EXP.TOPexp_Cell_opt,[2],'omitnan');

C2Pexp_adj_latavg_Cell_Tup5C = mean(EXP.TOCexp_Cell_Tup5C./EXP.TOPexp_Cell_Tup5C,[2],'omitnan');
C2Pexp_adj_latavg_Cell_Pdn20 = mean(EXP.TOCexp_Cell_Pdn20./EXP.TOPexp_Cell_Pdn20,[2],'omitnan');
C2Pexp_adj_latavg_Cell_CMIP = mean(EXP.TOCexp_Cell_CMIP./EXP.TOPexp_Cell_CMIP,[2],'omitnan');

% GM15 model
C2Pexp_adj_latavg_GM15 = mean(EXP.TOCexp_GM15_opt./EXP.TOPexp_GM15_opt,[2],'omitnan');
C2Pexp_adj_latavg_GM15_Tup5C = mean(EXP.TOCexp_GM15_Tup5C./EXP.TOPexp_GM15_Tup5C,[2],'omitnan');
C2Pexp_adj_latavg_GM15_Pdn20 = mean(EXP.TOCexp_GM15_Pdn20./EXP.TOPexp_GM15_Pdn20,[2],'omitnan');
C2Pexp_adj_latavg_GM15_CMIP = mean(EXP.TOCexp_GM15_CMIP./EXP.TOPexp_GM15_CMIP,[2],'omitnan');

% Tonly model
C2Pexp_adj_latavg_Tz = mean(EXP.TOCexp_Tz_opt./EXP.TOPexp_Tz_opt,[2],'omitnan');
C2Pexp_adj_latavg_Tz_Tup5C = mean(EXP.TOCexp_Tz_Tup5C./EXP.TOPexp_Tz_Tup5C,[2],'omitnan');
C2Pexp_adj_latavg_Tz_Pdn20 = mean(EXP.TOCexp_Tz_Pdn20./EXP.TOPexp_Tz_Pdn20,[2],'omitnan');
C2Pexp_adj_latavg_Tz_CMIP = mean(EXP.TOCexp_Tz_CMIP./EXP.TOPexp_Tz_CMIP,[2],'omitnan');

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
colors.mblue 		= [51/255 51/255 255/255];
colors.navy 		= [ 0 0 128/255];
colors.darkmagenta 	= [139/255 0 139/255];

postercolors.aqua2 = [88/255, 182/255, 192/255];
postercolors.aqua2dark25 = [58/255, 143/255, 152/255];
postercolors.aqua2dark50 = [38/255, 96/255, 102/255];
postercolors.aqua1 = [52/255, 148/255, 186/255];
postercolors.teal = [117/255, 189/255, 167/255];

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


load po4obs_91x180x24.mat % WOA PO4 observation

% load temperature
load tempobs_91x180x24.mat
tempobs(tempobs(:)<-2.0) = -2.0 ;
% normalize temperature
for ji = 1:24
	t2d = tempobs(:,:,ji);
	tempobs(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
end
vT = tempobs(iwet) ;
Tz01 = (vT - min(vT))./(max(vT) - min(vT)) ;
% Tz = zscore(vT)  ;
vT0min = min(vT);
vT0max = max(vT);
clear t2d vT tempobs


%load('/Users/megansullivan/Documents/UC Irvine/DATASETS/weilei_gp_DATA/M3d90x180x24v2.mat')
%iwet = find(M3d(:));
%nwet = length(iwet);
%addpath('/Users/megansullivan/Documents/UC Irvine/GitHub/WeiLei_code/my_func')

%% ------ compute export C:P
% weight
dVtwet = M3d*nan;
dVtwet(iwet) = dVt(iwet);
Wexp = dVtwet(:,:,1)/sum(dVtwet(:,:,1),'all','omitnan');

% Cell Model
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_cell = xhat_cell.kdC*model_cell.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_cell.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_cell = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('Cell Model integrated DOC below the Euphotic zone is %3.4f Pg C /yr \n',Sum_DOCremin_cell);
	DOCremin_cell = sum(DOCremin_cell,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_cell = xhat_cell.kdP*model_cell.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_cell.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_cell =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('Cell Model integrated DOP below the Euphotic zone is %3.4f Pg P /yr \n',Sum_DOPremin_cell);
    DOPremin_cell = sum(DOPremin_cell,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_cell = (DOCremin_cell)./(DOPremin_cell);
	C2Pexp_avg_cell = sum(C2Pexp_cell.*Wexp,'all','omitnan');
	fprintf('Cell Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_cell);

% GM15
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_GM15 = xhat_GM15.kdC*model_GM15.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_GM15.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_GM15 = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('GM15 Model integrated DOC below the Euphotic zone is %3.4f Pg C /yr \n',Sum_DOCremin_GM15);
	DOCremin_GM15 = sum(DOCremin_GM15,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_GM15 = xhat_GM15.kdP*model_GM15.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_GM15.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_GM15 =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('GM15 Model integrated DOP below the Euphotic zone is %3.4f Pg P /yr \n',Sum_DOPremin_GM15);
    DOPremin_GM15 = sum(DOPremin_GM15,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_GM15 = (DOCremin_GM15)./(DOPremin_GM15);
    C2Pexp_GM15(C2Pexp_GM15>1000) = NaN;

	C2Pexp_avg_GM15 = sum(C2Pexp_GM15.*Wexp,'all','omitnan');
	fprintf('GM15 Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_GM15);

% Tz
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_Tz = xhat_Tz.kdC*model_Tz.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_Tz.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_Tz = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('Tz Model integrated DOC below the Euphotic zone is %3.4f Pg C /yr \n',Sum_DOCremin_Tz);
	DOCremin_Tz = sum(DOCremin_Tz,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_Tz = xhat_Tz.kdP*model_Tz.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_Tz.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_Tz =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('Tz Model integrated DOP below the Euphotic zone is %3.4f Pg P /yr \n',Sum_DOPremin_Tz);
    DOPremin_Tz = sum(DOPremin_Tz,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_Tz = (DOCremin_Tz)./(DOPremin_Tz);
    C2Pexp_Tz(C2Pexp_Tz>1000) = NaN;

	C2Pexp_avg_Tz = sum(C2Pexp_Tz.*Wexp,'all','omitnan');
	fprintf('Tz Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_Tz);


%%%%%%

%% plot axes
lon_model = grd.xt;
lat_model = grd.yt;

%% C2P uptake lat average

C2P_cell     = model_cell.CellOut.C2P(:,:,1:2);
C2P_cell(C2P_cell==0) = NaN;

C2P_latavg_cell = mean(C2P_cell,[2 3],'omitnan');
C2P_latstd_cell = std(C2P_cell,0,[2 3],'omitnan');


C2P_GM15     = M3d + nan  ;
C2P_GM15(iwet)= 1./(xhat_GM15.cc*po4obs(iwet) + xhat_GM15.dd);
C2P_GM15(C2P_GM15==0) = NaN;
C2P_GM15 = C2P_GM15(:,:,1:2);

C2P_latavg_GM15 = mean(C2P_GM15,[2 3],'omitnan');
C2P_latstd_GM15 = std(C2P_GM15,0,[2 3],'omitnan');


% Tz model
C2P_Tz     = M3d + nan  ;
C2P_Tz(iwet)= 1./(xhat_Tz.ccT*Tz01 + xhat_Tz.ddT);
C2P_Tz(C2P_Tz==0) = NaN;
C2P_Tz = C2P_Tz(:,:,1:2);

C2P_latavg_Tz = mean(C2P_Tz,[2 3],'omitnan');
C2P_latstd_Tz = std(C2P_Tz,0,[2 3],'omitnan');

%% C2P export lat avg

% extract model C2P in Surface Atlantic
%C2Pexp_cell(C2Pexp_cell<0) = NaN;
%C2Pexp_GM15(C2Pexp_GM15<0) = NaN;

% M.ATL = MSKS.ATL(:,:,1);
% M.ATL(M.ATL==0)=NaN;
% C2Pexp_cell_ATL = C2Pexp_cell.*M.ATL;
% C2Pexp_GM15_ATL = C2Pexp_GM15.*M.ATL;

% zonal average C:P in the atlantic basin
C2Pexp_latavg_cell = mean(C2Pexp_cell,[2],'omitnan');
C2Pexp_std_cell = std(C2Pexp_cell,0,[2],'omitnan');

C2Pexp_latavg_GM15 = mean(C2Pexp_GM15,[2],'omitnan');
C2Pexp_std_GM15 = std(C2Pexp_GM15,0,[2],'omitnan');

C2Pexp_latavg_Tz = mean(C2Pexp_Tz,[2],'omitnan');
C2Pexp_std_Tz = std(C2Pexp_Tz,0,[2],'omitnan');


%% Future scenario: Temperature increased 5 C
C2P_cell_Tup5C     = model_cell_Tup5C.CellOut.C2P(:,:,1:2);
C2P_cell_Tup5C(C2P_cell_Tup5C==0) = NaN;

C2P_latavg_cell_Tup5C = mean(C2P_cell_Tup5C,[2 3],'omitnan');

% Cell Model
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_cell_Tup5C = xhat_cell.kdC*model_cell_Tup5C.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_cell_Tup5C.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_cell_Tup5C = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('T+5C : Cell Model integrated DOC below the Euphotic zone is %3.4f Pg C /yr \n',Sum_DOCremin_cell_Tup5C);
	DOCremin_cell_Tup5C = sum(DOCremin_cell_Tup5C,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_cell_Tup5C = xhat_cell.kdP*model_cell_Tup5C.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_cell_Tup5C.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_cell_Tup5C =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('T+5C : Cell Model integrated DOP below the Euphotic zone is %3.4f Pg P /yr \n',Sum_DOPremin_cell_Tup5C);
    DOPremin_cell_Tup5C = sum(DOPremin_cell_Tup5C,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_cell_Tup5C = (DOCremin_cell_Tup5C)./(DOPremin_cell_Tup5C);
	C2Pexp_avg_cell_Tup5C = sum(C2Pexp_cell_Tup5C.*Wexp,'all','omitnan');
	fprintf('T+5C : Cell Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_cell_Tup5C);

C2Pexp_latavg_cell_Tup5C = mean(C2Pexp_cell_Tup5C,[2],'omitnan');

%------ Tz model -----------
% ----overwrite temperature obs from SetUp -----
load tempobs_91x180x24.mat
tempobs(tempobs(:)<-2.0) = -2.0 ;
% Uniform  +5 degree C temperature increase
tempobs = tempobs+5;
Tempup5C    = tempobs ;

%-- normalize temperature --
for ji = 1:24
	t2d = Tempup5C(:,:,ji);
	Tempup5C(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
end
vT = Tempup5C(iwet) ;
Tz01_Tup5C = (vT - vT0min)./(vT0max - vT0min) ;
% Tz = zscore(vT)  ;
clear t2d vT tempobs

C2P_Tz_Tup5C     = M3d + nan  ;
C2P_Tz_Tup5C(iwet)     = 1./(xhat_Tz.ccT*Tz01_Tup5C + xhat_Tz.ddT);

C2P_Tz_Tup5C(C2P_Tz_Tup5C==0) = NaN;
C2P_Tz_Tup5C = C2P_Tz_Tup5C(:,:,1:2);

C2P_latavg_Tz_Tup5C = mean(C2P_Tz_Tup5C,[2 3],'omitnan');

% Tz Model
	%integrated DOC remineralization below the euphotic zone. (equal the TOC export)
	DOCremin_Tz_Tup5C = xhat_Tz.kdC*model_Tz_Tup5C.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_Tz_Tup5C.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_Tz_Tup5C = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('T+5C : Tz Model integrated DOC below the Euphotic zone is %3.4f Pg C /yr \n',Sum_DOCremin_Tz_Tup5C);
	DOCremin_Tz_Tup5C = sum(DOCremin_Tz_Tup5C,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

	%integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
	DOPremin_Tz_Tup5C = xhat_Tz.kdP*model_Tz_Tup5C.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
	tem_DOPremin = sum(DOPremin_Tz_Tup5C.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOPremin_Tz_Tup5C =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
	fprintf('T+5C : Tz Model integrated DOP below the Euphotic zone is %3.4f Pg P /yr \n',Sum_DOPremin_Tz_Tup5C);
	DOPremin_Tz_Tup5C = sum(DOPremin_Tz_Tup5C,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_Tz_Tup5C = (DOCremin_Tz_Tup5C)./(DOPremin_Tz_Tup5C);
	C2Pexp_avg_Tz_Tup5C = sum(C2Pexp_Tz_Tup5C.*Wexp,'all','omitnan');
	fprintf('T+5C : Tz Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_Tz_Tup5C);

C2Pexp_latavg_Tz_Tup5C = mean(C2Pexp_Tz_Tup5C,[2],'omitnan');


% ---------------------------------------------------
%% Future scenario: P down 20 percent

C2P_cell_Pdn20     = model_cell_Pdn20.CellOut.C2P(:,:,1:2);
C2P_cell_Pdn20(C2P_cell_Pdn20==0) = NaN;

C2P_latavg_cell_Pdn20 = mean(C2P_cell_Pdn20,[2 3],'omitnan');

% Cell Model
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_cell_Pdn20 = xhat_cell.kdC*model_cell_Pdn20.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_cell_Pdn20.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_cell_Tup5C = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('P-20%% : Cell Model integrated DOC below the Euphotic zone is %3.4f Pg C /yr \n',Sum_DOCremin_cell_Tup5C);
	DOCremin_cell_Pdn20 = sum(DOCremin_cell_Pdn20,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_cell_Pdn20 = xhat_cell.kdP*model_cell_Pdn20.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_cell_Pdn20.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_cell_Tup5C =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('P-20%% : Cell Model integrated DOP below the Euphotic zone is %3.4f Pg P /yr \n',Sum_DOPremin_cell_Tup5C);
    DOPremin_cell_Pdn20 = sum(DOPremin_cell_Pdn20,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_cell_Pdn20 = (DOCremin_cell_Pdn20)./(DOPremin_cell_Pdn20);
	C2Pexp_avg_cell_Pdn20 = sum(C2Pexp_cell_Pdn20.*Wexp,'all','omitnan');
	fprintf('P-20%% : Cell Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_cell_Pdn20);

C2Pexp_latavg_cell_Pdn20 = mean(C2Pexp_cell_Pdn20,[2],'omitnan');


%% GM15
C2P_GM15_Pdn20     = M3d + nan  ;
C2P_GM15_Pdn20(iwet)     = 1./(xhat_GM15.cc*(0.8*po4obs(iwet)) + xhat_GM15.dd);

C2P_GM15_Pdn20(C2P_GM15_Pdn20==0) = NaN;
C2P_GM15_Pdn20 = C2P_GM15_Pdn20(:,:,1:2);

C2P_latavg_GM15_Pdn20 = mean(C2P_GM15_Pdn20,[2 3],'omitnan');

% GM15 Model
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_GM15_Pdn20 = xhat_GM15.kdC*model_GM15_Pdn20.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_GM15_Pdn20.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_GM15_Pdn20 = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('P -20%% : GM15 Model integrated DOC below the Euphotic zone is %3.4f Pg C /yr \n',Sum_DOCremin_GM15_Pdn20);
	DOCremin_GM15_Pdn20 = sum(DOCremin_GM15_Pdn20,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_GM15_Pdn20 = xhat_GM15.kdP*model_GM15_Pdn20.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_GM15_Pdn20.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_GM15_Pdn20 =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('P -20%% : GM15  Model integrated DOP below the Euphotic zone is %3.4f Pg P /yr \n',Sum_DOPremin_GM15_Pdn20);
    DOPremin_GM15_Pdn20 = sum(DOPremin_GM15_Pdn20,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_GM15_Pdn20 = (DOCremin_GM15_Pdn20)./(DOPremin_GM15_Pdn20);
	C2Pexp_avg_GM15_Pdn20 = sum(C2Pexp_GM15_Pdn20.*Wexp,'all','omitnan');
	fprintf('P -20%% : GM15  Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_GM15_Pdn20);

C2Pexp_latavg_GM15_Pdn20 = mean(C2Pexp_GM15_Pdn20,[2],'omitnan');

%% --------- Tz model Pdn20pct -----------
%% Tz Pdn20pct
C2P_Tz_Pdn20     = C2P_Tz ;
C2P_latavg_Tz_Pdn20 = mean(C2P_Tz_Pdn20,[2 3],'omitnan');

% Tz Model
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_Tz_Pdn20 = xhat_Tz.kdC*model_Tz_Pdn20.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_Tz_Pdn20.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_Tz_Pdn20 = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('P -20%% : Tz Model integrated DOC below the Euphotic zone is %3.4f Pg C /yr \n',Sum_DOCremin_Tz_Pdn20);
	DOCremin_Tz_Pdn20 = sum(DOCremin_Tz_Pdn20,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_Tz_Pdn20 = xhat_Tz.kdP*model_Tz_Pdn20.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_Tz_Pdn20.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_Tz_Pdn20 =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('P -20%% : Tz  Model integrated DOP below the Euphotic zone is %3.4f Pg P /yr \n',Sum_DOPremin_Tz_Pdn20);
    DOPremin_Tz_Pdn20 = sum(DOPremin_Tz_Pdn20,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_Tz_Pdn20 = (DOCremin_Tz_Pdn20)./(DOPremin_Tz_Pdn20);
	C2Pexp_avg_Tz_Pdn20 = sum(C2Pexp_Tz_Pdn20.*Wexp,'all','omitnan');
	fprintf('P -20%% : Tz  Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_Tz_Pdn20);

C2Pexp_latavg_Tz_Pdn20 = mean(C2Pexp_Tz_Pdn20,[2],'omitnan');


%% -------------------------------------------------
%% Future scenario: CMIP projection to 2100

load('/DFS-L/DATA/primeau/meganrs/DATA/CMIP5/CMIP5mean_no3_thetao_91x180x24.mat')
% overwrite po4obs
po4proj_CMIP  = PO4_CMIP5mean  ;

% overwrite temperature obs from SetUp -----
T_CMIP5mean(T_CMIP5mean(:)<-2.0) = -2.0 ;
vT = T_CMIP5mean(iwet) ;
Tz01_CMIP = (vT - vT0min)./(vT0max - vT0min) ;


%% Cell model
C2P_cell_CMIP     = model_cell_CMIP.CellOut.C2P(:,:,1:2);
C2P_cell_CMIP(C2P_cell_CMIP==0) = NaN;
C2P_latavg_cell_CMIP = mean(C2P_cell_CMIP,[2 3],'omitnan');

%% GM15 CMIP2100
C2P_GM15_CMIP     = M3d + nan  ;
C2P_GM15_CMIP(iwet)     = 1./(xhat_GM15.cc*po4proj_CMIP(iwet) + xhat_GM15.dd);
C2P_GM15_CMIP = C2P_GM15_CMIP(:,:,1:2);
C2P_latavg_GM15_CMIP = mean(C2P_GM15_CMIP,[2 3],'omitnan');

%% Tz model CMIP2100
C2P_Tz_CMIP     = M3d + nan  ;
C2P_Tz_CMIP(iwet)     = 1./(xhat_Tz.ccT*Tz01_CMIP + xhat_Tz.ddT);
C2P_Tz_CMIP(C2P_Tz_CMIP==0) = NaN;
C2P_Tz_CMIP = C2P_Tz_CMIP(:,:,1:2);
C2P_latavg_Tz_CMIP = mean(C2P_Tz_CMIP,[2 3],'omitnan');


%--------------------- MAKE FIGURES -----------------------
%% -----Temp Plot

%% plot current C:P uptake and export
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')
h(1) = plot(grd.yt,C2P_latavg_cell,'-o','color',colors.darkgreen,'LineWidth',1.5,'MarkerSize',8);
hold on;
h(2) = plot(grd.yt,C2P_latavg_GM15,'-^','color',colors.lblue,'LineWidth',1.5,'MarkerSize',8); hold on
h(3) = plot(grd.yt,C2P_latavg_Tz,'-s','color',colors.darkmagenta,'LineWidth',1.5,'MarkerSize',8);hold on;

h(4) = plot(grd.yt,C2Pexp_adj_latavg_Cell,'--','color',colors.darkgreen,'LineWidth',2); hold on
h(5) = plot(grd.yt,C2Pexp_adj_latavg_GM15,'--','color',colors.lblue,'LineWidth',2); hold on
h(6) = plot(grd.yt,C2Pexp_adj_latavg_Tz,'--','color',colors.darkmagenta,'LineWidth',2); hold on

l=legend(h([1 2 3 4 5 6]),'Cell C:P_{uptake}','P-only C:P_{uptake}','T-only C:P_{uptake}','Cell C:P_{export}','P-only C:P_{export}','T-only C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P (mol C / mol P)')
title('optimal model $C:P_{uptake}$ and $C:P_{export}$')

figTitle = 'C2Platavg_upVexpadj_CellGM15Tz';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')
h(1) = plot(grd.yt,C2Pexp_adj_latavg_Cell,'-o','color',colors.darkgreen,'LineWidth',1.5,'MarkerSize',8);
hold on;
h(2) = plot(grd.yt,C2Pexp_adj_latavg_GM15,'-^','color',colors.lblue,'LineWidth',1.5,'MarkerSize',8); hold on
h(3) = plot(grd.yt,C2Pexp_adj_latavg_Tz,'-s','color',colors.darkmagenta,'LineWidth',1.5,'MarkerSize',8);hold on;

h(4) = plot(grd.yt,C2Pexp_latavg_cell,'--','color',colors.darkgreen,'LineWidth',2); hold on
h(5) = plot(grd.yt,C2Pexp_latavg_GM15,'--','color',colors.lblue,'LineWidth',2); hold on
h(6) = plot(grd.yt,C2Pexp_latavg_Tz,'--','color',colors.darkmagenta,'LineWidth',2); hold on

l=legend(h([1 2 3 4 5 6]),'Cell C:P_{export}','P-only C:P_{export}','T-only C:P_{export}','Cell C:P_{W.C.remin}','P-only C:P_{W.C.remin}','T-only C:P_{W.C.remin}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P (mol C / mol P)')
title('optimal model $C:P$ of export and of water column integrated remineralization')

figTitle = 'C2Platavg_expVexpadj_CellGM15Tz';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


% ------------- PLOTS -------------------------------------
%% C2P uptake for all
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,C2P_latavg_cell,'-o','color',colors.darkgreen,'LineWidth',1.5,'MarkerSize',8); hold on
h(2) = plot(grd.yt,C2P_latavg_cell_Tup5C,'--','color',colors.darkgreen,'LineWidth',2); hold on
h(3) = plot(grd.yt,C2P_latavg_cell_Pdn20,'-.','color',colors.darkgreen,'LineWidth',2,'MarkerSize',8); hold on

h(4) = plot(grd.yt,C2P_latavg_GM15,'-^','color',colors.lblue,'LineWidth',1.5,'MarkerSize',8); hold on
h(5) = plot(grd.yt,C2P_latavg_GM15,'--','color',colors.lblue,'LineWidth',2); hold on
h(6) = plot(grd.yt,C2P_latavg_GM15_Pdn20,'-.','color',colors.lblue,'LineWidth',2,'MarkerSize',8); hold on

h(7) = plot(grd.yt,C2P_latavg_Tz,'-^','color',colors.darkmagenta,'LineWidth',1.5,'MarkerSize',8); hold on
h(8) = plot(grd.yt,C2P_latavg_Tz_Tup5C,'--','color',colors.darkmagenta,'LineWidth',2); hold on
h(9) = plot(grd.yt,C2P_latavg_Tz_Pdn20,'-.','color',colors.darkmagenta,'LineWidth',2,'MarkerSize',8); hold on

l=legend(h([1 2 3 4 5 6 7 8 9]),'opt Cell','T+5 Cell','P-20% Cell','opt P-only','T+5 P-only','P-20% P-only','opt T-only','T+5 T-only','P-20% T-only','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('Cell model $C:P_{uptake}$ for Future projections')

figTitle = 'C2Platavg_uptake_all';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);



 %% T up 5 degrees for all models
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,C2P_latavg_cell_Tup5C,'-','color',colors.darkgreen,'LineWidth',2);
hold on;
h(2) = plot(grd.yt,C2Pexp_latavg_cell_Tup5C,'--','color',colors.darkgreen,'LineWidth',2); hold on
h(3) = plot(grd.yt,C2P_latavg_GM15,'-','color',colors.lblue,'LineWidth',2); hold on
h(4) = plot(grd.yt,C2Pexp_latavg_GM15,'--','color',colors.lblue,'LineWidth',2); hold on

h(5) = plot(grd.yt,C2P_latavg_Tz_Tup5C,'-','color',colors.darkmagenta,'LineWidth',2); hold on
h(6) = plot(grd.yt,C2Pexp_latavg_Tz_Tup5C,'--','color',colors.darkmagenta,'LineWidth',2); hold on

l=legend(h([1 2 3 4 5 6]),'Cell C:P_{uptake}','Cell C:P_{export}','P-only C:P_{uptake}','P-only C:P_{export}','T-only C:P_{uptake}','T-only C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('$C:P_{uptake}$ and $C:P_{export}$ for Future projection: Temperature +5C')

figTitle = 'C2Platavg_upVexp_Tup5C_all';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


%% DIP down 20% for all models

figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')
h(1) = plot(grd.yt,C2P_latavg_cell_Pdn20,'-','color',colors.darkgreen,'LineWidth',2);
hold on;
h(2) = plot(grd.yt,C2P_latavg_GM15_Pdn20,'-','color',colors.lblue,'LineWidth',2); hold on
h(3) = plot(grd.yt,C2P_latavg_Tz_Pdn20,'-','color',colors.darkmagenta,'LineWidth',2); hold on
h(4) = plot(grd.yt,C2Pexp_latavg_cell_Pdn20,'--','color',colors.darkgreen,'LineWidth',2); hold on
h(5) = plot(grd.yt,C2Pexp_latavg_GM15_Pdn20,'--','color',colors.lblue,'LineWidth',2); hold on
h(6) = plot(grd.yt,C2Pexp_latavg_Tz_Pdn20,'--','color',colors.darkmagenta,'LineWidth',2); hold on

l=legend(h([1 2 3 4 5 6]),'Cell C:P_{uptake}','P-only C:P_{uptake}','T-only C:P_{uptake}','Cell C:P_{export}','P-only C:P_{export}','T-only C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('$C:P_{uptake}$ and $C:P_{export}$ for Future projection: Phosphate-20\%')

figTitle = 'C2Platavg_upVexp_Pdn20_all';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


%% cell model only
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')
h(1) = plot(grd.yt,C2P_latavg_cell,'-o','color','k','LineWidth',2); hold on
h(2) = plot(grd.yt,C2Pexp_adj_latavg_Cell,'--','color',[0.3 0.3 0.3],'LineWidth',2); hold on
h(3) = plot(grd.yt,C2Pexp_latavg_cell,'-.','color','k','LineWidth',1); hold on

h(4) = plot(grd.yt,C2P_latavg_cell_Tup5C,'-r','LineWidth',1.5,'MarkerSize',8); hold on
h(5) = plot(grd.yt,C2Pexp_adj_latavg_Cell_Tup5C,'--','color',colors.tomato,'LineWidth',2); hold on
h(6) = plot(grd.yt,C2Pexp_latavg_cell_Tup5C,'-.','color',colors.tomato,'LineWidth',1); hold on

h(7) = plot(grd.yt,C2P_latavg_cell_Pdn20,'-b','LineWidth',2,'MarkerSize',8); hold on
h(8) = plot(grd.yt,C2Pexp_adj_latavg_Cell_Pdn20,'--','color',colors.lblue,'LineWidth',2); hold on
h(9) = plot(grd.yt,C2Pexp_latavg_cell_Pdn20,'-.','color',colors.lblue,'LineWidth',1); hold on


l=legend(h([1 2 3 4 5 6 7 8 9]),'C:P_{uptake}','C:P_{export}','C:P_{WCremin}','T+5 C:P_{uptake}','T+5 C:P_{export}','T+5 C:P_{WCremin}','P-20% C:P_{uptake}','P-20% C:P_{export}','P-20% C:P_{WCremin}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('Cell model $C:P_{uptake}$ and $C:P_{export}$ for Future projections')

figTitle = 'C2Platavg_upVexp_Cell';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);

%% cell model only
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')
h(1) = plot(grd.yt,C2P_latavg_cell,'-o','color','k','LineWidth',2); hold on
h(2) = plot(grd.yt,C2Pexp_adj_latavg_Cell,'--','color',[0.3 0.3 0.3],'LineWidth',1); hold on
%h(3) = plot(grd.yt,C2Pexp_latavg_cell,'-.','color','k','LineWidth',1); hold on

h(3) = plot(grd.yt,C2P_latavg_cell_Tup5C,'-r','LineWidth',2,'MarkerSize',8); hold on
h(4) = plot(grd.yt,C2Pexp_adj_latavg_Cell_Tup5C,'--','color',colors.tomato,'LineWidth',1); hold on
%h(6) = plot(grd.yt,C2Pexp_latavg_cell_Tup5C,'-.','color',colors.tomato,'LineWidth',1); hold on

h(5) = plot(grd.yt,C2P_latavg_cell_Pdn20,'-b','LineWidth',2,'MarkerSize',8); hold on
h(6) = plot(grd.yt,C2Pexp_adj_latavg_Cell_Pdn20,'--','color',colors.mblue,'LineWidth',1); hold on
%h(9) = plot(grd.yt,C2Pexp_latavg_cell_Pdn20,'-.','color',colors.lblue,'LineWidth',1); hold on

h(7) = plot(grd.yt,C2P_latavg_cell_CMIP,'-','color','m','LineWidth',2,'MarkerSize',8); hold on
h(8) = plot(grd.yt,C2Pexp_adj_latavg_Cell_CMIP,'--','color',colors.darkmagenta,'LineWidth',1); hold on

l=legend(h([1 2 3 4 5 6 7 8]),'C:P_{uptake}','C:P_{export}','T+5 C:P_{uptake}','T+5 C:P_{export}','P-20% C:P_{uptake}','P-20% C:P_{export}','CMIP2100 C:P_{uptake}','CMIP2100 C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('Cell model $C:P_{uptake}$ and $C:P_{export}$ for Future projections')

figTitle = 'C2Platavg_upVexpadj_Cell';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


%% GM15 model only
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,C2P_latavg_GM15,'-','color','k','LineWidth',2); hold on
h(2) = plot(grd.yt,C2Pexp_latavg_GM15,'--ok','color','k','LineWidth',2); hold on

h(3) = plot(grd.yt,C2P_latavg_GM15,'-r','LineWidth',1,'MarkerSize',8); hold on
hold on;
h(4) = plot(grd.yt,C2Pexp_latavg_GM15,'-.or','LineWidth',1.5,'MarkerSize',8); hold on

h(5) = plot(grd.yt,C2P_latavg_GM15_Pdn20,'-b','LineWidth',1.5,'MarkerSize',8); hold on
h(6) = plot(grd.yt,C2Pexp_latavg_GM15_Pdn20,'--ob','LineWidth',1.5,'MarkerSize',8); hold on

l=legend(h([1 2 3 4 5 6]),'C:P_{uptake}','C:P_{export}','T+5 C:P_{uptake}','T+5 C:P_{export}','P-20% C:P_{uptake}','P-20% C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('GM15 (phosphate-only) model $C:P_{uptake}$ and $C:P_{export}$ for Future projections')

figTitle = 'C2Platavg_upVexp_GM15';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);

figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')
h(1) = plot(grd.yt,C2P_latavg_GM15,'-^','color','k','LineWidth',2); hold on
h(2) = plot(grd.yt,C2Pexp_adj_latavg_GM15,'--','color',[0.3 0.3 0.3],'LineWidth',1); hold on
%h(3) = plot(grd.yt,C2Pexp_latavg_cell,'-.','color','k','LineWidth',1); hold on

h(3) = plot(grd.yt,C2P_latavg_GM15,'-r','LineWidth',2,'MarkerSize',8); hold on
h(4) = plot(grd.yt,C2Pexp_adj_latavg_GM15_Tup5C,'-.','color',colors.tomato,'LineWidth',1); hold on
%h(6) = plot(grd.yt,C2Pexp_latavg_cell_Tup5C,'-.','color',colors.tomato,'LineWidth',1); hold on

h(5) = plot(grd.yt,C2P_latavg_GM15_Pdn20,'-b','LineWidth',2,'MarkerSize',8); hold on
h(6) = plot(grd.yt,C2Pexp_adj_latavg_GM15_Pdn20,'--','color',colors.mblue,'LineWidth',1); hold on
%h(9) = plot(grd.yt,C2Pexp_latavg_cell_Pdn20,'-.','color',colors.lblue,'LineWidth',1); hold on

h(7) = plot(grd.yt,C2P_latavg_GM15_CMIP,'-','color','m','LineWidth',2,'MarkerSize',8); hold on
h(8) = plot(grd.yt,C2Pexp_adj_latavg_GM15_CMIP,'--','color',colors.darkmagenta,'LineWidth',1); hold on

l=legend(h([1 2 3 4 5 6 7 8]),'C:P_{uptake}','C:P_{export}','T+5 C:P_{uptake}','T+5 C:P_{export}','P-20% C:P_{uptake}','P-20% C:P_{export}','CMIP2100 C:P_{uptake}','CMIP2100 C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('[PO4]-only model $C:P_{uptake}$ and $C:P_{export}$ for Future projections')

figTitle = 'C2Platavg_upVexpadj_GM15';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


%% Temp model only
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,C2P_latavg_Tz,'-','color','k','LineWidth',2); hold on
h(2) = plot(grd.yt,C2Pexp_latavg_Tz,'--ok','color','k','LineWidth',2); hold on

h(3) = plot(grd.yt,C2P_latavg_Tz_Tup5C,'-r','LineWidth',1,'MarkerSize',8); hold on
hold on;
h(4) = plot(grd.yt,C2Pexp_latavg_Tz_Tup5C,'--or','LineWidth',1.5,'MarkerSize',8); hold on

h(5) = plot(grd.yt,C2P_latavg_Tz_Pdn20,'-b','LineWidth',1.5,'MarkerSize',8); hold on
h(6) = plot(grd.yt,C2Pexp_latavg_Tz_Pdn20,'-.ob','LineWidth',1.5,'MarkerSize',8); hold on


l=legend(h([1 2 3 4 5 6]),'C:P_{uptake}','C:P_{export}','T+5 C:P_{uptake}','T+5 C:P_{export}','P-20% C:P_{uptake}','P-20% C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('Tz (Temperature-only) model $C:P_{uptake}$ and $C:P_{export}$ for Future projections')

figTitle = 'C2Platavg_upVexp_Tz';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);

figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')
h(1) = plot(grd.yt,C2P_latavg_Tz,'-o','color','k','LineWidth',2); hold on
h(2) = plot(grd.yt,C2Pexp_adj_latavg_Tz,'--','color',[0.3 0.3 0.3],'LineWidth',1); hold on
%h(3) = plot(grd.yt,C2Pexp_latavg_Tz,'-.','color','k','LineWidth',1); hold on

h(3) = plot(grd.yt,C2P_latavg_Tz_Tup5C,'-r','LineWidth',2,'MarkerSize',8); hold on
h(4) = plot(grd.yt,C2Pexp_adj_latavg_Tz_Tup5C,'--','color',colors.tomato,'LineWidth',1); hold on
%h(6) = plot(grd.yt,C2Pexp_latavg_Tz_Tup5C,'-.','color',colors.tomato,'LineWidth',1); hold on

h(5) = plot(grd.yt,C2P_latavg_Tz_Pdn20,'-b','LineWidth',2,'MarkerSize',8); hold on
h(6) = plot(grd.yt,C2Pexp_adj_latavg_Tz_Pdn20,'-.','color',colors.mblue,'LineWidth',1); hold on
%h(9) = plot(grd.yt,C2Pexp_latavg_Tz_Pdn20,'-.','color',colors.lblue,'LineWidth',1); hold on

h(7) = plot(grd.yt,C2P_latavg_Tz_CMIP,'-','color','m','LineWidth',2,'MarkerSize',8); hold on
h(8) = plot(grd.yt,C2Pexp_adj_latavg_Tz_CMIP,'--','color',colors.darkmagenta,'LineWidth',1); hold on

l=legend(h([1 2 3 4 5 6 7 8]),'C:P_{uptake}','C:P_{export}','T+5 C:P_{uptake}','T+5 C:P_{export}','P-20% C:P_{uptake}','P-20% C:P_{export}','CMIP2100 C:P_{uptake}','CMIP2100 C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P')
title('T-only model $C:P_{uptake}$ and $C:P_{export}$ for Future projections')

figTitle = 'C2Platavg_upVexpadj_Tz';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


%%

%%
% figure;
% ind1 = ~isnan(C2P_latavg_cell);
% ind2 = ~isnan(C2P_latavg_GM15);
%
% h(1) = fill([lat(ind1),fliplr(lat(ind1))],[(C2P_latavg_now(ind1)-C2P_latstd_now(ind1))', fliplr((C2P_latavg_now(ind1)+C2P_latstd_now(ind1))')],colors.lblue,'LineStyle','none'); alpha(0.1); hold on
%
% h(2) = fill([lat(ind2),fliplr(lat(ind2))],[(C2P_latavg_future(ind2)-C2P_latstd_future(ind2))', fliplr((C2P_latavg_future(ind2)+C2P_latstd_future(ind2))')],'r','LineStyle','none'); alpha(0.1);
%
% h(3) = plot(lat,C2P_latavg_now,'-o','Color',colors.lblue); hold on
% %plot(lat,C2P_latmedian1,'-.b','linewidth',2);
% h(4) = plot(lat,C2P_latavg_future,'-ro');
%
% %h(5) = plot(lat,C2P_latavg_now1,'b--'); hold on
% %h(6) = plot(lat,C2P_latavg_now2,'b-.'); hold on
% %h(7) = plot(lat,C2P_latavg_future1,'--','Color',colors.maroon); hold on
% %h(8) = plot(lat,C2P_latavg_future2,'-.','Color',colors.maroon); hold on
%
% legend(h([3 4]),'Modern','T+5^oC');
% xlabel('Latitude')
% ylabel('C:P [molC:molP]')
% title('Zonal Average Cellular C:P')
% axis tight; grid on
% ylim([0 400])
% clear h;
