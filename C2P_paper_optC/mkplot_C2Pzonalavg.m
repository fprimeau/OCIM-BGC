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
figDir = strcat(outputDir,'FIGS_optC_Tz/');

% load GM15 model output
fname = strcat(outputDir,'optC_Tz_CTL_He_PC_DOC0.25_DOP0.mat');
load(fname);
model_Tz = data;
fxhat = strcat(outputDir,'optC_Tz_CTL_He_PC_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_Tz = xhat;
clear data xhat


% load cell model output fields
fname = strcat(outputDir,'optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0.mat');
load(fname);
model_cell = data;
% load optimal parameter values
fxhat = strcat(outputDir,'optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_cell = xhat;
clear data xhat


% load GM15 model output
fname = strcat(outputDir,'optC_GM15_CTL_He_PC_DOC0.25_DOP0.mat');
load(fname);
model_GM15 = data;
fxhat = strcat(outputDir,'optC_GM15_CTL_He_PC_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_GM15 = xhat;
clear data xhat


% load Tonly model output

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

postercolors.aqua2 = [88/255, 182/255, 192/255];
postercolors.aqua2dark25 = [58/255, 143/255, 152/255];
postercolors.aqua2dark50 = [38/255, 96/255, 102/255];
postercolors.aqua1 = [52/255, 148/255, 186/255];
postercolors.teal = [117/255, 189/255, 167/255];

%% -------------load data and set up parameters---------------------
% load grid
%addpath('/DFS-L/DATA/primeau/weilewang/DATA/')

load('/DFS-L/DATA/primeau/weilewang/DATA/M3d91x180x24.mat','MSKS')
OperName = sprintf('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2/OCIM2_CTL_He.mat');
    load(OperName,'output') ;
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
	par.TRdiv = -TR     ;
    %
iwet = find(M3d(:)) ;
dAt  = grd.DXT3d.*grd.DYT3d;
dVt  = dAt.*grd.DZT3d;

% -------- load observed PO4 and T fields -------------
load('/DFS-L/DATA/primeau/weilewang/DATA/po4obs_91x180x24.mat') % WOA PO4 observation

load('/DFS-L/DATA/primeau/weilewang/DATA/tempobs_91x180x24.mat')
tempobs(tempobs(:)<-2.0) = -2.0 ;
par.Temp    = tempobs ;
% --- normalize temperature -----
for ji = 1:24
    t2d = par.Temp(:,:,ji);
    par.Temp(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
end
vT = par.Temp(iwet) ;
Tz01 = (vT - min(vT))./(max(vT) - min(vT)) ;


%load('/Users/megansullivan/Documents/UC Irvine/DATASETS/weilei_gp_DATA/M3d90x180x24v2.mat')
%iwet = find(M3d(:));
%nwet = length(iwet);
%addpath('/Users/megansullivan/Documents/UC Irvine/GitHub/WeiLei_code/my_func')

%% plot axes
lon_model = grd.xt;
lat_model = grd.yt;


%% -------  C2P uptake lat average Tz
	C2P_Tz     = M3d + nan  ;
	C2P_Tz(iwet)= 1./(xhat_Tz.ccT*Tz01 + xhat_Tz.ddT);
	%C2P_Tz(C2P_Tz==0) = NaN;
	C2P_Tz = C2P_Tz(:,:,1:2);

	C2P_latavg_Tz = mean(C2P_Tz,[2 3],'omitnan');
	C2P_latstd_Tz = std(C2P_Tz,0,[2 3],'omitnan');

%% ------ compute export C:P for Tz model
% weight
dVtwet = M3d*nan;
dVtwet(iwet) = dVt(iwet);
Wexp = dVtwet(:,:,1)/sum(dVtwet(:,:,1),'all','omitnan');

% Tz
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_Tz = xhat_Tz.kdC*model_Tz.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_Tz.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_Tz = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('Tz Model integrated DOC below the Euphotic zone is %3.3e Pg C /yr \n',Sum_DOCremin_Tz);
	DOCremin_Tz = sum(DOCremin_Tz,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_Tz = xhat_Tz.allparams.kdP*model_Tz.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_Tz.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_Tz =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('Tz Model integrated DOP below the Euphotic zone is %3.3e Pg P /yr \n',Sum_DOPremin_Tz);
    DOPremin_Tz = sum(DOPremin_Tz,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_Tz = (DOCremin_Tz)./(DOPremin_Tz);
	ibadC2Pexp = find(C2Pexp_Tz>1000);
	fprintf('Tz Model: removed %i extreme C2Pexp values \n',length(ibadC2Pexp))
    C2Pexp_Tz(C2Pexp_Tz>1000) = NaN;

	C2Pexp_avg_Tz = sum(C2Pexp_Tz.*Wexp,'all','omitnan');
	fprintf('Tz Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_Tz);

	% ------ zonal average of C:P export ---------------
	C2Pexp_latavg_Tz = mean(C2Pexp_Tz,[2],'omitnan');
	C2Pexp_std_Tz = std(C2Pexp_Tz,0,[2],'omitnan');


%% -------- plot current C:P uptake and export --------
	figure;
	clear h
	set(gcf,'position',[61 614 1331 403],'units','pixels')
	h(1) = plot(grd.yt,C2P_latavg_Tz,'-','color',colors.darkmagenta,'LineWidth',2);
	hold on;
	h(2) = plot(grd.yt,C2Pexp_latavg_Tz,'--','color',colors.darkmagenta,'LineWidth',2); hold on

	l=legend(h([1 2]),'C:P_{uptake}','C:P_{export}','Location','NorthEast');
	grid on; axis tight
	xlabel('Latitude')
	ylabel('C:P')
	title('optimal T-only model $C:P_{uptake}$ and $C:P_{export}$')

	figTitle = 'C2Platavg_upVexp_Tz';
	exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


%% ----------------------------------------------------------------------------
%% ------------- Compare to other models --------------------------------------
%compute export C:P
% Cell Model
    %integrated DOC remineralization below the euphotic zone. (equal the TOC export)
    DOCremin_cell = xhat_cell.kdC*model_cell.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCremin = sum(DOCremin_cell.*dAt(:,:,3:end),3,'omitnan');
	Sum_DOCremin_cell = sum(tem_DOCremin(:),'omitnan')*365*1e-18;
	fprintf('Cell Model integrated DOC below the Euphotic zone is %3.3e Pg C /yr \n',Sum_DOCremin_cell);
	DOCremin_cell = sum(DOCremin_cell,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_cell = xhat_cell.allparams.kdP*model_cell.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_cell.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_cell =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('Cell Model integrated DOP below the Euphotic zone is %3.3e Pg P /yr \n',Sum_DOPremin_cell);
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
	fprintf('GM15 Model integrated DOC below the Euphotic zone is %3.3e Pg C /yr \n',Sum_DOCremin_GM15);
	DOCremin_GM15 = sum(DOCremin_GM15,3,'omitnan')/12/1000*365; % [mol C/m^2 /yr]

    %integrated DOP remineralization below the euphotic zone. (equal the TOP export calculated by the adjoint method)
    DOPremin_GM15 = xhat_GM15.allparams.kdP*model_GM15.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
    tem_DOPremin = sum(DOPremin_GM15.*dAt(:,:,3:end),3,'omitnan');
    Sum_DOPremin_GM15 =sum(tem_DOPremin(:),'omitnan')*365*1e-18;
    fprintf('GM15 Model integrated DOP below the Euphotic zone is %3.3e Pg P /yr \n',Sum_DOPremin_GM15);
    DOPremin_GM15 = sum(DOPremin_GM15,3,'omitnan')/31/1000*365; % [mol P/m^2 /yr]

	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexp_GM15 = (DOCremin_GM15)./(DOPremin_GM15);
	ibadC2Pexp = find(C2Pexp_GM15>1000);
	fprintf('GM15 Model: removed %i extreme C2Pexp values \n',length(ibadC2Pexp))
	C2Pexp_GM15(C2Pexp_GM15>1000) = NaN;

	C2Pexp_avg_GM15 = sum(C2Pexp_GM15.*Wexp,'all','omitnan');
	fprintf('GM15 Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexp_avg_GM15);
%


%% C2P uptake lat average
% cell model
C2P_cell     = model_cell.CellOut.C2P(:,:,1:2);
C2P_cell(C2P_cell==0) = NaN;

C2P_latavg_cell = mean(C2P_cell,[2 3],'omitnan');
C2P_latstd_cell = std(C2P_cell,0,[2 3],'omitnan');

% GM15
C2P_GM15     = M3d + nan  ;
C2P_GM15(iwet)= 1./(xhat_GM15.cc*po4obs(iwet) + xhat_GM15.dd);
%C2P_GM15(C2P_GM15==0) = NaN;
C2P_GM15 = C2P_GM15(:,:,1:2);

C2P_latavg_GM15 = mean(C2P_GM15,[2 3],'omitnan');
C2P_latstd_GM15 = std(C2P_GM15,0,[2 3],'omitnan');


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


%% plot current C:P uptake and export
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')
h(1) = plot(grd.yt,C2P_latavg_cell,'-o','color',colors.darkgreen,'LineWidth',1.5,'MarkerSize',8);
hold on;
h(2) = plot(grd.yt,C2P_latavg_GM15,'-^','color',colors.lblue,'LineWidth',1.5,'MarkerSize',8); hold on
h(3) = plot(grd.yt,C2P_latavg_Tz,'-s','color',colors.darkmagenta,'LineWidth',1.5,'MarkerSize',8);hold on;

h(4) = plot(grd.yt,C2Pexp_latavg_cell,'--','color',colors.darkgreen,'LineWidth',2); hold on
h(5) = plot(grd.yt,C2Pexp_latavg_GM15,'--','color',colors.lblue,'LineWidth',2); hold on
h(6) = plot(grd.yt,C2Pexp_latavg_Tz,'--','color',colors.darkmagenta,'LineWidth',2); hold on

l=legend(h([1 2 3 4 5 6]),'Cell C:P_{uptake}','P-only C:P_{uptake}','T-only C:P_{uptake}','Cell C:P_{export}','P-only C:P_{export}','T-only C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P (mol C / mol P)')
title('optimal model $C:P_{uptake}$ and $C:P_{export}$')

figTitle = 'C2Platavg_upVexp_CellGM15Tz';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);
