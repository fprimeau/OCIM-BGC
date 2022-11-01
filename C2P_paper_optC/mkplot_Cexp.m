%% plot Carbon export from adjoint method
clear all; close all
%% load in export fields from all runs
spd = 60*60*24;
spa = spd*365;
if ismac
    %load('/Users/megansullivan/Documents/UC Irvine/GitHub/TraitModel_output/C2P_paper/EXPORT_adjointv2.mat')
    outputDir = '/Users/megansullivan/Documents/UC Irvine/GitHub/TraitModel_output/C2P_paper/';

    load('/Users/megansullivan/Documents/UC Irvine/DATASETS/weilei_gp_DATA/M3d91x180x24.mat');
    load('/Users/megansullivan/Documents/UC Irvine/DATASETS/weilei_gp_DATA/OCIM2_Grid.mat');

else
    %load('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper/EXPORT_adjointv2.mat')
    outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/';
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


    load('/DFS-L/DATA/primeau/weilewang/DATA/M3d91x180x24.mat','MSKS')
    OperName = sprintf('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2/OCIM2_CTL_He.mat');
    load(OperName,'output') ;
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
    %

end
% ---------------- choose which Phosphate dn 20pct version -----
affectonlyC2P = true; % on if PO4 input only affects C2P and not production
figDir = strcat(outputDir,'FIGS_CellGM15Tzv2/');
%figDirCell = strcat(outputDir,'FIGS_optC_Cell/');
%figDirGM15 = strcat(outputDir,'FIGS_optC_GM15/');

iwet = find(M3d(:)) ;
dAt  = grd.DXT3d.*grd.DYT3d;
dVt  = dAt.*grd.DZT3d;


%% ------ set up figures ---------

set(groot,'defaultAxesFontName','Times',...
    'defaultAxesFontSize',14,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','on',...
    'defaultAxesYMinorTick','on');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Times',...
    'defaultTextInterpreter','latex');

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
colors.darkgrey 	= [0.3 0.3 0.3];

%% ------- compute zonal average C export ------------------
% Cell model
Cexp_latavg_Cell = mean(EXP.TOCexp_Cell_opt,[2],'omitnan');
Cexp_std_Cell = std(EXP.TOCexp_Cell_opt,0,[2],'omitnan');

Cexp_latavg_Cell_Tup5C = mean(EXP.TOCexp_Cell_Tup5C,[2],'omitnan');
Cexp_latavg_Cell_Pdn20 = mean(EXP.TOCexp_Cell_Pdn20,[2],'omitnan');
Cexp_latavg_Cell_CMIP = mean(EXP.TOCexp_Cell_CMIP,[2],'omitnan');

% GM15 model
Cexp_latavg_GM15 = mean(EXP.TOCexp_GM15_opt,[2],'omitnan');
Cexp_latavg_GM15_Tup5C = mean(EXP.TOCexp_GM15_Tup5C,[2],'omitnan');
Cexp_latavg_GM15_Pdn20 = mean(EXP.TOCexp_GM15_Pdn20,[2],'omitnan');
Cexp_latavg_GM15_CMIP = mean(EXP.TOCexp_GM15_CMIP,[2],'omitnan');

% Tonly model
Cexp_latavg_Tz = mean(EXP.TOCexp_Tz_opt,[2],'omitnan');
Cexp_latavg_Tz_Tup5C = mean(EXP.TOCexp_Tz_Tup5C,[2],'omitnan');
Cexp_latavg_Tz_Pdn20 = mean(EXP.TOCexp_Tz_Pdn20,[2],'omitnan');
Cexp_latavg_Tz_CMIP = mean(EXP.TOCexp_Tz_CMIP,[2],'omitnan');


%% Plot
% ----- plot cell model zonal average C export ---------
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,Cexp_latavg_Cell,'-ok','color','k','LineWidth',2,'MarkerSize',8); hold on
h(2) = plot(grd.yt,Cexp_latavg_Cell_Tup5C,'--or','LineWidth',1.5,'MarkerSize',8); hold on
h(3) = plot(grd.yt,Cexp_latavg_Cell_Pdn20,'--ob','LineWidth',1.5,'MarkerSize',8); hold on
h(4) = plot(grd.yt,Cexp_latavg_Cell_CMIP,'--om','LineWidth',1.5,'MarkerSize',8); hold on

l=legend(h([1 2 3 4]),'opt C_{export}','T+5 C_{export}','P-20% C_{export}','CMIP2100 C_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C export [$mol C/m^2/yr$]')
title('Cell model zonal average $C_{export}$ for optimal fit and future projections')

figTitle = 'Cexp_adj_future_latavg_Cell';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);

%% ------- plot GM15 model zonal average C export ---------
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,Cexp_latavg_GM15,'-^k','color','k','LineWidth',2,'MarkerSize',8); hold on
h(2) = plot(grd.yt,Cexp_latavg_GM15_Tup5C,'--^r','LineWidth',1.5,'MarkerSize',8); hold on
h(3) = plot(grd.yt,Cexp_latavg_GM15_Pdn20,'--^b','LineWidth',1.5,'MarkerSize',8); hold on
h(4) = plot(grd.yt,Cexp_latavg_GM15_CMIP,'--^m','LineWidth',1.5,'MarkerSize',8); hold on

l=legend(h([1 2 3 4]),'opt C_{export}','T+5 C_{export}','P-20% C_{export}','CMIP2100 C_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C export [$mol C/m^2/yr$]')
title('P-only model zonal average $C_{export}$ for optimal fit and future projections')

figTitle = 'Cexp_adj_future_latavg_GM15';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);

%% ----- plot Tz model zonal average C export ---------
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,Cexp_latavg_Tz,'-sk','color','k','LineWidth',2,'MarkerSize',8); hold on
h(2) = plot(grd.yt,Cexp_latavg_Tz_Tup5C,'--sr','LineWidth',1.5,'MarkerSize',8); hold on
h(3) = plot(grd.yt,Cexp_latavg_Tz_Pdn20,'--sb','LineWidth',1.5,'MarkerSize',8); hold on
h(4) = plot(grd.yt,Cexp_latavg_Tz_CMIP,'--sb','LineWidth',1.5,'MarkerSize',8); hold on

l=legend(h([1 2 3 4]),'opt C_{export}','T+5 C_{export}','P-20% C_{export}','CMIP2100 C_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C export [$mol C/m^2/yr$]')
title('T-only model zonal average $C_{export}$ for optimal fit and future projections')

figTitle = 'Cexp_adj_future_latavg_Tz';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


%% -------- plot all models opt zonal average C export ----------
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,Cexp_latavg_Cell,'-o','color',colors.darkgreen,'LineWidth',1.5,'MarkerSize',8); hold on
h(2) = plot(grd.yt,Cexp_latavg_Cell_Tup5C,'-.','color',colors.darkgreen,'LineWidth',2); hold on
h(3) = plot(grd.yt,Cexp_latavg_Cell_Pdn20,'-.','color',colors.darkgreen,'LineWidth',2); hold on
h(4) = plot(grd.yt,Cexp_latavg_Cell_CMIP,'-^','color',colors.darkgreen,'LineWidth',1,'MarkerSize',8); hold on

h(5) = plot(grd.yt,Cexp_latavg_GM15,'-o','color',colors.lblue,'LineWidth',1.5,'MarkerSize',8); hold on
h(6) = plot(grd.yt,Cexp_latavg_GM15_Tup5C,'-.','color',colors.lblue,'LineWidth',2); hold on
h(7) = plot(grd.yt,Cexp_latavg_GM15_Pdn20,'-','color',colors.lblue,'LineWidth',2); hold on
h(8) = plot(grd.yt,Cexp_latavg_GM15_CMIP,'-^','color',colors.lblue,'LineWidth',1,'MarkerSize',8); hold on

h(9) = plot(grd.yt,Cexp_latavg_Tz,'-o','color',colors.darkmagenta,'LineWidth',1.5,'MarkerSize',8); hold on
h(10) = plot(grd.yt,Cexp_latavg_Tz_Tup5C,'-.','color',colors.darkmagenta,'LineWidth',2); hold on
h(11) = plot(grd.yt,Cexp_latavg_Tz_Pdn20,'-','color',colors.darkmagenta,'LineWidth',2); hold on
h(12) = plot(grd.yt,Cexp_latavg_Tz_CMIP,'-^','color',colors.darkmagenta,'LineWidth',1,'MarkerSize',8); hold on

l=legend(h([1 2 3 4 5 6 7 8 9 10 11 12]),'opt Cell','T+5C Cell','P-20% Cell','CMIP2100 Cell','opt P-only','T+5C P-only','P-20% P-only','CMIP2100 P-only','opt T-only','T+5C T-only','P-20% T-only','CMIP2100 T-only','Location','EastOutside');

grid on; axis tight
xlabel('Latitude')
ylabel('C export [$mol C/m^2/yr$]')
title('Zonal Average $C_{export}$ for all models')

figTitle = 'Cexp_adj_future_latavg_all';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);

%{

%% ----------- C export zonal sum -------------------
Cexp_latsum_Cell = sum(EXPORT.TOCexp_Cell_opt.*dAt(:,:,1)*12,[2],'omitnan')*1e-15;

Cexp_latsum_GM15 = sum(EXPORT.TOCexp_GM15_opt.*dAt(:,:,1)*12,[2],'omitnan')*1e-15;

Cexp_latsum_Tz = sum(EXPORT.TOCexp_Tz_opt.*dAt(:,:,1)*12,[2],'omitnan')*1e-15;

%% plot Cell and GM15 total C export
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,Cexp_latsum_Cell,'-o','color',colors.darkgreen,'LineWidth',1.5); hold on
h(2) = plot(grd.yt,Cexp_latsum_GM15,'-^','color',colors.lblue,'LineWidth',1.5); hold on
h(3) = plot(grd.yt,Cexp_latsum_Tz,'-s','color',colors.darkmagenta,'LineWidth',1.5); hold on

l=legend(h([1 2 3]),'Cell opt C_{export}','GM15 opt C_{export}','Tonly opt C_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C export [$Pg C/yr$]')
title('Zonal Total $C_{export}$ for optimal fit')

figTitle = 'Cexp_latsum_opt_all';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);

%%


%% --------- compute C2P lat averages ----------
% zonal average C export
C2Pexp_latavg_Cell = mean(EXPORT.TOCexp_Cell_opt./EXPORT.TOPexp_Cell_opt,[2],'omitnan');
C2Pexp_latavg_GM15 = mean(EXPORT.TOCexp_GM15_opt./EXPORT.TOPexp_GM15_opt,[2],'omitnan');
C2Pexp_latavg_Tz = mean(EXPORT.TOCexp_Tz_opt./EXPORT.TOPexp_Tz_opt,[2],'omitnan');


%% plot C2P lat average
figure;
clear h
set(gcf,'position',[61 614 1331 403],'units','pixels')

h(1) = plot(grd.yt,C2Pexp_latavg_Cell,'-o','color',colors.darkgreen,'LineWidth',1.5); hold on
h(2) = plot(grd.yt,C2Pexp_latavg_GM15,'-^','color',colors.lblue,'LineWidth',1.5); hold on
h(3) = plot(grd.yt,C2Pexp_latavg_Tz,'-s','color',colors.darkmagenta,'LineWidth',1.5); hold on

l=legend(h([1 2 3]),'Cell opt C:P_{export}','GM15 opt C:P_{export}','Tonly opt C:P_{export}','Location','EastOutside');
grid on; axis tight
xlabel('Latitude')
ylabel('C:P export [$mol C/ molP$]')
title('Zonal Average $C:P_{export}$ for optimal fit')

figTitle = 'C2Pexp_adjoint_latavg_opt_all';
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);

%}


% %% plot Cell model zonal average C export
% figure;
% clear h
% set(gcf,'position',[61 614 1331 403],'units','pixels')
%
% h(1) = plot(grd.yt,Cexp_latavg_Cell,'--ok','color','k','LineWidth',2); hold on
% h(2) = plot(grd.yt,Cexp_latavg_Cell_Tup5C,'--or','LineWidth',1.5,'MarkerSize',8); hold on
%
% h(3) = plot(grd.yt,Cexp_latavg_Cell_Pdn20,'--ob','LineWidth',1.5,'MarkerSize',8); hold on
%
%
% l=legend(h([1 2 3]),'opt C_{export}','T+5 C_{export}','P-20% C_{export}','Location','EastOutside');
% grid on; axis tight
% xlabel('Latitude')
% ylabel('C export [$mol C/m^2/yr$]')
% title('Cell model $C_{export}$ for optimal fit and future projections')
%
% figTitle = 'Cexp_latavg_Cell';
% exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);


% %% plot Cell and GM15 zonal average C export
% figure;
% clear h
% set(gcf,'position',[61 614 1331 403],'units','pixels')
%
% h(1) = plot(grd.yt,Cexp_latavg_Cell,'-ok','color','k','LineWidth',2); hold on
% h(2) = plot(grd.yt,Cexp_latavg_Cell_Tup5C,'-o','color',colors.tomato,'LineWidth',1.5); hold on
%
% h(3) = plot(grd.yt,Cexp_latavg_Cell_Pdn20,'-ob','LineWidth',1.5); hold on
% h(4) = plot(grd.yt,Cexp_latavg_GM15,'--^k','color',colors.darkgrey,'LineWidth',2); hold on
%
% h(5) = plot(grd.yt,Cexp_latavg_GM15_Pdn20,'--^','color',colors.lblue,'LineWidth',1.5); hold on
%
% l=legend(h([1 2 3 4 5]),'Cell opt C_{export}','Cell T+5 C_{export}','Cell P-20% C_{export}','GM15 opt C_{export}','GM15 P-20% C_{export}','Location','EastOutside');
% grid on; axis tight
% xlabel('Latitude')
% ylabel('C export [$mol C/m^2/yr$]')
% title('Zonal Average $C_{export}$ for optimal fit and future projections')
%
% figTitle = 'Cexp_latavg_CellvGM15';
% %exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);
