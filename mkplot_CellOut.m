clc; clear all; close all
on   = true  ;
off  = false ;

%ver = datestr(now,'mmmdd');
RunVer = 'v6b_duplicate/Tv4_PCCellv6b_DOC0.25_DOP0';

%model output directory
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/';
figDir = strcat(outputDir,'FIGS_PCCellv6b_DOC0.25_DOP0/');
outPath = figDir;

% load model output fields
fname = strcat(outputDir, RunVer, '.mat');
load(fname);
model = data;

% load optimal parameter values
fxhat = strcat(outputDir, RunVer,'_xhat.mat');
load(fxhat);

GridVer  = 90  ;
operator = 'A' ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

%-------------load data and set up parameters---------------------
SetUp ;
xhat

%-------------set up figure properties------------------
% SET DEFAULT FIGURE SIZE
%set(groot, 'defaultFigureUnits','normalized',...
%    'defaultFigurePosition', [0 0.2 0.6 0.4],...
%    'defaultFigurePaperPositionMode','auto');
% SET DEFAULT AXES PROPERTIES
%   FontName
%   FontSize
%   Grid = on
set(groot,'defaultAxesFontName','Times',...
    'defaultAxesFontSize',14,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','on',...
    'defaultAxesYMinorTick','on');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Times',...
    'defaultTextInterpreter','latex');

%Colors
aqua = [0.2 0.8 0.8];
teal = [0 128/255 128/255];
darkgreen = [0 100/255 0];

lblue = [0 191/255 255/255];
navy = [ 0 0 128/255];

%------------------------------------------------------
%% rename variables
C2P     = model.CellOut.C2P;
N2P     = model.CellOut.N2P;
C2N     = model.CellOut.C2N;
radius  = model.CellOut.r;
LimType = model.CellOut.LimType;
mu      = model.CellOut.mu;

% set land points to NaN
C2P(find(C2P==0)) = NaN;
C2N(find(C2N==0)) = NaN;
N2P(find(N2P==0)) = NaN;
radius(find(radius==0)) = NaN;

%% plot axes
if GridVer == 90;
    lat = [-89:2:89];
    lon = [1:2:360];
end

if GridVer == 91;
    lat = [-90:2:90];
    lon = [1:2:360];
end

%% annotation for parameter values
dim = [0.0 0.6 0.1 0.2];
dim = [0.1199 0.695 0.1 0.2]
%parstr = 'Cell Model Parameters';
parstr = {'Cell Model Parameters'}
kk =1;
if isfield(xhat,'Q10Photo')
	kk=kk+1;
    parstr{kk,1} = ['Q10Photo=' num2str(xhat.Q10Photo)];
end
if isfield(xhat,'fStorage')
	kk=kk+1;
    parstr{kk,1} = ['fStorage=' num2str(xhat.fStorage)] ;
end
if isfield(xhat,'fRibE')
	kk=kk+1;
    parstr{kk,1} = ['fRibE=' num2str(xhat.fRibE)] ;
end
if isfield(xhat,'kST0')
	kk=kk+1;
    parstr{kk,1} = ['kST0=' num2str(xhat.kST0)] ;
end
if isfield(xhat,'PLip_PCutoff')
	kk=kk+1;
    parstr{kk,1} = ['PLip PCutoff=' num2str(xhat.PLip_PCutoff)] ;
end
if isfield(xhat,'PLip_scale')
	kk=kk+1;
    parstr{kk,1} = ['PLip scale=' num2str(xhat.PLip_scale)] ;
end
if isfield(xhat,'PStor_rCutoff')
	kk=kk+1;
    parstr{kk,1} = ['PStor rCutoff=' num2str(xhat.PStor_rCutoff)] ;
end
if isfield(xhat,'PStor_scale')
	kk=kk+1;
    parstr{kk,1} = ['PStor scale=' num2str(xhat.PStor_scale)];
end
if isfield(xhat,'alphaS')
	kk=kk+1;
    parstr{kk,1} = ['alphaS=' num2str(xhat.alphaS)]
end

%% model POC/POP
POC_POP = M3d*NaN;
indx1=find(model.POC(iwet)>0 & model.POP(iwet)>0);
POC_POP(iwet) = -1;
POC_POP(iwet(indx1)) = model.POC(iwet(indx1))./model.POP(iwet(indx1));
for ii = 1:4
	figure
	contourf(lon,lat,POC_POP(:,:,ii)); hold on
	cb=colorbar;
	colormap(flipud(summer))
	title(['Cell Model POC:POP concentration: Z level ' num2str(ii)],'Fontsize',18);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'POC:POP [ ]');
	figTitle = ['POC2POP_Z' num2str(ii)];
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
end
%why are there negative values?
%need to make dry poits nan

%% POC:POP & DOC:DOP scatter plot
%% plot Pstorage & PLip?

%% C2P surface plot
%Zlevs = [100:20:300];
figure;
contourf(lon,lat,C2P(:,:,1)); hold on
cb=colorbar;
colormap(flipud(summer))
title('Cell Model C:P Uptake Ratio: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'C:P [molC/molP]');

%keyboard;
%ax = gca;
%ax.Position(1) = 0.2;
%cb.Position(4) = 0.6;
%dim = [0.8490 0.95 0.15 0.25]

annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'C2Psurface';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% C2P lower EZ
%Zlevs =
figure; hold on;
%contourf(lon,lat,C2P(:,:,2)); hold on
imAlpha = ones(size(C2P(:,:,2)));
imAlpha(isnan(C2P(:,:,2))) =0;
imagesc(lon,lat,C2P(:,:,2),'AlphaData',imAlpha)

c=colorbar;
colormap(flipud(summer));
title('Cell Model C:P Uptake Ratio: Lower EZ','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'C:P [molC/molP]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'C2P_Z2';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% C2P as a function of latitude
C2P_latavg1 = mean(C2P(:,:,1),2,'omitnan');
C2P_latavg2 = mean(C2P(:,:,2),2,'omitnan');
C2P_latavg = mean(C2P(:,:,1:2),[2 3],'omitnan');
C2P_latmedian1 = median(C2P(:,:,1),2,'omitnan');
C2P_latmedian2 = median(C2P(:,:,2),2,'omitnan');
C2P_latstd1 = std(C2P(:,:,1),0,2,'omitnan');
C2P_latstd2 = std(C2P(:,:,2),0,2,'omitnan');
ind1 = ~isnan(C2P_latavg1);
ind2 = ~isnan(C2P_latavg2);

figure;
% plot +/-1 standard deviation
% boundedline(lat,C2P_latavg1,C2P_latstd1,'-b*',lat,C2P_latavg2,C2P_latstd2,'-m*','alpha','nan','gap')
h(1) = fill([lat(ind1),fliplr(lat(ind1))],[(C2P_latavg1(ind1)-C2P_latstd1(ind1))', fliplr((C2P_latavg1(ind1)+C2P_latstd1(ind1))')],'b','LineStyle','none'); alpha(0.1); hold on
h(2) = fill([lat(ind2),fliplr(lat(ind2))],[(C2P_latavg2(ind2)-C2P_latstd2(ind2))', fliplr((C2P_latavg2(ind2)+C2P_latstd2(ind2))')],'m','LineStyle','none'); alpha(0.1);
h(3) = plot(lat,C2P_latavg1,'-bo'); hold on
%plot(lat,C2P_latmedian1,'-.b','linewidth',2);
h(4) = plot(lat,C2P_latavg2,'-mo')
%plot(lat,C2P_latmedian2,'-.m','linewidth',2);
legend(h([3 4]),'surface','lower EZ');
xlabel('Latitude')
ylabel('C:P [molC:molP]')
title('Latitudinal Average Cellular C:P')
axis tight;

figTitle = 'C2P_lat_avg';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h;


%% C2P vs temperature

%C2P1 = C2P(:,:,1);
%C2P2 = C2P(:,:,2);
indx = ~isnan(C2P);
%indx1 = ~isnan(C2P1);
%indx2 = ~isnan(C2P2);

indx1 = M3d(:,:,1);
indx2 = M3d(:,:,2);

%tempobs1 = tempobs(:,:,1);
%tempobs2 = tempobs(:,:,2);

% seperate data by basin
%iATL = ATL(:,:,1:2) %ATL(indx)
iC2P_ATL = find(~isnan(C2P) & ATL);
iC2P_PAC = find(~isnan(C2P) & PAC);
iC2P_IND = find(~isnan(C2P) & IND);
iC2P_ARC = find(~isnan(C2P) & ARC);
iC2P_MED = find(~isnan(C2P) & MED);

figure;
h1 = plot(tempobs(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
h2 = plot(tempobs(:,:,2),C2P(:,:,2),'.','Color',navy);
% plot(tempobs(iC2P_ATL), C2P(iC2P_ATL),'ro'); hold on
% plot(tempobs(iC2P_PAC), C2P(iC2P_PAC),'ks');
% plot(tempobs(iC2P_IND), C2P(iC2P_IND),'b^');
% plot(tempobs(iC2P_ARC), C2P(iC2P_ARC),'g*');
% plot(tempobs(iC2P_MED), C2P(iC2P_MED),'c>');
% legend('ATL','PAC','IND','ARC','MED')
legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','northwest')
xlabel('Temperature [degC]')
ylabel('C:P [molC:molP]')
title('Cellular C:P vs Temperature')
%axis tight;
grid on
figTitle = 'C2PvsTemp';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h1 h2

%% C2P vs Light (PAR)
figure;
h1 = plot(par.PARobs(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
h2 = plot(par.PARobs(:,:,2),C2P(:,:,2),'.','Color',navy);
legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
xlabel('Photosynthetically Active Radiation (PAR) [$\mu mol photon /m^2 /s$]')
ylabel('C:P [molC:molP]')
title('Cellular C:P vs Light')
%axis tight;
figTitle = 'C2PvsPAR';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h1 h2

%% C2P vs modeled PO4-
figure;
h1 = plot(model.DIP(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
h2 = plot(model.DIP(:,:,2),C2P(:,:,2),'.','Color',navy);
legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
xlabel('Phosphate [$mmol/m^3$]')
ylabel('C:P [molC:molP]')
title('Cellular C:P vs DIP')
grid on
figTitle = 'C2PvsDIP';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h1 h2

%% C2P vs observed nitrate (GLODAP)
figure;
h1 = plot(par.no3raw(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
h2 = plot(par.no3raw(:,:,2),C2P(:,:,2),'.','Color',navy);
legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
xlabel('Nitrate [$mmol/m^3$]')
ylabel('C:P [molC:molP]')
title('Cellular C:P vs DIN')
grid on
figTitle = 'C2PvsDIN';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h1 h2

%% C2P vs cell radius
figure;
%plot(radius(indx),C2P(indx),'k.')
h1 = plot(radius(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
h2 = plot(radius(:,:,2),C2P(:,:,2),'.','Color',navy);
legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
xlabel('radius [$\mu m$]')
ylabel('C:P [molC:molP]')
title('optimal cellular C:P vs cell radius')
grid on
figTitle = 'C2Pvsradius';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h1 h2

%% C2P vs growh rate
figure;
h1 = plot(mu(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
h2 = plot(mu(:,:,2),C2P(:,:,2),'.','Color',navy);
legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
xlabel('Growth rate [1/hr]')  		% need to check this unit
ylabel('C:P [molC:molP]')
title('C:P vs cellular growth rate')
grid on
figTitle = 'C2PvsMu';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h1 h2

%% growth rate vs cell radius
figure;
%plot(radius(indx),C2P(indx),'k.')
h1 = plot(radius(:,:,1),mu(:,:,1),'.','Color',lblue); hold on
h2 = plot(radius(:,:,2),mu(:,:,2),'.','Color',navy);
legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
xlabel('radius [$\mu m$]')
ylabel('growth rate')
title('optimal cellular growth rate vs cell radius')
grid on
figTitle = 'Muvsradius';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h1 h2

iNlim1 = find(LimType(:,:,1) == 0);
iPlim1 = find(LimType(:,:,1) == 1);
iColim1 = find(LimType(:,:,1) == 2 | LimType(:,:,1) == 3);
r1 = radius(:,:,1);
mu1 = mu(:,:,1);
figure
h2 = plot(r1(iColim1),mu1(iColim1),'.','Color',aqua); hold on
h0 = plot(r1(iNlim1),mu1(iNlim1),'.r'); hold on
h1 = plot(r1(iPlim1),mu1(iPlim1),'.b');
legend([h0(1) h1(1) h2(1)],{'N-limited','P-Limited','Co-Limited'},'Location','best')
xlabel('radius [$\mu m$]')
ylabel('growth rate')
title('optimal cellular growth rate vs cell radius for surface layer')
grid on
figTitle = 'Muvsradius_LimType';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h0 h1 h2
% would it be meaningful to multiply cell model's growth rate by biomass in the grid cell?

%%-----------------C2P vs cellular allocation -------------
%% C2P vs PStor
if isfield(model.CellOut,'PStor')
	figure;
	h1 = plot(model.CellOut.PStor(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.PStor(:,:,2),C2P(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	xlabel('PStorage [gP/gCell]')  		% need to check this unit
	ylabel('C:P [molC:molP]')
	title('C:P vs inorganic phosphorus storage')
	grid on
	figTitle = 'C2PvsPStor';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf

	figure;
	h1 = plot(model.DIP(:,:,1),model.CellOut.PStor(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.DIP(:,:,2),model.CellOut.PStor(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('PStorage [gP/gCell]')  		% need to check this unit
	xlabel('DIP [mmol/L]')
	title('inorganic phosphorus storage vs DIP')
	grid on
	figTitle = 'PStorvsDIP';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf
	%PStor vs radius
	figure;
	h1 = plot(model.CellOut.r(:,:,1),model.CellOut.PStor(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.r(:,:,2),model.CellOut.PStor(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('PStorage [gP/gCell]')  		% need to check this unit
	xlabel('radius [$\mu m$]')
	title('Inorganic Phosphorus Storage vs Cell Radius')
	grid on
	figTitle = 'PStorvsradius';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf
end
if isfield(model.CellOut,'PLip')
	figure;
	h1 = plot(model.CellOut.PLip(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.PLip(:,:,2),C2P(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	xlabel('Phospholipid mass fraction')  		% need to check this unit
	ylabel('C:P [molC:molP]')
	title('C:P vs cell phospholipid content')
	grid on
	figTitle = 'C2PvsPLip';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf
end
%% C2P vs E
if isfield(model.CellOut,'E')
	figure;
	h1 = plot(model.CellOut.E(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.E(:,:,2),C2P(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	xlabel('mass fraction of biosynthetic apparatus')  		% need to check this unit
	ylabel('C:P [molC:molP]')
	title('C:P vs Cellular Allocation to Biosynthesis')
	grid on
	figTitle = 'C2PvsE';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
end
%% C2P vs L
if isfield(model.CellOut,'L')
	figure;
	h1 = plot(model.CellOut.L(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.L(:,:,2),C2P(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	xlabel('mass fraction of photosynthetic apparatus')  		% need to check this unit
	ylabel('C:P [molC:molP]')
	title('C:P vs Cellular Allocation to Photosynthesis')
	grid on
	figTitle = 'C2PvsL';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf

	figure;
	h1 = plot(par.PARobs(:,:,1),model.CellOut.L(:,:,1),'.','Color',lblue); hold on
	h2 = plot(par.PARobs(:,:,2),model.CellOut.L(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('mass fraction of photosynthetic apparatus')  		% need to check this unit
	xlabel('Photosynthetically Active Radiation (PAR) [$ \mu mol photon /m^2 /s$]')
	title('Cellular Allocation to Photosynthesis vs Light')
	grid on
	figTitle = 'LvsPAR';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf
end
%% C2P vs A
if isfield(model.CellOut,'A')
	figure;
	h1 = plot(model.CellOut.A(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.A(:,:,2),C2P(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	xlabel('mass fraction of A (or M)')
	ylabel('C:P [molC:molP]')
	title('C:P vs Cellular Allocation to Nutrient Uptake')
	grid on
	figTitle = 'C2PvsA';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf
	figure;
	h1 = plot(model.DIP(:,:,1),model.CellOut.A(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.DIP(:,:,2),model.CellOut.A(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('mass fraction of A (or M)')
	xlabel('DIP [$\mu mol/kg$]') %CHECK UNITS!!
	title('Cellular Allocation to Nutrient Uptake vs DIP')
	grid on
	figTitle = 'AvsDIP';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf
end
%---------------------C:N---------------------------------
%% C2N
figure;
contourf(lon,lat,C2N(:,:,1)); hold on
c=colorbar;
colormap(flipud(summer));
title('Cell Model C:N Uptake Ratio: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'C:N [molC/molN]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'C2Nsurface';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%%% C2N lower EZ
figure; hold on;
imAlpha = ones(size(C2N(:,:,2)));
imAlpha(isnan(C2N(:,:,2))) =0;
imagesc(lon,lat,C2N(:,:,2),'AlphaData',imAlpha)
%contourf(lon,lat,C2N(:,:,2)); hold on
c=colorbar;
colormap(flipud(summer));
title('Cell Model C:N Uptake Ratio: Lower EZ','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'C:N [molC/molN]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'C2N_Z2';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% N2P surface
figure;
contourf(1:2:360,-89:2:89,N2P(:,:,1)); hold on
c=colorbar;
colormap(flipud(summer));
title('Cell Model N:P Uptake Ratio: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'N:P [molN/molP]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'N2Psurface';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%%%% N2P lowwer EZ
figure;
contourf(lon,lat,N2P(:,:,2)); hold on
c=colorbar;
colormap(flipud(summer));
title('Cell Model N:P Uptake Ratio: Lower EZ','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'N:P [molN/molP]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'N2P_Z2';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% LimType
% Surface Limitation Type: 0=N-Lim; 1=P-Lim; 2&3=Co-Lim
Lim_cmap = [1, 0, 0; ...
    0, 0, 1; ...
    0, 1, 1; ...
    0, 1, 1];
figure; hold on
Llevs= [0,1,2,3];
%plt=pcolor(lon,lat,LimType(:,:,1));
%set(plt,'EdgeColor','none');

imAlpha = ones(size(LimType(:,:,1)));
imAlpha(isnan(LimType(:,:,1))) =0;
imagesc(lon,lat,LimType(:,:,1),'AlphaData',imAlpha)
cb=colorbar('Ticks',[0,1,2,3],'TickLabels',{'N-Lim','P-Lim','Co-Lim','Co-Lim-alt'});
colormap(Lim_cmap);
ylabel(cb,'Limitation Type');
axis tight
%colormap(parula(length(Llevs)))

title('Phytoplankton Nutrient Limitation Type: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
grid off

figTitle = 'LimType_surf';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%%%% LimType lower EZ
%plt = pcolor(lon,lat,LimType(:,:,2));
%set(plt,'EdgeColor','none');

imAlpha = ones(size(LimType(:,:,2)));
imAlpha(isnan(LimType(:,:,2))) =0;
imagesc(lon,lat,LimType(:,:,2),'AlphaData',imAlpha)
cb=colorbar('Ticks',[0,1,2,3],'TickLabels',{'N-Lim','P-Lim','Co-Lim','Co-Lim-alt'});
colormap(Lim_cmap);
ylabel(cb,'Limitation Type');
axis tight

title('Phytoplankton Nutrient Limitation Type: Lower EZ','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
grid off

figTitle = 'LimType_Z2';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')


%----------radius-----------------------------
%% radius
Zlevs = [0:0.2:6];
figure;
contourf(lon,lat,radius(:,:,1),Zlevs); hold on
c=colorbar;
colormap(flipud(summer));
[CC,hh] = contour(lon,lat,radius(:,:,1),[0.2,2],'k');
clabel(CC,hh,'FontName','Times');
%caxis([Zlevs(1) Zlevs(end)]);
%cmocean('matter',length(Zlevs)-1);

title('Cell radius: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'radius [\mu m]');
axis tight; grid off

figTitle = 'radius_surface';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')


%Zlevs = [0.25:0.25:2.75];
figure;
contourf(lon,lat,radius(:,:,2),Zlevs); hold on
c=colorbar;
colormap(flipud(summer));
[CC,hh] = contour(lon,lat,radius(:,:,2),[0.2,2],'k');
clabel(CC,hh,'FontName','Times');
%caxis([Zlevs(1) Zlevs(end)]);
%cmocean('matter',length(Zlevs)-1);
title('Cell radius: Lower EZ','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'radius [\mu m]');
axis tight; grid off

figTitle = 'radius_Z2';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')


%% radius as a function of latitude
r_latavg1 = mean(radius(:,:,1),2,'omitnan');
r_latavg2 = mean(radius(:,:,2),2,'omitnan');
r_latavg = mean(radius(:,:,1:2),[2 3],'omitnan');
r_latstd1 = std(radius(:,:,1),0,2,'omitnan');
r_latstd2 = std(radius(:,:,2),0,2,'omitnan');
ind1 = ~isnan(r_latavg1);
ind2 = ~isnan(r_latavg2);

figure;
% plot +/-1 standard deviation
h(1) = fill([lat(ind1),fliplr(lat(ind1))],[(r_latavg1(ind1)-r_latstd1(ind1))', fliplr((r_latavg1(ind1)+r_latstd1(ind1))')],'b','LineStyle','none'); alpha(0.1); hold on
h(2) = fill([lat(ind2),fliplr(lat(ind2))],[(r_latavg2(ind2)-r_latstd2(ind2))', fliplr((r_latavg2(ind2)+r_latstd2(ind2))')],'m','LineStyle','none'); alpha(0.1);
h(3) = plot(lat,r_latavg1,'-bo'); hold on
h(4) = plot(lat,r_latavg2,'-mo')
legend(h([3 4]),'surface','lower EZ');
xlabel('Latitude')
ylabel('Radius [um]')
title('Latitudinal Average Modeled Cell Radius')
axis tight;

figTitle = 'radius_lat_avg';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% Growth rate
figure;
contourf(lon,lat,mu(:,:,1)); hold on
c=colorbar;
colormap(flipud(summer));
title('Cell growth rate: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'growth rate [ ]');
axis tight; grid off

figTitle = 'mu_surf';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

close all;
