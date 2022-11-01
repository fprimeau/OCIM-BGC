clc; clear all; close all
on   = true  ;
off  = false ;

%ver = datestr(now,'mmmdd');
%RunVer = 'Tv4_PCCellv9c_onlyStor_DOC0.25_DOP0';
%RunVer = 'testNPP_Tv4_PCCella1e-8b1e-3_DOC0.25_DOP0';
%RunVer = 'testPobs_CTL_He_PCCella1b1_DOC0.25_DOP0';
%RunVer = 'testPobs_Tv4_PCCella1b1_DOC0.25_DOP0'
%RunVer = 'testNPP_CTL_He_PCCella1e-4bfix_DOC0.25_DOP0'
RunVer = 'optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0'

GridVer  = 91  ;
operator = 'A' ;

%model output directory
%outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/';
%outputDir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/', GridVer);
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/';
figDir = strcat(outputDir,'FIGS_optC_Cell/');
%figDir = strcat(outputDir,'FIGS_PCCellv9_DOC0.25_DOP0/v9c_onlyPStor_');
%figDir = strcat(outputDir,'FIGS_PCCell_Pobs/a1b1_');
%figDir = strcat(outputDir,'FIGS_testPobs_PCCell/a1b1_');
%figDir = strcat(outputDir,'FIGS_testNPP_PCCellfixb/a1e-4_');
outPath = figDir;

% load model output fields
fname = strcat(outputDir, RunVer, '.mat');
load(fname);
model = data;

% load optimal parameter values
fxhat = strcat(outputDir, RunVer,'_xhat.mat');
par.fxhat = fxhat;
par.fxhatload = fxhat;
load(fxhat);

par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function
par.LoadOpt = on;
par.dynamicP = off ;

%-------------load data and set up parameters---------------------
SetUp ;
xhat

par = SetPar(par) ;

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

%Colors ROYGBIV
%gradient: aqua, teal, darkgreen
lblue = [0 191/255 255/255];
navy = [ 0 0 128/255];

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
lon = grd.xt;
lat = grd.yt;
% if GridVer == 90;
%     lat = [-89:2:89];
%     lon = [1:2:360];
% end
% if GridVer == 91;
%     lat = [-90:2:90];
%     lon = [1:2:360];
% end

%% annotation for parameter values
dim = [0.0 0.6 0.1 0.2];
dim = [0.1199 0.695 0.1 0.2]
%parstr = 'Cell Model Parameters';
parstr = {'Cell Model Parameters'};
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
    parstr{kk,1} = ['alphaS=' num2str(xhat.alphaS)];
end
if isfield(xhat,'gammaDNA')
	kk=kk+1;
    parstr{kk,1} = ['gammaDNA=' num2str(xhat.gammaDNA)];
end
parstr

nzo = par.nzo; % 2 eupotic layers


%% ------ Limitation types --------
iprod = find(M3d(:,:,1:2));
uLimTypes = unique(model.CellOut.LimType(iprod));

fprintf('# of N limited points : %d \n',length(find(LimType(iprod) ==0)));
fprintf('# of P limited points : %d \n',length(find(LimType(iprod) ==1)));
fprintf('# of Colimited2 points: %d \n',length(find(LimType(iprod) ==2)));
fprintf('# of Colimited3 points: %d \n',length(find(LimType(iprod) ==3)));

fprintf('range of C:P values is %6.1f to %6.1f \n',min(C2P(iprod)),max(C2P(iprod)))

fprintf('mean of C:P is %6.1f \n',mean(C2P(iprod)))
fprintf('std of C:P  is %6.1f \n\n',std(C2P(iprod)))

fprintf('range of radius is      %6.2f to %6.2f um \n',min(radius(iprod)),max(radius(iprod)))

fprintf('range of growth rate is %6.3f to %6.3f hr^-1 \n',min(mu(iprod)),max(mu(iprod)))
fprintf('range of growth rate is %6.3f to %6.3f day^-1 \n',min(mu(iprod)*24),max(mu(iprod)*24))



%% ---- average C2P  production weighted--------
DIP  = model.DIP(iwet) ;

LAM        = 0*M3d;
LAM(:,:,1) = (par.npp1.^par.beta).*par.Lambda(:,:,1);
LAM(:,:,2) = (par.npp2.^par.beta).*par.Lambda(:,:,2);
Lam          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
clear LAM;
%--------------- calculate primary production --------------------
G        = M3d*0        ;
G(iwet)  = par.alpha*Lam*DIP  ; % primary production [unit: mmol P/m^3/s]

Int_CNPP = G(:,:,1:nzo).*grd.DZT3d(:,:,1:nzo).*C2P(:,:,1:nzo)*12;
CNPP = Int_CNPP*spa*1e-3 ; % convert production from mg C/m^3/s to gC/m^2/year;

% prod weighted C:P
WeightNPP = CNPP.*dAt(:,:,1:nzo)/nansum(CNPP.*dAt(:,:,1:nzo),'all'); %gC/yr
C2Pavg = nansum(C2P(:,:,1:nzo).*WeightNPP,'all');
fprintf('Average C:P uptake in Euphotic Layers (NPP weighted) is %4.2f \n',C2Pavg);

%% ----- cell model based NPP ----------
% POC*cell model growth rate (mu)
POC = model.POC; 		% [unit??: mmol C/m^3]
cellNPP = POC(:,:,1:nzo).*mu(:,:,1:nzo); %[unit??: mmol C/m^3/hr]
cellNPP = cellNPP.*grd.DZT3d(:,:,1:nzo)./1000.*24*365; % unit = mol C/m^3/yr

figure; hold on
imAlpha = ones(size(cellNPP(:,:,1)));
imAlpha(isnan(cellNPP(:,:,1))) =0;
imagesc(lon,lat,cellNPP(:,:,1),'AlphaData',imAlpha)
c=colorbar;
ax1 = gca;
%ax1.CLim(1) = 0;
%cmap = colormap(parula);
%cmap(1,:) = [1 0 0];
%colormap(cmap)
%[CC,hh] = contour(lon,lat,mu(:,:,1),[1 1],'k');
%clabel(CC,hh,'FontName','Times');
title('Cell Model growth rate * POC: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'production [mol C/m^2/yr]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'cellNPP_surf_molperm2yr';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')


% plot sum of surface and box 2
cellNPPtotal = sum(cellNPP,3);

figure; hold on
imAlpha = ones(size(cellNPPtotal(:,:,1)));
imAlpha(isnan(cellNPPtotal(:,:,1))) =0;
imagesc(lon,lat,cellNPPtotal(:,:,1),'AlphaData',imAlpha)
c=colorbar;
ax1 = gca;
%ax1.CLim(1) = 0;
%cmap = colormap(parula);
%cmap(1,:) = [1 0 0];
%colormap(cmap)
%[CC,hh] = contour(lon,lat,cellNPPtotal,[0:25:max(cellNPPtotal,[],'all')],'k');
[CC,hh] = contour(lon,lat,cellNPPtotal,[0 10 30 50 100],'k');
clabel(CC,hh,'FontName','Times');
title('Cell Model growth rate * POC: EZ total','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'production [mol C/m^2/yr]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'cellNPP_total_molperm2yr_contourlow';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')


%% -----------
% make latitude zone masks: for 0-15, 15-30, 30-45, and 45-60 degrees latitude (inc. both hemispheres)
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

junk4 =  M3d(:,:,1)*0;
junk4(1:length(grd.yt),:) = 1;
msk_polar = junk4-junk3;

clear junk1 junk2 junk3 junk4



%% model POC/POP
POC_POP = M3d*NaN;
indx1=find(model.POC(iwet)>0 & model.POP(iwet)>0);
POC_POP(iwet) = -1;
POC_POP(iwet(indx1)) = model.POC(iwet(indx1))./model.POP(iwet(indx1));
figure;
for ii = 1:4
	subplot(2,2,ii)
	contourf(lon,lat,POC_POP(:,:,ii)); hold on
	cb=colorbar;
	colormap(flipud(summer))
	title(['Cell Model POC:POP - Z level ' num2str(ii)],'Fontsize',14);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'POC:POP [ ]');
end
figTitle = 'POC2POP_z1-4' % ['POC2POP_Z' num2str(ii)];
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
%why are there negative values?
%need to make dry poits nan

%% POC:POP & DOC:DOP scatter plot
%% plot Pstorage & PLip?

%% C2P surface plot
%Zlevs = [100:20:300];
figure; hold on;
%contourf(lon,lat,C2P(:,:,1)); hold on
imAlpha = ones(size(C2P(:,:,1)));
imAlpha(isnan(C2P(:,:,1))) =0;
imagesc(lon,lat,C2P(:,:,1),'AlphaData',imAlpha)
cb=colorbar;
colormap(flipud(summer))
%cmap = cmocean('-curl','pivot',0);
%colormap(cmap)
[CC,hh] = contour(lon,lat,C2P(:,:,1),[106 106],'k');
clabel(CC,hh,'FontName','Times');
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
[CC,hh] = contour(lon,lat,C2P(:,:,2),[106 106],'k');
clabel(CC,hh,'FontName','Times');
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
h(4) = plot(lat,C2P_latavg2,'-mo');
%plot(lat,C2P_latmedian2,'-.m','linewidth',2);
legend(h([3 4]),'surface','lower EZ');
xlabel('Latitude')
ylabel('C:P [molC:molP]')
title('Latitudinal Average Cellular C:P')
axis tight;
grid on;

figTitle = 'C2P_lat_avg';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h;


%% C2P as a function of latitude ATL v PAC
% seperate data by basin
%iATL = ATL(:,:,1:2) %ATL(indx)
iC2P_ATL = find(~isnan(C2P) & ATL);
iC2P_PAC = find(~isnan(C2P) & PAC);
iC2P_IND = find(~isnan(C2P) & IND);
iC2P_ARC = find(~isnan(C2P) & ARC);
iC2P_MED = find(~isnan(C2P) & MED);

M.ATL = ATL(:,:,1:2);
M.ATL(M.ATL==0)=NaN;
M.PAC = PAC(:,:,1:2);
M.PAC(M.PAC==0)=NaN;
M.IND = IND(:,:,1:2);
M.IND(M.IND==0)=NaN;
M.ARC = ARC(:,:,1:2);
M.ARC(M.ARC==0)=NaN;

C2P_ATL = C2P(:,:,1:2).*M.ATL;
C2P_ATL_latavg = mean(C2P_ATL,[2 3],'omitnan');
C2P_ATL_latstd = std(C2P_ATL,0,[2 3],'omitnan');

C2P_PAC = C2P(:,:,1:2).*M.PAC;
C2P_PAC_latavg = mean(C2P_PAC,[2 3],'omitnan');
C2P_PAC_latstd = std(C2P_PAC,0,[2 3],'omitnan');

C2P_IND = C2P(:,:,1:2).*M.IND;
C2P_IND_latavg = mean(C2P_IND,[2 3],'omitnan');
C2P_IND_latstd = std(C2P_PAC,0,[2 3],'omitnan');

C2P_ARC = C2P(:,:,1:2).*M.ARC;
C2P_ARC_latavg = mean(C2P_ARC,[2 3],'omitnan');
C2P_ARC_latstd = std(C2P_ARC,0,[2 3],'omitnan');

ind1 = ~isnan(C2P_ATL_latavg);
ind2 = ~isnan(C2P_PAC_latavg);
ind3 = ~isnan(C2P_IND_latavg);
ind4 = ~isnan(C2P_ARC_latavg);


figure;
h(1) = fill([lat(ind1),fliplr(lat(ind1))],[(C2P_ATL_latavg(ind1)-C2P_ATL_latstd(ind1))', fliplr((C2P_ATL_latavg(ind1)+C2P_ATL_latstd(ind1))')],'m','LineStyle','none'); alpha(0.1); hold on
h(2) = fill([lat(ind2),fliplr(lat(ind2))],[(C2P_PAC_latavg(ind2)-C2P_PAC_latstd(ind2))', fliplr((C2P_PAC_latavg(ind2)+C2P_PAC_latstd(ind2))')],'b','LineStyle','none'); alpha(0.1);
h(3) = fill([lat(ind3),fliplr(lat(ind3))],[(C2P_IND_latavg(ind3)-C2P_IND_latstd(ind3))', fliplr((C2P_IND_latavg(ind3)+C2P_IND_latstd(ind3))')],colors.limegreen,'LineStyle','none'); alpha(0.1);
h(4) = fill([lat(ind4),fliplr(lat(ind4))],[(C2P_ARC_latavg(ind4)-C2P_ARC_latstd(ind4))', fliplr((C2P_ARC_latavg(ind4)+C2P_ARC_latstd(ind4))')],colors.lblue,'LineStyle','none'); alpha(0.1);

h(5) = plot(lat,C2P_ATL_latavg,'-mo'); hold on
%plot(lat,C2P_latmedian1,'-.b','linewidth',2);
h(6) = plot(lat,C2P_PAC_latavg,'-bo');
h(7) = plot(lat,C2P_IND_latavg,'-o','Color',colors.limegreen);
h(8) = plot(lat,C2P_ARC_latavg,'-o','Color',colors.lblue);
legend(h([5 6 7 8]),'Atlantic','Pacific','Indian','Arctic');
xlabel('Latitude')
ylabel('C:P [molC:molP]')
title('Zonal Average Cellular C:P')
axis tight; grid on
ylim([0 400])
figTitle = 'C2P_lat_avg_basin';
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
%{
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
%}

figure;
h5 = plot(radius(:,:,1).*msk_polar,C2P(:,:,1).*msk_polar,'.c'); hold on
h5 = plot(radius(:,:,2).*msk_polar,C2P(:,:,2).*msk_polar,'.c'); hold on
h1 = plot(radius(:,:,1).*msk_tropical,C2P(:,:,1).*msk_tropical,'.','Color','red'); hold on
h1 = plot(radius(:,:,2).*msk_tropical,C2P(:,:,2).*msk_tropical,'.','Color','red'); hold on
h2 = plot(radius(:,:,1).*msk_subtro,C2P(:,:,1).*msk_subtro,'.m'); hold on
h2 = plot(radius(:,:,2).*msk_subtro,C2P(:,:,2).*msk_subtro,'.m'); hold on
h3 = plot(radius(:,:,1).*msk_subtro_subpo,C2P(:,:,1).*msk_subtro_subpo,'.','Color',colors.limegreen); hold on
h3 = plot(radius(:,:,2).*msk_subtro_subpo,C2P(:,:,2).*msk_subtro_subpo,'.','Color',colors.limegreen); hold on
h4 = plot(radius(:,:,1).*msk_subpolar,C2P(:,:,1).*msk_subpolar,'.b'); hold on
h4 = plot(radius(:,:,2).*msk_subpolar,C2P(:,:,2).*msk_subpolar,'.b'); hold on
legend([h1(1) h2(1) h3(1) h4(1) h5(1)],{'Tropical (0-15)','Subtropical (15-30)','midlat (30-45)','subpolar (45-60)','polar (60-90)'},'Location','best')
xlabel('radius [$\mu m$]')
ylabel('C:P [molC:molP]')
title('optimal cellular C:P vs cell radius')
grid on
figTitle = 'C2Pvsradius_latzones';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h1 h2 h3 h4 h5

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

if uLimTypes > 2
	iNlim1 = find(LimType(:,:,1) == 0);
	iPlim1 = find(LimType(:,:,1) == 1);
	iColim1 = find(LimType(:,:,1) == 2 | LimType(:,:,1) == 3);
	r1 = radius(:,:,1);
	mu1 = mu(:,:,1);
	figure
	h2 = plot(r1(iColim1),mu1(iColim1),'.','Color',colors.aqua); hold on
	h0 = plot(r1(iNlim1),mu1(iNlim1),'.r'); hold on
	h1 = plot(r1(iPlim1),mu1(iPlim1),'.b');
	legend([h0(1) h1(1) h2(1)],{'N-limited','P-Limited','Co-Limited'},'Location','best')
	%legend([h0(1) h1(1)],{'N-limited','P-Limited'},'Location','best')
	xlabel('radius [$\mu m$]')
	ylabel('growth rate')
	title('optimal cellular growth rate vs cell radius for surface layer')
	grid on
	figTitle = 'Muvsradius_LimType';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h0 h1 h2
	% would it be meaningful to multiply cell model's growth rate by biomass in the grid cell?
end

%%-----------------C2P vs cellular allocation -------------
%% C2P vs PStor
if isfield(model.CellOut,'PStor')
	figure;
	h1 = plot(model.CellOut.PStor(:,:,1),C2P(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.PStor(:,:,2),C2P(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	xlabel('PStorage [mol P]')  		% need to check this unit
	ylabel('C:P [molC:molP]')
	title('C:P vs inorganic phosphorus storage')
	grid on
	figTitle = 'C2PvsPStor';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf

%{
	figure;
	h1 = plot(model.DIP(:,:,1),model.CellOut.PStor(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.DIP(:,:,2),model.CellOut.PStor(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('PStorage [mol P]')  		% need to check this unit
	xlabel('DIP [mmol/L]')
	title('inorganic phosphorus storage vs DIP')
	grid on
	figTitle = 'PStorvsDIP';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf

	% MSKS.ATL(MSKS.ATL==0)=NaN;
	% MSKS.PAC(MSKS.PAC==0)=NaN;
	% MSKS.IND(MSKS.IND==0)=NaN;
	% MSKS.MED(MSKS.MED==0)=NaN;
	% MSKS.ARC(MSKS.ARC==0)=NaN;
	figure;
	h1 = plot(model.DIP(:,:,1).*MSKS.ATL(:,:,1),model.CellOut.PStor(:,:,1).*MSKS.ATL(:,:,1),'.','Color','red'); hold on
	h1 = plot(model.DIP(:,:,2).*MSKS.ATL(:,:,2),model.CellOut.PStor(:,:,2).*MSKS.ATL(:,:,2),'.','Color','red'); hold on
	h2 = plot(model.DIP(:,:,1).*MSKS.PAC(:,:,1),model.CellOut.PStor(:,:,1).*MSKS.PAC(:,:,1),'.','Color','blue'); hold on
	h2 = plot(model.DIP(:,:,2).*MSKS.PAC(:,:,2),model.CellOut.PStor(:,:,2).*MSKS.PAC(:,:,2),'.','Color','blue'); hold on
	h3 = plot(model.DIP(:,:,1).*MSKS.IND(:,:,1),model.CellOut.PStor(:,:,1).*MSKS.IND(:,:,1),'.','Color',colors.limegreen); hold on
	h3 = plot(model.DIP(:,:,2).*MSKS.IND(:,:,2),model.CellOut.PStor(:,:,2).*MSKS.IND(:,:,2),'.','Color',colors.limegreen); hold on
	h4 = plot(model.DIP(:,:,1).*MSKS.ARC(:,:,1),model.CellOut.PStor(:,:,1).*MSKS.ARC(:,:,1),'.','Color',lblue); hold on
	h4 = plot(model.DIP(:,:,2).*MSKS.ARC(:,:,2),model.CellOut.PStor(:,:,2).*MSKS.ARC(:,:,2),'.','Color',lblue); hold on
	h5 = plot(model.DIP(:,:,1).*MSKS.MED(:,:,1),model.CellOut.PStor(:,:,1).*MSKS.MED(:,:,1),'.','Color','m'); hold on
	h5 = plot(model.DIP(:,:,2).*MSKS.MED(:,:,2),model.CellOut.PStor(:,:,2).*MSKS.MED(:,:,2),'.','Color','m'); hold on
	legend([h1(1) h2(1) h3(1) h4(1) h5(1)],{'Atlantic','Pacific','Indian','Arctic','Mediterranean'},'Location','best')
	ylabel('PStorage [molP]')
	xlabel('DIP [mmol/L]')
	title('inorganic phosphorus storage vs DIP')
	grid on
	figTitle = 'PStorvsDIP_oceans';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2 h3 h4 h5
	clf
%}
	figure;
	h1 = plot(model.DIP(:,:,1).*msk_tropical,model.CellOut.PStor(:,:,1).*msk_tropical,'.','Color','red'); hold on
	h1 = plot(model.DIP(:,:,2).*msk_tropical,model.CellOut.PStor(:,:,2).*msk_tropical,'.','Color','red'); hold on
	h2 = plot(model.DIP(:,:,1).*msk_subtro,model.CellOut.PStor(:,:,1).*msk_subtro,'.m'); hold on
	h2 = plot(model.DIP(:,:,2).*msk_subtro,model.CellOut.PStor(:,:,2).*msk_subtro,'.m'); hold on
	h3 = plot(model.DIP(:,:,1).*msk_subtro_subpo,model.CellOut.PStor(:,:,1).*msk_subtro_subpo,'.','Color',colors.aqua); hold on
	h3 = plot(model.DIP(:,:,2).*msk_subtro_subpo,model.CellOut.PStor(:,:,2).*msk_subtro_subpo,'.','Color',colors.aqua); hold on
	h4 = plot(model.DIP(:,:,1).*msk_subpolar,model.CellOut.PStor(:,:,1).*msk_subpolar,'.b'); hold on
	h4 = plot(model.DIP(:,:,2).*msk_subpolar,model.CellOut.PStor(:,:,2).*msk_subpolar,'.b'); hold on
	legend([h1(1) h2(1) h3(1) h4(1)],{'Tropical','Subtropical','subtrop-subpolar','subpolar'},'Location','best')
	ylabel('PStorage [molP]')  		% need to check this unit
	xlabel('DIP [mmol/L]')
	title('inorganic phosphorus storage vs DIP')
	grid on
	figTitle = 'PStorvsDIP_latzones';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2 h3 h4
	clf

	%PStor vs radius
	figure;
	h1 = plot(model.CellOut.r(:,:,1),model.CellOut.PStor(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.r(:,:,2),model.CellOut.PStor(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('PStorage [mol P]')  		% need to check this unit
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
	xlabel('Phospholipid content [mol P]')  		% need to check this unit
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
figure; hold on
imAlpha = ones(size(N2P(:,:,1)));
imAlpha(isnan(N2P(:,:,1))) =0;
imagesc(lon,lat,N2P(:,:,1),'AlphaData',imAlpha)
%contourf(1:2:360,-89:2:89,N2P(:,:,1)); hold on
c=colorbar;
colormap(flipud(summer));
[CC,hh] = contour(lon,lat,N2P(:,:,1),[16 16],'k');
clabel(CC,hh,'FontName','Times');
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
imAlpha = ones(size(N2P(:,:,2)));
imAlpha(isnan(N2P(:,:,2))) =0;
imagesc(lon,lat,N2P(:,:,2),'AlphaData',imAlpha)
%contourf(lon,lat,N2P(:,:,2)); hold on
c=colorbar;
colormap(flipud(summer));
title('Cell Model N:P Uptake Ratio: Lower EZ','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'N:P [molN/molP]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis xy
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
Clevs= [0,1,2,3]+1;
%plt=pcolor(lon,lat,LimType(:,:,1));
%set(plt,'EdgeColor','none');

imAlpha = ones(size(LimType(:,:,1)));
imAlpha(isnan(LimType(:,:,1))) =0;
image(lon,lat,LimType(:,:,1)+1,'AlphaData',imAlpha)
cb=colorbar('Ticks',Clevs,'TickLabels',{'N-Lim','P-Lim','Co-Lim','Co-Lim-alt'});
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
image(lon,lat,LimType(:,:,2)+1,'AlphaData',imAlpha)
cb=colorbar('Ticks',Clevs,'TickLabels',{'N-Lim','P-Lim','Co-Lim','Co-Lim-alt'});
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
%Zlevs = [0:0.5:20];
Zlevs = 20;
figure;
contourf(lon,lat,radius(:,:,1),Zlevs); hold on
c=colorbar;
colormap(flipud(summer));
[CC,hh] = contour(lon,lat,radius(:,:,1),[0.2,1,10],'k');
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
[CC,hh] = contour(lon,lat,radius(:,:,2),[0.2,1,10],'k');
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
grid on;

figTitle = 'radius_lat_avg';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% ------ Growth rate --------------
figure;
contourf(lon,lat,mu(:,:,1)); hold on
c=colorbar;
colormap(flipud(summer));
title('Cell growth rate: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'growth rate [1/hr]');
axis tight; grid off

figTitle = 'mu_surf';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

% surface (image)
figure; hold on
imAlpha = ones(size(mu(:,:,1)));
imAlpha(isnan(mu(:,:,1))) =0;
imagesc(lon,lat,mu(:,:,1),'AlphaData',imAlpha)
c=colorbar;
ax1 = gca;
ax1.CLim(1) = 0;
cmap = colormap(parula);
cmap(1,:) = [1 0 0];
colormap(cmap)
%[CC,hh] = contour(lon,lat,mu(:,:,1),[1 1],'k');
%clabel(CC,hh,'FontName','Times');
title('Cell Model growth rate: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'growth rate [1/hr]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'mu_surf_image_parula';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

% growth rate z2
figure; hold on
imAlpha = ones(size(mu(:,:,2)));
imAlpha(isnan(mu(:,:,2))) =0;
imagesc(lon,lat,mu(:,:,2),'AlphaData',imAlpha)
c=colorbar;
ax1 = gca;
ax1.CLim(1) =0;
cmap = colormap(parula);
cmap(1,:) = [1 0 0];
colormap(cmap);
%[CC,hh] = contour(lon,lat,mu(:,:,1),[1 1],'k');
%clabel(CC,hh,'FontName','Times');
title('Cell Model growth rate: lower EZ','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'growth rate [1/hr]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'mu_z2_image_parula';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

close all



%% ------ PStor map --------------
molarP = 31.0;
figure;
contourf(lon,lat,model.CellOut.PStor(:,:,1)./molarP); hold on
c=colorbar;
colormap(flipud(summer));
title('P Storage: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'P Storage [mol P]');
axis tight; grid off

figTitle = 'PStor_surf';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

figure;
contourf(lon,lat,model.CellOut.PStor(:,:,2)./molarP); hold on
c=colorbar;
colormap(flipud(summer));
title('P Storage: Z2','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'P Storage [mol P]');
axis tight; grid off

figTitle = 'PStor_Z2';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

close all;

%% QP
if isfield(model.CellOut,'QP')
	figure;
	h1 = plot(model.DIP(:,:,1),model.CellOut.QP(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.DIP(:,:,2),model.CellOut.QP(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('Cellular P quota [mol P]')  		% need to check this unit
	xlabel('DIP [mmol/L]')
	title('Cellular P content vs DIP')
	grid on
	figTitle = 'QPvsDIP';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf

	figure;
	PStorfrac = (model.CellOut.PStor)./model.CellOut.QP;
	h1 = plot(model.DIP(:,:,1),PStorfrac(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.DIP(:,:,2),PStorfrac(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('PStorage / total PQuota [molP/molP]')
	xlabel('DIP [mmol/L]')
	title('fraction of cellular P quota in P storage vs DIP')
	grid on
	figTitle = 'PStor_QPvsDIP';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf

	figure;
	h1 = plot(model.CellOut.QP(:,:,1),model.CellOut.QC(:,:,1),'.','Color',lblue); hold on
	h2 = plot(model.CellOut.QP(:,:,2),model.CellOut.QC(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('Cellular C quota [mol C]')
	xlabel('Cellular P Quota [mol P]')
	title('Cellular C:P')
	grid on
	figTitle = 'QCvsQP';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf

	figure;
	QPnostor = model.CellOut.QP - model.CellOut.PStor;
	h1 = plot(QPnostor(:,:,1),model.CellOut.QC(:,:,1),'.','Color',lblue); hold on
	h2 = plot(QPnostor(:,:,2),model.CellOut.QC(:,:,2),'.','Color',navy);
	legend([h1(1) h2(1)],{'Surface','Lower EZ'},'Location','best')
	ylabel('Cellular C quota [mol C]')
	xlabel('Cellular P Quota excluding P storage [mol P]')
	title('Cellular C:P excluding P storage')
	grid on
	figTitle = 'QCvsQPnostor';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h1 h2
	clf

end
