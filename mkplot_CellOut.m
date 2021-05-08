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

%%
close all;

GridVer = 90;
RunVer = 'Tv4_PCCellv3b_DOC0.25_DOP0';

%ver = datestr(now,'mmmdd');
%model output directory
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/';
figDir = strcat(outputDir,'FIGS_PCCellv3b_DOC0.25_DOP0/');

%outPath='/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/FIGS_PCCell_DOC0.25_DOP0/';
outPath = figDir;

% load model output fields
%fname='/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/Tv4_PCCell_DOC0.25_DOP0.mat';
fname = strcat(outputDir, RunVer, '.mat');
load(fname);

% load optimal parameter values
%fxhat = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/Tv4_PCCell_DOC0.25_DOP0_xhat.mat';
fxhat = strcat(outputDir, RunVer,'_xhat.mat');
load(fxhat);


%% rename variables
C2P = data.CellOut.C2P;
N2P = data.CellOut.N2P;
C2N = data.CellOut.C2N;
radius = data.CellOut.r;
LimType = data.CellOut.LimType;

% set land points to NaN
C2P(find(C2P==0))=NaN;
C2N(find(C2N==0))=NaN;
N2P(find(N2P==0))=NaN;
radius(find(radius==0))=NaN;

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


%% C2P surface plot
%Zlevs = [100:20:300];
figure;
contourf(lon,lat,C2P(:,:,1)); hold on
cb=colorbar;
colormap(flipud(summer))
title('Cell Model C:P Uptake Ratio: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'C:P [gC/gP]');

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
ylabel(c,'C:P [gC/gP]');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'C2P_Z2';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% C2N
figure;
contourf(lon,lat,C2N(:,:,1)); hold on
c=colorbar;
colormap(flipud(summer));
title('Cell Model C:N Uptake Ratio: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'C:N [gC/gN]');
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
ylabel(c,'C:N [gC/gN]');
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
ylabel(c,'N:P [gN/gP]');
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
ylabel(c,'N:P [gN/gP]');
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
