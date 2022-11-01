clc; clear all; close all
on = true;      off = false;
spd  = 24*60^2; spa  = 365*spd;

%RunVer = 'Tv4_PCCellv8_DOC0.25_DOP0';
%RunVer = 'optGM15_CTL_He_PC_DOC0.25_DOP0'
%RunVer = 'testC2PTempP_CTL_He_PC_DOC0.25_DOP0'
RunVer = 'optC_GM15_CTL_He_PC_DOC0.25_DOP0' ;

%model output directory
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/';
%outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK91/';
%figDir = strcat(outputDir,'FIGS_testC2PTempP/v1_lcc_');
figDir = strcat(outputDir,'FIGS_optC_GM15/');
outPath = figDir;

% load model output fields
fname = strcat(outputDir, RunVer, '.mat');
load(fname);
model = data;

% load optimal parameter values
fxhat = strcat(outputDir, RunVer,'_xhat.mat');
load(fxhat);

GridVer  = 91  ;
operator = 'A' ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = off; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function
par.dynamicP = off;

%-------Define some colors ------
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

set(groot,'defaultAxesFontName','Times',...
    'defaultAxesFontSize',14,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','on',...
    'defaultAxesYMinorTick','on');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Times',...
    'defaultTextInterpreter','latex');

%-------------load data and set up parameters---------------------
SetUp ;
xhat
%{
% addpath according to opterating system
if ismac
    addpath('~/Dropbox/myfunc'     )
    addpath('~/Documents/DATA/'    )
    addpath('~/Documents/DATA/OCIM')
else
    addpath('/DFS-L/DATA/primeau/weilewang/DATA')
    addpath('/DFS-L/DATA/primeau/weilewang/my_func')
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
end
format long
%
Cmodel  = on ;
Omodel  = off ;
Simodel = off ;
Cellmodel = on;
pscale  = 0.0 ;
cscale  = 0.25 ;
%
GridVer   = 90 ;
operator = 'A' ;
if GridVer == 90
    TRdivVer = 'Tv4' ;
elseif GridVer == 91
    switch(operator)
      case 'A'
        TRdivVer = 'CTL_He'   ;
      case 'B'
        TRdivVer = 'CTL_noHe' ;
      case 'C'
        TRdivVer = 'KiHIGH_He'   ;
      case 'D'
        TRdivVer = 'KiHIGH_noHe' ;
      case 'E'
        TRdivVer = 'KvHIGH_KiLOW_He'  ;
      case 'F'
        TRdivVer = 'KvHIGH_KiLOW_noHe';
      case 'G'
        TRdivVer = 'KiLOW_He'   ;
      case 'H'
        TRdivVer = 'KiLOW_noHe' ;
      case 'I'
        TRdivVer = 'KvHIGH_He'  ;
      case 'J'
        TRdivVer = 'KvHIGH_noHe';
      case 'K'
        TRdivVer = 'KvHIGH_KiHIGH_noHe';
    end
end

% save results
% ATTENTION: please change this directory to where you wanna
if ismac
    input_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer);
elseif isunix
    % input_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
    % 'MSK%2d/'],GridVer);
    %input_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
    %                    'MSK%2d/PME4DICALK/'],GridVer);
	input_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/', GridVer);
    % input_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/' ...
                        % 'TempSensi/MSK91/Zscore/'], GridVer);
    % input_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
    % 'MSK%2d/'],GridVer);
end
% directory to save figures to
figPath = strcat(input_dir,'FIGS_PCCellv3b_DOC0.25_DOP0/')


VER   = strcat(input_dir,TRdivVer);
catDOC = sprintf('_DOC%0.2g_DOP%0.2g',cscale,pscale); % used to add scale factors to file names
if (Cmodel == off & Omodel == off & Simodel == off & Cellmodel == off)
	fname = strcat(VER,'_P');
elseif (Cmodel == on & Omodel == off & Simodel == off & Cellmodel == off)
	base_name = strcat(VER,'_PC');
	fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == on & Simodel == off & Cellmodel == off)
	base_name = strcat(VER,'_PCO');
	fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == off & Simodel == on & Cellmodel == off)
	base_name = strcat(VER,'_PCSi');
	fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == on & Simodel == on & Cellmodel == off)
	base_name = strcat(VER,'_PCOSi');
	fname = strcat(base_name,catDOC);
elseif (Cmodel == off & Omodel == off & Simodel == off & Cellmodel == on) % cell model does nothing if C model is not on, so this case =Ponly
	base_name = strcat(VER,'_PCell');
	fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == off & Simodel == off & Cellmodel == on)
	base_name = strcat(VER,'_PCCellv3b');
	fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == on & Simodel == off & Cellmodel == on)
	base_name = strcat(VER,'_PCOCell');
	fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == on & Simodel == on & Cellmodel == on)
	base_name = strcat(VER,'_PCOSiCell');
	fname = strcat(base_name,catDOC);
% elseif (Cmodel == on & Omodel == off & Simodel == off)
%     base_name = strcat(VER,'_PCv1');
%     % catDOC = sprintf('_DOC%2.0e',cscale);
%     catDOC = sprintf('_DOC%2.0e_DOP%2.0e',cscale,pscale);
%     fname = strcat(base_name,catDOC);
% elseif (Cmodel == on & Omodel == on & Simodel == off)
%     base_name = strcat(VER,'_PCOv6');
%     catDOC = sprintf('_DOC%2.0e_DOP%2.0e',cscale,pscale);
%     fname = strcat(base_name,catDOC);
% elseif (Cmodel == on & par.Omodel == off & Simodel == on)
%     base_name = strcat(VER,'_PCSi');
%     catDOC = sprintf('_DOC%2.0e_DOP%2.0e',cscale,pscale);
%     fname = strcat(base_name,catDOC);
% elseif (Cmodel == on & Omodel == on & Simodel == on)
%     base_name = strcat(VER,'_PCOSi');
%     catDOC = sprintf('_DOC%2.0e_DOP%2.0e',cscale,pscale);
%     fname = strcat(base_name,catDOC);
end

% load optimal parameters if they exist
fxhat = strcat(fname,'_xhat.mat');
if GridVer == 90
    load transport_v4.mat
    load Sobs_90x180x24.mat     % woa2013 salinity data.
    load tempobs_90x180x24.mat
    load Siobs_90x180x24.mat Siobs
    load po4obs_90x180x24.mat
	load no3obs_90x180x24.mat		% WOA NO3 obs
    grd  = grid;
elseif GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    load Sobs_91x180x24.mat     % woa2013 salinity data.
    load tempobs_91x180x24.mat
    load po4obs_91x180x24.mat
    load Siobs_91x180x24.mat Siobs
    M3d = output.M3d    ;
    grd = output.grid   ;
    TR  = output.TR/spa ;
end
%
load(fxhat)
iwet = find(M3d(:))   ;
nwet = length(iwet)   ;
I    = speye(nwet)    ;
PO4  = po4obs(iwet)   ;
dAt  = grd.DXT3d.*grd.DYT3d ;
dVt  = dAt.*grd.DZT3d ;
dzt  = grd.dzt;

kappa_p    = 1/(720*60^2)   ;
par.taup   = 720*60^2    ; % (s) pic dissolution time-scale
par.tau_TA = 1./par.taup ;
par.tauPIC = (1/0.38)*spd;
par.M3d    = M3d     ;
par.grd    = grd     ;
par.iwet   = iwet    ;
par.nwet   = nwet    ;
par.TRdiv  = -TR     ;
par.dVt    = dVt     ;

tempobs(tempobs(:)<-2.0) = -2.0 ;
par.Temp   = tempobs ;

% ------------------- normalize temperature -------------
for ji = 1:24
    t2d = par.Temp(:,:,ji);
    par.Temp(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
end
vT = par.Temp(iwet) ;
Tz = (vT - min(vT))./(max(vT) - min(vT)) ;
% Tz = zscore(vT)  ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
aveT   = nanmean(Tz3d(:,:,1:2),3) ;
%}

aveT = par.aveT;
Tz = par.Tz;

lat = grd.yt;
lon = grd.xt;

% ----------------make figures---------------------
nfig = 0;
if isfield(xhat,'bP_T')
    nfig = nfig + 1;
    figure(nfig)
    bP   = xhat.bP   ;
    bP_T = xhat.bP_T ;
    bP2D = bP_T*aveT + bP   ;
    pcolor(bP2D) ; colorbar ; shading flat
    title('b4P')
    % saveas(gcf,'Figs91/b4P.png')
	figTitle = 'b4P';
	print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
end

if isfield(xhat,'kP_T')
    nfig = nfig + 1  ;
    figure(nfig)
    kP_T  = xhat.kP_T ;
    kdP   = xhat.kdP  ;

    kP3d = M3d + nan ;
    kP3d(iwet) = (1./(kP_T * Tz + kdP))/spd ;
    pcolor(kP3d(:,:,2));colorbar;shading flat
    % caxis([80 140])
    title('ka4P')
    % saveas(gcf,'Figs91/kappa4P.png')
	figTitle = 'kappa4P';
	print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
end

if par.Cmodel == on
    % DOP remineralization rate constant.
    if isfield(xhat,'bC_T')
        nfig = nfig + 1        ;
        figure(nfig)
        bC   = xhat.bC   ;
        bC_T = xhat.bC_T ;

        bC2D = bC_T*aveT + bC  ;
        pcolor(bC2D); colorbar ; shading flat
        title('b4C')
        % saveas(gcf,'Figs91/b4C.png')
		figTitle = 'b4C';
		print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
    end

    if isfield(xhat,'R_Si')
        nfig = nfig + 1  ;
        figure(nfig)
        par.DSi  = Siobs     ;
        par.R_Si = xhat.R_Si ;
        par.rR   = xhat.rR   ;
        par.opt_R_Si = on    ;
        vout = mkPIC2P(par)  ;
        RR3d = M3d + nan     ;
        RR3d(iwet) = diag(vout.RR) ;
        pcolor(RR3d(:,:,1));colorbar;shading flat
        % caxis([100 600])
        title('rain ratio')
        % saveas(gcf,'Figs91/R_Si.png')
		figTitle = 'RainRatio';
		print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
    end

    if isfield(xhat,'kC_T')
        nfig = nfig + 1  ;
        figure(nfig)
        kC_T = xhat.kC_T ;
        kdC  = xhat.kdC  ;
        kC3d = M3d + nan ;
        kC3d(iwet) = (1./(kC_T * Tz + kdC))/spd;
        pcolor(kC3d(:,:,2));colorbar;shading flat
        % caxis([100 600])
        title('kd4C')
        % saveas(gcf,'Figs91/kappa4C.png')
		figTitle = 'kappa4C';
		print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
    end

    if isfield(xhat,'cc')
        nfig = nfig + 1  ;
        figure(nfig); hold on
        cc   = xhat.cc   ;
        dd   = xhat.dd   ;
		%DIP  = model.DIP(iwet) ;
		DIP = par.po4obs(iwet) ;
        C2P = M3d + nan  ;
        C2P(iwet)  = 1./(cc*DIP + dd) ;

        %pcolor(C2P(:,:,1)); colorbar;shading flat
		imAlpha = ones(size(C2P(:,:,1)));
		imAlpha(isnan(C2P(:,:,1))) =0;
		imagesc(lon,lat,C2P(:,:,1),'AlphaData',imAlpha); hold on
		cb=colorbar;
		colormap(flipud(summer))
		[CC,hh] = contour(lon,lat,C2P(:,:,1),[106 106],'k');
		clabel(CC,hh,'FontName','Times');
		xlabel('Longitude');
		ylabel('Latitude');
		ylabel(cb,'C:P [molC/molP]');
        title('Surface C:P uptake ratio')
		figTitle = 'C2Pratio_po4obs';
		print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
        % saveas(gcf,'Figs91/CP ratio.png')

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

		nfig = nfig + 1  ;
		figure(nfig);
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
		ylabel(['C:P [molC:molP] = 1/(' num2str(cc) '*DIP + ' num2str(dd) ')'])
		title('Zonal Average Biological C:P')
		axis tight; grid on
		ylim([0 400])
		figTitle = 'C2P_lat_avg_po4obs';
		print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
		clear h;


		% --- Pacific vs Atlantic
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
		title('Zonal Average C:P in EZ')
		axis tight; grid on
		ylim([0 400])
		figTitle = 'C2P_lat_avg_basin_po4obs';
		print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
		clear h;

    end
end

if par.Omodel == on
    if isfield(xhat,'O2C_T')
        nfig = nfig + 1    ;
        figure(nfig)
        O2C_T = xhat.O2C_T ;
        rO2C  = xhat.rO2C  ;
        O2C   = M3d + nan  ;
        O2C(iwet) = (O2C_T * Tz + rO2C) ;
        pcolor(O2C(:,:,10)); colorbar; shading flat
        title('O2C consumption ratio')
    end

    if isfield(xhat,'O2P_T')
        par.opt_O2P_T = on      ;
        nfig = nfig + 1         ;
        figure(nfig)
        par.O2P_T  = xhat.O2P_T ;
        par.rO2P   = xhat.rO2P  ;
        vout = mkO2P(par)       ;
        O2P  = M3d + nan        ;
        O2P(iwet) = vout.O2P    ;
        pcolor(nanmean(O2P(:,:,1:2),3)); colorbar; shading flat
        title('O2 production to P ratio')
    end
end

if par.Simodel == on
    nfig  = nfig + 1 ;
    aa    = xhat.aa  ;
    bb    = xhat.bb  ;
    Z     = Siobs(iwet)     ;
    mu    = surface_mean(Z) ;
    Delta = sqrt(surface_mean((Z-mu).^2)) ;

    % standardize the regressor variables
    ZR = 0.5+0.5*tanh((Z-mu)/Delta) ;
    %
    Si2C = (aa*ZR + bb)./Z ;
    figure(nfig)
    pcolor(Si2C);colorbar; shading flat;
end

%if par.Cellmodel ==on
	% CellOut.C2P=data.CellOut.C2P(:,:,1:2);
	% CellOut.N2P=data.CellOut.N2P(:,:,1:2);
	% CellOut.C2N=data.CellOut.C2N(:,:,1:2);
	% CellOut.LimType=data.CellOut.LimType(:,:,1:2);
	% CellOut.r=data.CellOut.r(:,:,1:2);
	% FIGdir = strcat(fname,'_FIGS/');
	% if ~exists(FIGdir)
	% 	mkdir FIGdir
	% end
	% fnamecell =strcat(FIGdir,'PCCell_CellOut_surf.mat');
	% save('fnamecell','CellOut')

%end
