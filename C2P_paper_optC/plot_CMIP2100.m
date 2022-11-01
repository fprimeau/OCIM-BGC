%% notes: imaginary values in the covariance matrix indicate:
% need to run cell model once more after the second deriv so the imaginary part isn't carried through?

clc; clear all; close all
on   = true  ;
off  = false ;

RunVer = 'futureproj/CMIP2100_Cell_CTL_He_PCCell';
%RunVer = 'futureproj/CMIP2100_GM15_CTL_He_PC'
%'optC_Tz_CTL_He_PC_DOC0.25_DOP0'

outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/';
figDir = strcat(outputDir,'futureproj/FIGS_CMIP2100_Cell/');
outPath = figDir;

% load optimal parameter values
fxhat = strcat(outputDir, 'optC_Cell_CTL_He_PCCell_DOC0.25_DOP0_xhat.mat');
%fxhat = strcat(outputDir, 'optC_GM15_CTL_He_PC_DOC0.25_DOP0_xhat.mat');
load(fxhat);

% load model output fields
fname = strcat(outputDir, 'optC_Cell_CTL_He_PCCell_DOC0.25_DOP0.mat');
%fname = strcat(outputDir, 'optC_GM15_CTL_He_PC_DOC0.25_DOP0.mat');
load(fname);
model = data;
clear data

%% Future run
fname_future = strcat(outputDir, RunVer, '.mat');
d = load(fname_future);
model_future = d.model;

% ----------------------------

GridVer  = 91  ;
operator = 'A' ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function
par.dynamicP = off;
par.useProjectionInputs = on;

%-------------set up plots---------------------

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

%-------------load data and set up parameters---------------------
cd ../
SetUp ;
cd mkFIGS/

% get original min and Max temperature for Normalization
vT0 = par.Temp(iwet) ;
vT0min = min(vT0);
vT0max = max(vT0);

% ----overwrite temperature obs from SetUp -----
if GridVer == 90
	%load tempobs_90x180x24.mat
    %load po4obs_90x180x24.mat       % WOA PO4 observation [units: umol/kg]
	%load no3obs_90x180x24.mat		% WOA NO3 obs [units: umol/kg]
elseif GridVer == 91
	load('/DFS-L/DATA/primeau/meganrs/DATA/CMIP5/CMIP5mean_no3_thetao_91x180x24.mat')
end
%iprod = find(M3d(:,:,1:2));

% overwrite po4obs
%par.po4obs = PO4_CMIP5mean;
par.po4proj  = PO4_CMIP5mean  ;

% ----overwrite temperature obs from SetUp -----
tempobs = T_CMIP5mean;
tempobs(tempobs(:)<-2.0) = -2.0 ;

par.Temp_proj    = tempobs ;
%-------------------- normalize temperature --------------------
% for ji = 1:24
%     t2d = par.Temp_proj(:,:,ji);
%     par.Temp_proj(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
% end
vT = par.Temp_proj(iwet) ;
Tz = (vT - vT0min)./(vT0max - vT0min) ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par.Tz_proj     = Tz*1e-8 ;
par.aveT_proj   = nanmean(Tz3d(:,:,1:2),3) ;



%% plot axes
lon = grd.xt;
lat = grd.yt;

%----------- display all parameters -----------------------
fprintf('All Parameter Values: \n')
params = xhat.allparams

% fnames= fieldnames(xhat.allparams);
% for ii=1:length(fnames)
%     paramname=fnames{ii};
%     fprintf('%12s : %.5e \n', paramname, xhat.allparams.(paramname))
% end

fprintf('\n')
fprintf('Optimized Parameter Values: \n')
pnames = fieldnames(xhat);
for ii=1:length(xhat.fx)
    fprintf('%16s : %.5e \n', pnames{ii}, xhat.(pnames{ii}));
end

% ---------- plot C2P -------------
if isfield(xhat,'cc')
	nfig = 1 ;
	figure(nfig); hold on
	cc   = xhat.cc   ;
	dd   = xhat.dd   ;
	%DIP  = model.DIP(iwet) ;
	%DIP = par.po4obs(iwet) ;
	C2P = M3d + nan  ;
	C2Popt = M3d + nan  ;
	C2P(iwet)  = 1./(cc*par.po4proj(iwet) + dd) ;
	C2Popt(iwet)  = 1./(cc*par.po4obs(iwet) + dd) ;

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
	figTitle = 'C2Pratio_CMIP2100';
	print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')

	%% C2P now vs future as a function of latitude
	C2P_latavg = mean(C2P(:,:,1:2),[2 3],'omitnan');
	C2P_latstd = std(C2P(:,:,1:2),0,[2 3],'omitnan');
	ind1 = ~isnan(C2P_latavg);

	C2Popt_latavg = mean(C2Popt(:,:,1:2),[2 3],'omitnan');
	C2Popt_latstd = std(C2Popt(:,:,1:2),0,[2 3],'omitnan');

	nfig = nfig + 1  ;
	figure(nfig);
	set(gcf,'position',[61 614 1331 403],'units','pixels')
	% plot +/-1 standard deviation
	% boundedline(lat,C2P_latavg1,C2P_latstd1,'-b*',lat,C2P_latavg2,C2P_latstd2,'-m*','alpha','nan','gap')
	h(1) = fill([lat(ind1),fliplr(lat(ind1))],[(C2Popt_latavg(ind1)-C2Popt_latstd(ind1))', fliplr((C2Popt_latavg(ind1)+C2Popt_latstd(ind1))')],'k','LineStyle','none'); alpha(0.1); hold on
	h(2) = fill([lat(ind1),fliplr(lat(ind1))],[(C2P_latavg(ind1)-C2P_latstd(ind1))', fliplr((C2P_latavg(ind1)+C2P_latstd(ind1))')],'m','LineStyle','none'); alpha(0.1); hold on

	h(3) = plot(lat,C2Popt_latavg,'-ko'); hold on
	%plot(lat,C2P_latmedian1,'-.b','linewidth',2);
	h(4) = plot(lat,C2P_latavg,'-mo');
	%plot(lat,C2P_latmedian2,'-.m','linewidth',2);
	legend(h([3 4]),'optimal','CMIP2100');
	xlabel('Latitude')
	ylabel(['C:P [molC:molP] = 1/(' num2str(cc) '*[PO4] + ' num2str(dd) ')'])
	title('Zonal Average Biological C:P')
	axis tight; grid on
	%ylim([0 300])
	figTitle = 'C2P_lat_avg_CMIP2100vopt';
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
	ylim([0 300])
	figTitle = 'C2P_lat_avg_basin_CMIP2100';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
	clear h;

end


% %------- Cell Output Stats -----------
if par.Cellmodel == on
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

	iprod = find(M3d(:,:,1:2));
	uLimTypes = unique(model.CellOut.LimType(iprod));
	fprintf('MODERN CONDITIONS: \n')
	fprintf('# of N limited points : %d \n',length(find(LimType(iprod) ==0)));
	fprintf('# of P limited points : %d \n',length(find(LimType(iprod) ==1)));
	fprintf('# of Colimited2 points: %d \n',length(find(LimType(iprod) ==2)));
	fprintf('# of Colimited3 points: %d \n\n',length(find(LimType(iprod) ==3)));

	fprintf('range of C:P values is %6.1f to %6.1f \n',min(C2P(iprod)),max(C2P(iprod)))

	fprintf('mean of C:P is %6.1f \n',mean(C2P(iprod)))
	fprintf('std of C:P  is %6.1f \n\n',std(C2P(iprod)))

	fprintf('range of radius is      %6.2f to %6.2f um \n',min(radius(iprod)),max(radius(iprod)))
	fprintf('mean of radius is %6.1f \n',mean(radius(iprod)))
	fprintf('std of radius  is %6.1f \n\n',std(radius(iprod)))

	fprintf('range of growth rate is %6.3f to %6.3f hr^-1 \n',min(mu(iprod)),max(mu(iprod)))
	fprintf('range of growth rate is %6.3f to %6.3f day^-1 \n\n',min(mu(iprod)*24),max(mu(iprod)*24))


	%% ------- FUTURE ----------

	C2P_future     = model_future.CellOut.C2P;
	N2P_future     = model_future.CellOut.N2P;
	C2N_future     = model_future.CellOut.C2N;
	radius_future  = model_future.CellOut.r;
	LimType_future = model_future.CellOut.LimType;
	mu_future      = model_future.CellOut.mu;

	% set land points to NaN
	C2P_future(find(C2P_future==0)) = NaN;
	C2N_future(find(C2N_future==0)) = NaN;
	N2P_future(find(N2P_future==0)) = NaN;
	%radius(find(radius==0)) = NaN;

	fprintf('UNDER FUTURE CONDITIONS: CMIP5 2100 projection: \n')
	fprintf('# of N limited points : %d \n',length(find(LimType_future(iprod) ==0)));
	fprintf('# of P limited points : %d \n',length(find(LimType_future(iprod) ==1)));
	fprintf('# of Colimited2 points: %d \n',length(find(LimType_future(iprod) ==2)));
	fprintf('# of Colimited3 points: %d \n\n',length(find(LimType_future(iprod) ==3)));

	fprintf('range of C:P values is %6.1f to %6.1f \n',min(C2P_future(iprod)),max(C2P_future(iprod)))

	fprintf('mean of C:P is %6.1f \n',mean(C2P_future(iprod)))
	fprintf('std of C:P  is %6.1f \n\n',std(C2P_future(iprod)))

	fprintf('range of radius is      %6.2f to %6.2f um \n',min(radius_future(iprod)),max(radius_future(iprod)))
	fprintf('mean of radius is %6.1f \n',mean(radius_future(iprod)))
	fprintf('std of radius  is %6.1f \n\n',std(radius_future(iprod)))

	fprintf('range of growth rate is %6.3f to %6.3f hr^-1 \n',min(mu_future(iprod)),max(mu_future(iprod)))
	fprintf('range of growth rate is %6.3f to %6.3f day^-1 \n\n',min(mu_future(iprod)*24),max(mu_future(iprod)*24))
end
%

%% --------------------- MAKE PLOTS ----------------------------------







%%------------- Calculate C export ------------------
% ----------------------------------------------
nn= 2; %nn = par.nzo ; %number fo verticle boxes in euphotic zone / export depth

% -------------- C:P uptake ratio --------------------------------
W = d0(dVt(iwet)) ;
dVtwet = M3d*nan;
dVtwet(iwet) = dVt(iwet);
Wiprod = dVtwet(:,:,1:nn)/nansum(dVtwet(:,:,1:nn),'all');

C2P3D_now = M3d + nan ;
C2P3D = M3d + nan ;
if par.Cellmodel==on
	C2P3D_now(iwet) = model.CellOut.C2P(iwet); %zero beneath the surface  layers
	C2P3D(iwet) = model_future.CellOut.C2P(iwet); %zero beneath the surface  layers
elseif par.Cmodel ==on
	C2P3D(iwet) = 1./(params.cc*par.po4proj(iwet) + params.dd) ;  % DIP or PO4?
	%C2P3D(iwet) = 1./(7e-4*PO4 + 5e-3); % WL
	%C2P3D(iwet) = 1000./(6.6*PO4 + 5.3); %Qian
else
	C2P3D = M3d +nan;
	C2P3D(iwet) = 106; 		% redfield C:P
end

% DIP assimilation
LAM        = 0*M3d;
LAM(:,:,1) = (par.npp1.^params.beta).*par.Lambda(:,:,1);
LAM(:,:,2) = (par.npp2.^params.beta).*par.Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];

%--------------- calculate primary production --------------------
fprintf('----------- Future Cexp ------------- \n')
G        = M3d*0        ;
G(iwet)  = params.alpha*L*model_future.DIP(iwet)  ; % primary production [unit: mmol P/m^3/s]

% inegG = find(G<0); % should negative production values be removed?
% G(inegG)=nan;

Int_PNPP = G(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31;
Int_CNPP = G(:,:,1:nn).*grd.DZT3d(:,:,1:nn).*C2P3D(:,:,1:nn)*12;

PNPP = Int_PNPP*spa*1e-3 ; % convert to g P/m^2/yr
CNPP = Int_CNPP*spa*1e-3 ; % convert production from mg C/m^3/s to gC/m^2/year;
tem_PNPP = PNPP.*dAt(:,:,1:nn)*1e-15 ;
tem_CNPP = CNPP.*dAt(:,:,1:nn)*1e-15 ;
Sum_CNPP = nansum(tem_CNPP(:))    ;
fprintf('Model NPP (P) is %3.3e Pg P/yr \n',nansum(tem_PNPP(:))) ; %Pg/yr
fprintf('Model NPP is %3.3e Pg C/yr \n\n',Sum_CNPP) ; %Pg/yr

clear tem_PNPP tem_CNPP

%----- C2P --------------
C2Pavg = nansum(C2P3D(:,:,1:nn).*Wiprod,'all');
fprintf('Average C:P uptake in Euphotic Layers is %4.2f \n',C2Pavg);

% prod weighted C:P
WeightNPP = CNPP.*dAt(:,:,1:nn)/nansum(CNPP.*dAt(:,:,1:nn),'all'); %gC/yr
C2Pavg = nansum(C2P3D(:,:,1:nn).*WeightNPP,'all');
fprintf('Average C:P uptake in Euphotic Layers (NPP weighted) is %4.2f \n\n',C2Pavg);


% ------- P export -----
%POP export: integrte POP beneath the euphotic zone
POPexp = nansum(params.kappa_p*model_future.POP(:,:,3:end).*dVt(:,:,3:end),3)*31*spd;
Sum_POPexp =nansum(POPexp(:))*365*1e-18;
fprintf('Model POP export is %3.3e Pg P /yr  (beneath EZ) \n',Sum_POPexp);

%integrated DOP remineralization below the euphotic zone. (should equal the TOP export calculated by the adjoint method)
DOPexpint = params.kdP*model_future.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
tem_DOPexpint = nansum(DOPexpint.*dAt(:,:,3:end),3);
Sum_DOPexpint =nansum(tem_DOPexpint(:))*365*1e-18;
fprintf('Model TOP export (integrated DOP below the Euphotic zone) is %3.3f Pg P /yr \n\n',Sum_DOPexpint);
DOPexpint = sum(DOPexpint,3,'omitnan');


%----- C export ---------------
%POC export: integrate POC beneath the euphotic zone ----------
%POCexp = nansum(par.kappa_p*model.POC(:,:,3:end).*dVt(:,:,3:end),3)*12*spd;  % [1/s]*[mmol/m^3]*[m^3]*[mg/mmol]*[s/day] =  [mg/day]
POCexp = params.kappa_p*model_future.POC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
tem_POCexp = nansum(POCexp.*dAt(:,:,3:end),3);
Sum_POCexp =nansum(tem_POCexp(:))*365*1e-18;
fprintf('Model POC export is %3.3e Pg C /yr  (beneath EZ) \n',Sum_POCexp);
POCexp = sum(POCexp,3,'omitnan');

%integrated DOC remineralization below the euphotic zone. (should equal the TOC export calculated by the adjoint method)
DOCexpint = params.kdC*model_future.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
tem_DOCexpint = nansum(DOCexpint.*dAt(:,:,3:end),3);
Sum_DOCexpint =nansum(tem_DOCexpint(:))*365*1e-18;
fprintf('Model TOC export (integrated DOC below the Euphotic zone) is %7.5f Pg C /yr \n\n',Sum_DOCexpint);
DOCexpint = sum(DOCexpint,3,'omitnan');

%% DOCexp int future vs Now
DOCexpint_now = params.kdC*model.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
DOCexpint_now = sum(DOCexpint_now,3,'omitnan');

%
TOCexpdiff = (DOCexpint - DOCexpint_now)./12/1000*365;
[TOCmin, TOCmax] = bounds(TOCexpdiff,'all');

figure
imAlpha = M3d(:,:,1); %~isnan(TOCexpdiff);
imAlpha(isnan(TOCexpdiff)) = 0;
imagesc(grd.xt,grd.yt,TOCexpdiff,'AlphaData',imAlpha); hold on
cb = colorbar; colormap(gca,darkb2r(TOCmin, TOCmax));
%colormap(gca,darkb2r(-1.4, 7));
%[CC,hh] = contour(grd.xt,grd.yt,TOCexp/12/1000*365,[0:2:8],'k');
%clabel(CC,hh,'FontName','Times','FontSize',6);
axis xy;
title('Future - Modern differece in Total Organic Carbon Export')
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'\Delta TOC export [mol C/m^2/yr]');
figTitle = ['CMIP2100_deltaTOCexp'];
set(gca, 'color',[0.8 0.8 0.8])
set(gcf,'InvertHardcopy','off','Color',[1 1 1])
exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']) %,'BackgroundColor','none','ContentType','vector')
%print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')


% -------Cell model figures ------------------
if par.Cellmodel == on
	dim = [0.1199 0.695 0.1 0.2];
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

	%% difference between current and +5C MAPS
	C2P_latavg_now = mean(C2P(:,:,1:2),[2 3],'omitnan');
	C2P_latstd_now = std(C2P(:,:,1:2),0,[2 3],'omitnan');

	C2P_latavg_now1 = mean(C2P(:,:,1),[2],'omitnan');
	C2P_latavg_now2 = mean(C2P(:,:,2),[2],'omitnan');

	C2P_latavg_future = mean(C2P_future(:,:,1:2),[2 3],'omitnan');
	C2P_latstd_future = std(C2P_future(:,:,1:2),0,[2 3],'omitnan');

	C2P_latavg_future1 = mean(C2P_future(:,:,1),[2],'omitnan');
	C2P_latavg_future2 = mean(C2P_future(:,:,2),[2],'omitnan');

	ind1 = ~isnan(C2P_latavg_now);
	ind2 = ~isnan(C2P_latavg_future);

	figure;
	set(gcf,'position',[61 614 1331 403],'units','pixels')
	h(1) = fill([lat(ind1),fliplr(lat(ind1))],[(C2P_latavg_now(ind1)-C2P_latstd_now(ind1))', fliplr((C2P_latavg_now(ind1)+C2P_latstd_now(ind1))')],'k','LineStyle','none'); alpha(0.1); hold on

	h(2) = fill([lat(ind2),fliplr(lat(ind2))],[(C2P_latavg_future(ind2)-C2P_latstd_future(ind2))', fliplr((C2P_latavg_future(ind2)+C2P_latstd_future(ind2))')],'m','LineStyle','none'); alpha(0.1);


	h(3) = plot(lat,C2P_latavg_now,'-o','Color','k'); hold on
	%plot(lat,C2P_latmedian1,'-.b','linewidth',2);
	h(4) = plot(lat,C2P_latavg_future,'-mo');

	%h(5) = plot(lat,C2P_latavg_now1,'b--'); hold on
	%h(6) = plot(lat,C2P_latavg_now2,'b-.'); hold on
	%h(7) = plot(lat,C2P_latavg_future1,'--','Color',colors.maroon); hold on
	%h(8) = plot(lat,C2P_latavg_future2,'-.','Color',colors.maroon); hold on

	legend(h([3 4]),'Modern','2100');
	xlabel('Latitude')
	ylabel('C:P [molC:molP]')
	title('Zonal Average Cellular C:P')
	axis tight; grid on
	%ylim([0 400])
	figTitle = 'C2P_lat_avg_nowvCMIP2100';
	print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
	clear h;



	%% make difference map
	C2Pdiff = C2P_future(:,:,1:2)-C2P(:,:,1:2);

	figure; hold on;

	imAlpha = ones(size(C2Pdiff(:,:,1)));
	imAlpha(isnan(C2Pdiff(:,:,1))) =0;
	imagesc(lon,lat,C2P(:,:,1),'AlphaData',imAlpha)
	cb=colorbar;
	cmocean('balance','pivot',0);
	%cmap = cmocean('-curl','pivot',0);
	%colormap(cmap)
	%[CC,hh] = contour(lon,lat,C2P(:,:,1),[106 106],'k');
	%clabel(CC,hh,'FontName','Times');
	tstr = sprintf('Future Scenario: CMIP5 2100 projection \nCell Model C:P Change from Modern: Surface');
	tstr = {'Future Scenario: CMIP5 2100 projection', 'Cell Model C:P Change from Modern: Surface'};
	title(tstr,'Fontsize',14);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'\Delta C:P (future - modern) [molC/molP]');

	annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
	axis tight; grid off

	figTitle = 'deltaC2Psurface_CMIP2100';
	print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')


	%% Future Allocations
	L = model_future.CellOut.L;
	E = model_future.CellOut.E;
	S = 2.*model_future.CellOut.A;
	Zlevs = [0:0.1:1];


	%Fallocation = figure('units','normalized','outerposition',[0 0 1 1])
	Fallocation = figure('outerposition', [1 1 1600 1000]);
	tl = tiledlayout(2,3);

	% L surface
	nexttile
	imAlpha = M3d(:,:,1);
	imagesc(lon,lat,L(:,:,1),'AlphaData',imAlpha);
	hold on;
	cb=colorbar;
	caxis([Zlevs(1) Zlevs(end)]);
	cmocean('matter',2*(length(Zlevs)-1));
	[CC,hh] = contour(lon,lat,L(:,:,1),[0:0.1:0.7],'k');
	clabel(CC,hh,'FontName','Times');
	title('Allocation to L: Surface','Fontsize',14);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'L [volume fraction of cell]');
	axis xy; axis tight;
	%annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');

	% E surface
	nexttile;
	imAlpha = M3d(:,:,1);
	imagesc(lon,lat,E(:,:,1),'AlphaData',imAlpha)
	hold on;
	cb=colorbar;
	caxis([Zlevs(1) Zlevs(end)]);
	cmocean('matter',2*(length(Zlevs)-1));
	[CC,hh] = contour(lon,lat,E(:,:,1),[0:0.1:0.7],'k');
	clabel(CC,hh,'FontName','Times');
	title('Allocation to E: Surface','Fontsize',14);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'E [volume fraction of cell]');
	axis tight; axis xy

	% S surface
	nexttile;
	imAlpha = ones(size(S(:,:,1)));
	imAlpha(isnan(S(:,:,1))) =0;
	imagesc(lon,lat,S(:,:,1),'AlphaData',imAlpha)
	hold on
	cb=colorbar;
	caxis([Zlevs(1) Zlevs(end)]);
	cmocean('matter',2*(length(Zlevs)-1));
	[CC,hh] = contour(lon,lat,S(:,:,1),[0:0.1:1],'k');
	clabel(CC,hh,'FontName','Times','Color','k'); %[0.8 0.8 0.8]
	title('Allocation to A+M: Surface','Fontsize',14);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'Structure [volume fraction of cell]');
	axis xy; axis tight

	% L lower EZ
	nexttile;
	imAlpha = M3d(:,:,2);
	imagesc(lon,lat,L(:,:,2),'AlphaData',imAlpha)
	hold on;
	cb=colorbar;
	caxis([Zlevs(1) Zlevs(end)]);
	cmocean('matter',2*(length(Zlevs)-1));
	[CC,hh] = contour(lon,lat,L(:,:,2),[0:0.1:0.7],'k');
	clabel(CC,hh,'FontName','Times');
	title('Allocation to L: Lower EZ','Fontsize',14);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'L [volume fraction of cell]');
	axis xy; axis tight;

	% E lower EZ
	nexttile;
	imAlpha = M3d(:,:,2);
	imagesc(lon,lat,E(:,:,2),'AlphaData',imAlpha)
	hold on;
	cb=colorbar;
	caxis([Zlevs(1) Zlevs(end)]);
	cmocean('matter',2*(length(Zlevs)-1));
	[CC,hh] = contour(lon,lat,E(:,:,2),[0:0.1:0.7],'k');
	clabel(CC,hh,'FontName','Times');
	title('Allocation to E: Lower EZ','Fontsize',14);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'E [volume fraction of cell]');
	axis xy; axis tight;

	% S lower EZ
	nexttile;
	imAlpha = ones(size(S(:,:,2)));
	imAlpha(isnan(S(:,:,2))) =0;
	imagesc(lon,lat,S(:,:,2),'AlphaData',imAlpha)
	hold on;
	cb=colorbar;
	caxis([Zlevs(1) Zlevs(end)]);
	cmocean('matter',2*(length(Zlevs)-1));
	[CC,hh] = contour(lon,lat,S(:,:,2),[0:0.1:1],'k');
	clabel(CC,hh,'FontName','Times','Color','k'); % [0.8 0.8 0.8]
	title('Allocation to A+M: Lower EZ','Fontsize',14);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'Structure [volume fraction of cell]');
	axis xy; axis tight

	title(tl,'Future Scenario: CMIP5 2100 projection')
	figTitle = 'Allocations_CMIP2100';
	%print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
	exportgraphics(gcf,[figDir 'FIG_' figTitle '.png']);
end


fprintf('-------DONE!------------ \n')


%{
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
print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
clear h;


% --- Limitation Type ----
Lim_cmap = [1, 0, 0; ...
    0, 0, 1; ...
    0, 1, 1; ...
    0, 1, 1];
figure; hold on

imAlpha = ones(size(LimType(:,:,1)));
imAlpha(isnan(LimType(:,:,1))) =0;
colormap(Lim_cmap)
%shifted so smallest value is 1, so can use direct mapping
image(lon,lat,LimType(:,:,1)+1,'AlphaData',imAlpha)
%colormap(Lim_cmap);
cb=colorbar('Ticks',[1,2,3,4],'TickLabels',{'N-Lim','P-Lim','Co-Lim','Co-Lim-alt'});
ylabel(cb,'Limitation Type');
axis tight
title('Phytoplankton Nutrient Limitation Type: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
grid off

figTitle = 'LimType_surf';
print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')


%%------ P storage ---
PStorpct = model.CellOut.PStor./model.CellOut.QP;
PLippct = model.CellOut.PLip./model.CellOut.QP;

figure
t1 = tiledlayout('flow');
nexttile
histogram(PStorpct,'Normalization','probability')
%labeltxt = sprintf('molar fraction of Pstorage / PQuota \n given fStorage = %5.3f ; rCutoff = %5.3f ; PStor-scale = %5.3f', params.fStorage,params.PStor_rCutoff, params.PStor_scale);
title('fraction of total P quota in Storage molecules')
%xlabel('molar fraction of Pstorage / PQuota')
%xlabel(labeltxt)
ylabel('probability')

nexttile
histogram(PStorpct,'Normalization','cdf')
labeltxt = sprintf('molar fraction of Pstorage / PQuota \n given fStorage = %5.3f ; rCutoff = %5.3f ; PStor-scale = %5.3f', params.fStorage,params.PStor_rCutoff, params.PStor_scale);
title('fraction of total P quota in Storage molecules')
%xlabel('molar fraction of Pstorage / PQuota')
xlabel(labeltxt)
ylabel('cdf')

figTitle= sprintf('PStor_hist');

print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')

%}
