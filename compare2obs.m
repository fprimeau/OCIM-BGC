clc; clear all; close all
spd  = 24*60^2 ; spa  = 365*spd ;
on = true; off = false;
%
set(groot,'defaultTextInterpreter','Tex');

%RunVer = 'testNPP_Tv4_PCCella1e-8b1e-3_DOC0.25_DOP0';
%RunVer = 'testPobs_CTL_He_PCCella1b1_DOC0.25_DOP0';
%RunVer = 'Tv4_PCCellv9c_fixedQ10_DOC0.25_DOP0'
%RunVer = 'Tv4_PC_DOC0.25_DOP0v8_onlyC2P'
%RunVer ='constC2P_Tv4_PCv8_DOC0.25_DOP0'
RunVer = 'ESS225_Tv4_PCCell_DOC0.25_DOP0'
%RunVer = 'Tv4_PCv9_DOC0.25_DOP0'
%RunVer = 'testPobs_Tv4_PCCella8e-4b6e-1_DOC0.25_DOP0';

%model output directory
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/';
%figPath = strcat(outputDir,'FIGS_testPobs_PCCell/a1b1_');
%figPath = strcat(outputDir,'FIGS_testNPP_PCCell/a1e-8b1e-3_');
%figPath = strcat(outputDir,'FIGS_PCv9_DOC0.25_DOP0/');
%figPath = strcat(outputDir,'FIGS_PCCellv9_DOC0.25_DOP0/v9c_fixedQ10_');
%figPath = strcat(outputDir,'FIGS_PC_DOC0.25_DOP0v8/onlyC2P/');
%figPath = strcat(outputDir,'FIGS_constC2P_PCv8_DOC0.25_DOP0/');
figPath = strcat(outputDir,'FIGS_ESS225/');


% load model output fields
fname = strcat(outputDir, RunVer, '.mat');
load(fname);
model = data;

% load optimal parameter values
fxhat = strcat(outputDir, RunVer,'_xhat.mat');
load(fxhat);

% PCCellv8 = load(strcat(outputDir,'Tv4_PCCellv8_DOC0.25_DOP0_xhat.mat'));
% finames = fieldnames(PCCellv8.xhat);
% for ii = 1:11
% 	xhat.(finames{ii}) = PCCellv8.xhat.(finames{ii});
% end

GridVer  = 90  ;
operator = 'A' ;
par.Pmodel  = on ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
%par.cscale  = 0.25 ; % factor to weigh DOC in the objective function
par.cscale = 0.25;
par.dynamicP = off;

%-------------load data and set up parameters---------------------
SetUp ;

%--------------------- prepare parameters ------------------
xhat
% if par.LoadOpt == on
%     % load optimal parameters from a file or set them to default values
%     par = SetPar(par) ;
%     % pack parameters into an array, assign them corresponding indices.
%     par = PackPar(par) ;
% 	x0    = par.p0 ;
% 	iter = 0;
% 	PrintPara(x0, par) ;
% end

%----------------surface indices---------------------------
iprod = find(par.M3d(:,:,1:2));
iprod1 = find(par.M3d(:,:,1)); % surface layer only
iprod2 = find(par.M3d(:,:,2)); % second EZ depth layer

%%% PLOT ONLY SURFACE VALUES
%iwet = iprod;
%figPath = strcat(input_dir,'FIGS_PCCellv3b_DOC0.25_DOP0/surf')

%----------------NPP---------------


%------------------ compare P model ---------------------------------
%------ DOP ---------
nfig = 0;
if (par.Pmodel == on)
	DOP = model.DOP ;
	DIP = model.DIP ;

    %par.dopraw = DOPobs - 0.03 ; % less refractory DOP
	dopraw = par.dopraw;
    idop   = find(dopraw(iwet) > 0 & DOP(iwet)>0) ;
    iDOP_ATL = find(dopraw(iwet)>0 & ATL(iwet)>0) ;
    iDOP_PAC = find(dopraw(iwet)>0 & PAC(iwet)>0) ;
    iDOP_IND = find(dopraw(iwet)>0 & IND(iwet)>0) ;
    iDOP_ARC = find(dopraw(iwet)>0 & ARC(iwet)>0) ;
    iDOP_MED = find(dopraw(iwet)>0 & MED(iwet)>0) ;
    fprintf('R^2 for DOP is %3.4f \n',rsquare(dopraw(iwet(idop)),DOP(iwet(idop))))

    nfig = nfig + 1;
    figure(nfig)
    plot(dopraw(iwet(iDOP_ATL)), DOP(iwet(iDOP_ATL)),'ro')
    hold on
    plot(dopraw(iwet(iDOP_PAC)), DOP(iwet(iDOP_PAC)),'ks')
    hold on
    plot(dopraw(iwet(iDOP_IND)), DOP(iwet(iDOP_IND)),'b^')
    hold on
    plot(dopraw(iwet(iDOP_ARC)), DOP(iwet(iDOP_ARC)),'g*')
    hold on
    plot(dopraw(iwet(iDOP_MED)), DOP(iwet(iDOP_MED)),'c>')
    hold on
	legend('ATL','PAC','IND','ARC','Location','northwest')
	hold on
    plot([0 0.75],[0 0.75],'r-','linewidth',3)
    xlim([0 0.75])
    ylim([0 0.75])

	xlabel('Observed DOP (mmol/m^3)');
    ylabel('Model DOP (mmol/m^3)');
	figTitle = 'DOPcompare2obs';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

%------ DIP ----------------
    nfig = nfig+1;
    figure(nfig)
    ipo4 = find(DIP(iwet) > 0 & po4raw(iwet) > 0.02);
    O = po4raw(iwet(ipo4));
    M = DIP(iwet(ipo4));
    fprintf('R^2 for DIP is %3.4f \n',rsquare(O,M))
	R2string = sprintf('R^2 = %3.4f   \n',rsquare(O,M));

    OvsM = [O,M];
    W = (dVt(iwet(ipo4))./sum(dVt(iwet(ipo4))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,100,[0 0],[4 4],W);
    cr = 5:5:95;
    dx = X(3,5)-X(3,4);
    dy = Y(4,2)-Y(3,2);
    [q,ii] = sort(density(:)*dx*dy,'descend');
    D = density;
    D(ii) = cumsum(q);
    subplot('position',[0.2 0.2 0.6 0.6])
    contourf(X,Y,100*(1-D),cr); hold on
    contour(X,Y,100*(1-D),cr);

    caxis([5 95])
    %set(gca,'FontSize',16);
    grid on
    axis square
    xlabel('Observed DIP (mmol/m^3)');  % GLODAP units = umol/kg
    ylabel('Model DIP (mmol/m^3)');
    title('Model V.S. Observed DIP')
    plot([0 4],[0 4],'r--','linewidth',2);
	text(4, 0, R2string,'VerticalAlignment','bottom','HorizontalAlignment','right','Units','data');

    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    % exportfig(gcf,'DIP_MvsO','fontmode','fixed','fontsize',12, ...
              % 'color','rgb','renderer','painters')
	figTitle = 'DIPcompare2obs';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

	%%%%%%% surface only %%%%%%%%%%%%%
	nfig = nfig+1;
    figure(nfig)
	iwetsurf = find(par.M3d(:,:,1:2)) ;
    %ipo4 = find(DIP(iwetsurf) > 0 & po4raw(iwetsurf) > 0.02);
	ipo4 = find(DIP(iwetsurf) & ~isnan(po4raw(iwetsurf)) );
	%O = po4raw(iwetsurf(ipo4));
    O = po4raw(iwetsurf(ipo4));
    M = DIP(iwetsurf(ipo4));
    fprintf('R^2 for surface DIP is %3.4f \n',rsquare(O,M))
	R2string = sprintf('R^2 = %3.4f   \n',rsquare(O,M)); %add spaces before if putting at top left of plot.

    OvsM = [O,M];
    W = (dVt(iwetsurf(ipo4))./sum(dVt(iwetsurf(ipo4))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,100,[-0.05 -0.05],[4 4],W);
    cr = 5:5:95;
    dx = X(3,5)-X(3,4);
    dy = Y(4,2)-Y(3,2);
    [q,ii] = sort(density(:)*dx*dy,'descend');
    D = density;
    D(ii) = cumsum(q);
    subplot('position',[0.2 0.2 0.6 0.6])
    contourf(X,Y,100*(1-D),cr); hold on
    contour(X,Y,100*(1-D),cr);
	axis equal

    caxis([5 95])
    %set(gca,'FontSize',16);
    grid on
    axis square
    xlabel('GLODAP Observed DIP (mmol/m^3)');  % GLODAP & WOA units = umol/kg
    ylabel('Model DIP (mmol/m^3)');
    title('Model vs Observed DIP (Euphotic Zone)')
    plot([0 4],[0 4],'r--','linewidth',2);
	%text(0, 4, R2string,'VerticalAlignment','top','Units','data');
	text(4, 0, R2string,'VerticalAlignment','bottom','HorizontalAlignment','right','Units','data');

    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    % exportfig(gcf,'DIP_MvsO','fontmode','fixed','fontsize',12, ...
              % 'color','rgb','renderer','painters')
	figTitle = 'DIPsurfcompare2obs';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')


	% simple scatter plot of surface DIP compare to obs
	minDIP = min(OvsM,[],'all');
	maxDIP = max(OvsM,[],'all');
	nfig = nfig+1;
	figure(nfig)
	%plot(O,M,'.k'); hold on
	plot(O(ATL(iwetsurf(ipo4))>0), M(ATL(iwetsurf(ipo4))>0),'.r'); hold on
	plot(O(PAC(iwetsurf(ipo4))>0), M(PAC(iwetsurf(ipo4))>0),'.b')
	plot(O(IND(iwetsurf(ipo4))>0), M(IND(iwetsurf(ipo4))>0),'.g')
	plot(O(ARC(iwetsurf(ipo4))>0), M(ARC(iwetsurf(ipo4))>0),'.c')
	plot(O(MED(iwetsurf(ipo4))>0), M(MED(iwetsurf(ipo4))>0),'.m')
	plot([minDIP maxDIP], [minDIP maxDIP],'k--','LineWidth',2);
	axis equal;
	axis([minDIP maxDIP minDIP maxDIP]);
	axis square; grid on;
	legend('ATL','PAC','IND','ARC','MED','1:1 line','Location','northwest')
	text(maxDIP, minDIP, R2string,'VerticalAlignment','bottom','HorizontalAlignment','right','Units','data');
	%text(minDIP, maxDIP, R2string,'VerticalAlignment','top','Units','data');
	xlabel('GLODAP Observed DIP (mmol/m^3)');  % GLODAP & WOA units = umol/kg
    ylabel('Model DIP (mmol/m^3)');
    title('DIP: Surface 2 layers')
	figTitle = 'DIPsurfcompare2obs_scatter';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

% ------------ SURFACE MAPS ------------
	% surface map of difference from obs
	lon = grd.xt;
	lat = grd.yt;
	figure; hold on;
	imAlpha = M3d(:,:,1);
	MOdiff = DIP(:,:,1)-par.po4obs(:,:,1);
	imagesc(lon,lat,MOdiff,'AlphaData',imAlpha)
	cb=colorbar;
	cmocean('balance','pivot',0);
	title('Model DIP - WOA obs: Surface','Fontsize',18);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'DIP [mmol/m^3]');
	axis tight;
	figTitle = 'DIP_MvsO_surface_map';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

	figure; hold on;
	imAlpha = M3d(:,:,1);
	MOrelerr = (DIP(:,:,1)-par.po4obs(:,:,1))./(par.po4obs(:,:,1));
	imagesc(lon,lat,MOrelerr,'AlphaData',imAlpha)
	cb=colorbar;
	cmocean('balance','pivot',0);
	[CC,hh] = contour(lon,lat,MOrelerr,[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1],'k');
	clabel(CC,hh,'FontName','Times','FontSize',8);
	title('DIP Model relative error to WOA obs: Surface','Fontsize',16);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'DIP relative error (model DIP - WOA obs)/(WOA obs)');
	axis tight;
	figTitle = 'DIP_MvsO_relerr_surface_map';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

end

% -----------------------------------------------------
if (par.Cmodel == on)
	DIC = model.DIC - par.dicant;

    nfig = nfig + 1;
    figure(nfig)
    DICobs = par.dicraw ;
    iDIC = find(DICobs(iwet)>0);
    %
    O = DICobs(iwet(iDIC));
    % already including anthropogenic CO2
    % M = DIC(iwet(iDIC)) ;
    % not include anthropogenic CO2
    M = DIC(iwet(iDIC))+par.dicant(iwet(iDIC));
    fprintf('R^2 for DIC is %3.4f \n',rsquare(O,M))
	R2string = sprintf('R^2 = %3.4f   \n',rsquare(O,M));
    %
    OvsM = [O, M];
    W = (dVt(iwet(iDIC))./sum(dVt(iwet(iDIC))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,100,[2000 2000],[2500 2500],W);
    cr = 5:5:95;
    dx = X(3,5)-X(3,4);
    dy = Y(4,2)-Y(3,2);
    [q,ii] = sort(density(:)*dx*dy,'descend');
    D  = density;
    D(ii) = cumsum(q);
    subplot('position',[0.2 0.2 0.6 0.6])
    contourf(X,Y,100*(1-D),cr); hold on
    contour(X,Y,100*(1-D),cr);

    caxis([5 95])
    %set(gca,'FontSize',16);
    grid on
    axis square
    xlabel('Observed DIC (mmol/m^3)');
    ylabel('Model DIC (mmol/m^3)');
    title('Model V.S. Observed DIC')
    plot([2000 2500],[2000 2500],'r--','linewidth',2); hold on;
	text(2500, 2000, R2string,'VerticalAlignment','bottom','HorizontalAlignment','right','Units','data');

    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 1.5],cr,[cr(:),cr(:)],cr); hold on
    contour([1 1.5],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    % exportfig(gcf,'DIC_MvsO','fontmode','fixed','fontsize',12,...
              % 'color','rgb','renderer','painters')
	figTitle = 'DICcompare2obs';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

	%----surface Only DIC
	dim = [0.1199 0.695 0.1 0.2];

		nfig = nfig + 1;
	    figure(nfig)
	    DIC = model.DIC - par.dicant ;
	    DICobs = par.dicraw ;
	    iDIC = find(DICobs(iwetsurf)>0);
	    %
	    O = DICobs(iwetsurf(iDIC));
	    % already including anthropogenic CO2
	    % M = DIC(iwet(iDIC)) ;
	    % not include anthropogenic CO2
	    M = DIC(iwetsurf(iDIC))+par.dicant(iwetsurf(iDIC));
	    fprintf('R^2 for surface DIC is %3.4f \n',rsquare(O,M))
    	R2string = sprintf('R^2 = %3.4f   \n',rsquare(O,M));

	    OvsM = [O, M];
	    W = (dVt(iwet(iDIC))./sum(dVt(iwet(iDIC))));
	    [bandwidth,density,X,Y] = mykde2d(OvsM,100,[2000 2000],[2500 2500],W);
	    cr = 5:5:95;
	    dx = X(3,5)-X(3,4);
	    dy = Y(4,2)-Y(3,2);
	    [q,ii] = sort(density(:)*dx*dy,'descend');
	    D  = density;
	    D(ii) = cumsum(q);
	    subplot('position',[0.2 0.2 0.6 0.6])
	    contourf(X,Y,100*(1-D),cr); hold on
	    contour(X,Y,100*(1-D),cr);

	    caxis([5 95])
	    %set(gca,'FontSize',16);
	    grid on
	    axis square
	    xlabel('Observed DIC (mmol/m^3)');
	    ylabel('Model DIC (mmol/m^3)');
	    title('Model V.S. Observed DIC (Euphotic Zone)')
	    plot([2000 2500],[2000 2500],'r--','linewidth',2);
		text(2500, 2000, R2string,'VerticalAlignment','bottom','HorizontalAlignment','right','Units','data');

	    subplot('position',[0.82 0.2 0.05 0.6]);
	    contourf([1 1.5],cr,[cr(:),cr(:)],cr); hold on
	    contour([1 1.5],cr,[cr(:),cr(:)],cr);
	    hold off
	    %set(gca,'FontSize',14);
	    set(gca,'XTickLabel',[]);
	    set(gca,'YAxisLocation','right');
	    set(gca,'TickLength',[0 0])
	    ylabel('(percentile)')
		%annotation('textbox',dim,'String',R2string,'FitBoxToText','on','EdgeColor','none');
	    % exportfig(gcf,'DIC_MvsO','fontmode','fixed','fontsize',12,...
	              % 'color','rgb','renderer','painters')
		figTitle = 'DICsurfcompare2obs';
		print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

%---ALK----------------
	ALK = model.ALK;
    % if ~exist('ALK')
    %     ALK = data.ALK ;
    % end
    % normalize alkalinity for precipitation and eveparation
    nfig = nfig + 1;
    figure(nfig)
    iALK = find(par.alkraw(iwet)>0);
    %
    O = par.alkraw(iwet(iALK));
    M = ALK(iwet(iALK)); % already including anthropogenic CO2
    fprintf('R^2 for ALK is %3.4f \n',rsquare(O,M))
	R2string = sprintf('R^2 = %3.4f   ',rsquare(O,M));
    %
    OvsM = [O, M];
    W = (dVt(iwet(iALK))./sum(dVt(iwet(iALK))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,100,[2300 2300],[2550 2550],W);
    cr = 5:5:95;
    dx = X(3,5)-X(3,4);
    dy = Y(4,2)-Y(3,2);
    [q,ii] = sort(density(:)*dx*dy,'descend');
    D  = density;
    D(ii) = cumsum(q);
    subplot('position',[0.2 0.2 0.6 0.6])
    contourf(X,Y,100*(1-D),cr); hold on
    contour(X,Y,100*(1-D),cr);

    caxis([5 95])
    %set(gca,'FontSize',16);
    grid on
    axis square
    xlabel('Observed ALK (mmol/m^3)');
    ylabel('Model ALK (mmol/m^3)');
    % title('model V.S. observation')
    plot([2300 2550],[2300 2550],'r--','linewidth',2);
	text(2550, 2300, R2string,'VerticalAlignment','bottom','HorizontalAlignment','right','Units','data');

    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    % exportfig(gcf,'ALK_MvsO','fontmode','fixed','fontsize',12,...
    % 'color','rgb','renderer','painters')
	figTitle = 'ALKcompare2obs';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

    if isfield(model,'DOC')
        DOC = model.DOC  ;
    end

    iDOC_ATL = find(DOCclean(iwet)>0 & ATL(iwet)>0) ;
    iDOC_PAC = find(DOCclean(iwet)>0 & PAC(iwet)>0) ;
    iDOC_IND = find(DOCclean(iwet)>0 & IND(iwet)>0) ;
    iDOC_ARC = find(DOCclean(iwet)>0 & ARC(iwet)>0) ;
    iDOC_MED = find(DOCclean(iwet)>0 & MED(iwet)>0) ;

    nfig = nfig + 1;
    figure(nfig)
    plot(DOCclean(iwet(iDOC_ATL)), DOC(iwet(iDOC_ATL)),'ro')
    hold on
    plot(DOCclean(iwet(iDOC_PAC)), DOC(iwet(iDOC_PAC)),'ks')
    hold on
    plot(DOCclean(iwet(iDOC_IND)), DOC(iwet(iDOC_IND)),'b^')
    hold on
    plot(DOCclean(iwet(iDOC_ARC)), DOC(iwet(iDOC_ARC)),'g*')
    hold on
    plot(DOCclean(iwet(iDOC_MED)), DOC(iwet(iDOC_MED)),'c>')
    hold on
    legend('ATL','PAC','IND','ARC','Location','northwest')
    hold on
    plot([0 60],[0 60],'r-','linewidth',3)
    xlim([0 60])
    ylim([0 60])
    hold off
    xlabel('Observed DOC (mmol/m^3)')
    ylabel('Model DOC (mmol/m^3)')
    iDOC = find(DOCclean(iwet)>0 & DOC(iwet)>0) ;
    O = DOCclean(iwet(iDOC)) ;
    M = DOC(iwet(iDOC)) ;
    fprintf('R^2 for DOC is %3.4f \n',rsquare(O,M))
	R2string = sprintf('R^2 = %3.4f   \n', rsquare(O,M));
	text(60, 0, R2string,'VerticalAlignment','bottom','HorizontalAlignment','right','Units','data');
    % OvsM = [O, M] ;

	figTitle = 'DOCcompare2obs';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')

    %%%%%%%% compare to sediment trap %%%%%%%%%%%%
    %POC = data.POC ;
	POC = model.POC ;
    par.kappa_p  = 1/(720*60^2) ;
    if isfield(xhat,'bC_T')
        par.bC   = xhat.bC   ;
        par.bC_T = xhat.bC_T ;
    else
        par.bC   = xhat.bC ;
        par.bC_T = 0 ;
    end

    [PFdiv,Gout] = buildPFD(par,'POC');
    w = Gout.w(:,:,2:25) ;
    fPOC = -w.*POC*spd*12 ;
    fPOC(iarc) = nan ;
    fPOC([1:15,75:end],:,:) = nan;
    ifPOC_ATL = find(POC_flux(iwet)>0 & ATL(iwet)>0) ;
    ifPOC_PAC = find(POC_flux(iwet)>0 & PAC(iwet)>0) ;
    ifPOC_IND = find(POC_flux(iwet)>0 & IND(iwet)>0) ;
    ifPOC_ARC = find(POC_flux(iwet)>0 & ARC(iwet)>0) ;
    ifPOC_MED = find(POC_flux(iwet)>0 & MED(iwet)>0) ;

    nfig = nfig + 1;
    figure(nfig)
    loglog(POC_flux(iwet(ifPOC_ATL)), fPOC(iwet(ifPOC_ATL)),'ro','MarkerFaceColor','red')
    hold on
    loglog(POC_flux(iwet(ifPOC_PAC)), fPOC(iwet(ifPOC_PAC)),'ks','MarkerFaceColor','black')
    hold on
    loglog(POC_flux(iwet(ifPOC_IND)), fPOC(iwet(ifPOC_IND)),'b^','MarkerFaceColor','blue')
    hold on
    loglog(POC_flux(iwet(ifPOC_ARC)), fPOC(iwet(ifPOC_ARC)),'g*','MarkerFaceColor','green')
    hold on
    loglog(POC_flux(iwet(ifPOC_MED)), fPOC(iwet(ifPOC_MED)),'c>','MarkerFaceColor','cyan')
    legend('ATL','PAC','IND','ARC','MED','Location','northwest')
    hold on
    plot([0.1 1000],[0.1 1000],'r','linewidth',2)
    hold off
    xlim([0.1 1000])
    ylim([0.1 1000])
    xlabel('Observed POC flux (mg/m^2/day)')
    ylabel('Model POC flux (mg/m^2/day)')
    %
    ikeep = find(POC_flux(:)>0 & fPOC(:)>0);
    fprintf('R^2 for fPOC is %3.4f \n', ...
            rsquare(log10(POC_flux(ikeep)),log10(fPOC(ikeep))))

    PIC = model.PIC ;
    if isfield(xhat,'d')
        par.d   = xhat.d  ;
    end
    par.tauPIC = 30*spd   ;
    [~,Gout] = buildPFD(par,'PIC') ;
    w = Gout.w(:,:,2:25)  ;
    fPIC = -w.*PIC*spd*12 ;
    % exportfig(gcf,'fPOC','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    rRatio = fPOC./fPIC ;
    load RainRatio.mat rR
    rR([1:15,75:end],:,:) = nan ;
    irRatio_ATL = find(rR(iwet)>0 & rR(iwet)>0 & ATL(iwet)>0) ;
    irRatio_PAC = find(rR(iwet)>0 & rR(iwet)>0 & PAC(iwet)>0) ;
    irRatio_IND = find(rR(iwet)>0 & rR(iwet)>0 & IND(iwet)>0) ;
    irRatio_ARC = find(rR(iwet)>0 & rR(iwet)>0 & ARC(iwet)>0) ;
    irRatio_MED = find(rR(iwet)>0 & rR(iwet)>0 & MED(iwet)>0) ;

    nfig = nfig + 1;
    figure(nfig)
    loglog(rR(iwet(irRatio_ATL)), rRatio(iwet(irRatio_ATL)),'ro','MarkerFaceColor','red')
    hold on
    loglog(rR(iwet(irRatio_PAC)), rRatio(iwet(irRatio_PAC)),'ks','MarkerFaceColor','black')
    hold on
    loglog(rR(iwet(irRatio_IND)), rRatio(iwet(irRatio_IND)),'b^','MarkerFaceColor','blue')
    % hold on
    % loglog(rR(iwet(irRatio_ARC)), rRatio(iwet(irRatio_ARC)),'g*','MarkerFaceColor','green')
    % hold on
    % loglog(rR(iwet(irRatio_MED)), rRatio(iwet(irRatio_MED)),'c>','MarkerFaceColor','cyan')
    legend('ATL','PAC','IND','Location','northwest')
    hold on
    plot([0.1 100],[0.1 100],'r','linewidth',2)
    hold off
    xlim([0.1 100])
    ylim([0.1 100])
    xlabel('Observed rain ratio')
    ylabel('Model rain ratio')
    ikp= find(rR(iwet)>0 & rR(iwet)<100);
    fprintf('R^2 for rain ratio is %3.4f \n', ...
            rsquare(log10(rR(iwet(ikp))),log10(rRatio(iwet(ikp)))))
    % exportfig(gcf,'rRatio','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
	figTitle = 'POCcompare2obs';
	print(gcf,[figPath 'FIG_' figTitle '.png'],'-dpng')
end

% ---------------------------------------------------
if (par.Omodel == on)
    if ~exist('O2')
        %O2 = data.O2 ;
		O2 = model.O2 ;
    end
    nfig = nfig + 1;
    figure(nfig)
    o2obs = o2raw ;

    io2 = find(O2(iwet)>0 & o2obs(iwet)>0);
    %
    M = O2(iwet(io2));
    O = o2obs(iwet(io2));
    fprintf('R^2 for O2 is %3.3f \n',rsquare(O,M))
    OvsM = [O, M];
    %
    W = (dVt(iwet(io2))./sum(dVt(iwet(io2))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,100,[0 0],[300 300],W);
    cr = 5:5:95;
    dx = X(3,5)-X(3,4);
    dy = Y(4,2)-Y(3,2);
    [q,ii] = sort(density(:)*dx*dy,'descend');
    D  = density;
    D(ii) = cumsum(q);
    subplot('position',[0.2 0.2 0.6 0.6])
    contourf(X,Y,100*(1-D),cr); hold on
    contour(X,Y,100*(1-D),cr);

    caxis([5 95])
    %set(gca,'FontSize',16);
    grid on
    axis square
    xlabel('Observed O2 (mmol/m^3)');
    ylabel('Model O2 (mmol/m^3)');
    % title('model V.S. observation')
    plot([0 300],[0 300],'r--','linewidth',2);

    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    % exportfig(gcf,'O2_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
end

% -----------------------------------------------------
if (par.Simodel == on)
    if ~exist(DSi)
        %DSi = data.DSi ;
		DSi = model.DSi ;
    end

    nfig = nfig + 1;
    figure(nfig)
    iDSi = find(DSi(iwet)>0 & sio4raw(iwet)>0);

    cr = 5:5:95;
    M = DSi(iwet(iDSi));
    O = sio4raw(iwet(iDSi));
    OvsM = [O, M];
    fprintf('R^2 for DSi is %3.3f \n',rsquare(O,M))
    W = (dVt(iwet(iDSi))./sum(dVt(iwet(iDSi))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,500,[0 0],[200 200],W);

    dx = X(3,5)-X(3,4);
    dy = Y(4,2)-Y(3,2);
    [q,ii] = sort(density(:)*dx*dy,'descend');
    D = density;
    D(ii) = cumsum(q);
    subplot('position',[0.2 0.2 0.6 0.6])
    contourf(X,Y,100*(1-D),cr); hold on
    contour(X,Y,100*(1-D),cr);

    caxis([5 95])
    %set(gca,'FontSize',16);
    grid on
    axis square
    xlabel('Observed DSi (mmol/m^3)');
    ylabel('Model DSi (mmol/m^3)');
    % title('model V.S. observation')
    plot([0 200],[0 200],'r--','linewidth',2);

    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    % exportfig(gcf,'DSi_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
end
