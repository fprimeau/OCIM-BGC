clc; clear all; close all
spd  = 24*60^2 ; spa  = 365*spd ;
on = true; off = false;
%
expFigs  = off ;
GridVer  = 91 ;
operator = 'A';

par.optim   = off ; % on: do optimization
par.Pmodel  = on ; % on: run P model ;
par.Cmodel  = on ; % on: run C model ;
par.Omodel  = on ; % on: run O model; 
par.Simodel = off ; % on: run Si model;
par.LoadOpt = off ; % on: load optimial parameters;
                    % factor to weigh DOP in the objective function
par.pscale  = 0.0 ;
% factor to weigh DOC in the objective function
par.cscale  = 1.0 ; 

%-------------load data and set up parameters---------------------
SetUp ;

if ismac
    input_dir = sprintf('../DATA/MSK%2d/',GridVer); 
elseif isunix
    input_dir = sprintf('~/rDOC-OP/MSK%2d/',GridVer) ; 
end
VER = strcat(input_dir,TRdivVer);

% Creat output file names based on which model(s) is(are) optimized
if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
    fname = strcat(VER,'_Pv3');
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
    base_name = strcat(VER,'_PCv1'); 
    catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
    fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta');
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta');
    base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta');
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_MODIS_CbPM_aveTeu_diffSig_O2C_uniEta');
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta_noArcMed');
    catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
    fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
    base_name = strcat(VER,'_PCSi');
    catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
    fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
    base_name = strcat(VER,'_PCOSi');
    catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
    fname = strcat(base_name,catDOC);
end
par.fname = strcat(fname,'.mat') ; 
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ; 
load(par.fxhat) ;
load(par.fname) ;
%------------------ compare DIP ---------------------------------
nfig = 0;
if (par.Pmodel == on)
    if ~exist('DOP')
        DOP = data.DOP ;
    end 
    
    dopraw = DOPobs - 0.03 ; % less refractory DOP 
    idop   = find(dopraw(:) > 0 & DOP(:)>0) ;
    iDOP_ATL = find(dopraw(iwet)>0 & ATL(iwet)>0) ;
    iDOP_PAC = find(dopraw(iwet)>0 & PAC(iwet)>0) ;
    iDOP_IND = find(dopraw(iwet)>0 & IND(iwet)>0) ;
    iDOP_ARC = find(dopraw(iwet)>0 & ARC(iwet)>0) ;
    iDOP_MED = find(dopraw(iwet)>0 & MED(iwet)>0) ;
    fprintf('R^2 for DOP is %3.3f \n',rsquare(dopraw(idop),DOP(idop)))
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
    plot([0 0.75],[0 0.75],'r-','linewidth',3)
    xlim([0 0.75])
    ylim([0 0.75])

    if ~exist('DIP')
        DIP = data.DIP ;
    end 
    nfig = nfig+1;
    figure(nfig)
    ipo4 = find(DIP(iwet) > 0 & po4raw(iwet) > 0.05);

    O = po4raw(iwet(ipo4));
    M = DIP(iwet(ipo4));

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
    xlabel('Observed DIP (mmol/m^3)');
    ylabel('Model DIP (mmol/m^3)');

    % title('model V.S. observation')
    plot([0 4],[0 4],'r--','linewidth',2);
    R2 = sprintf('%3.2f',rsquare(O,M));
    fprintf('valid DIP data points %3.0d \n', length(O))
    fprintf('R^2 for DIP is %3.3f \n',str2num(R2))
    text('string',['N = ' num2str(length(O))],'units','normalized','position',[0.25,0.95],'fontsize',12)
    text('string',['R^2 = ' R2],'units','normalized','position',[0.25,0.90],'fontsize',12)

    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    if expFigs == on 
        exportfig(gcf,'mkFigs/SeaWiFS/DIP_MvsO_1to100','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    end 
end

% -----------------------------------------------------
if (par.Cmodel == on)
    nfig = nfig + 1;
    figure(nfig)
    DIC = data.DIC - par.dicant ;
    DICobs = par.dicraw ;
    iDIC = find(DICobs(iwet)>0);
    %
    O = DICobs(iwet(iDIC));
    % already including anthropogenic CO2
    % M = DIC(iwet(iDIC)) ;
    % not include anthropogenic CO2
    M = DIC(iwet(iDIC))+par.dicant(iwet(iDIC)); 
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
    % title('model V.S. observation')
    plot([2000 2500],[2000 2500],'r--','linewidth',2);
    R2 = sprintf('%3.2f',rsquare(O,M));
    fprintf('valid DIC data points %3.0d \n', length(O))
    fprintf('R^2 for DIC is %3.3f \n',str2num(R2))
    text('string',['N = ' num2str(length(O))],'units','normalized','position',[0.25,0.95],'fontsize',12)
    text('string',['R^2 = ' R2],'units','normalized','position',[0.25,0.90],'fontsize',12)

    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 1.5],cr,[cr(:),cr(:)],cr); hold on
    contour([1 1.5],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    if expFigs == on 
        exportfig(gcf,'mkFigs/SeaWiFS/DIC_MvsO_1to100','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    end
    
    if ~exist('ALK')
        ALK = data.ALK ;
    end
    % normalize alkalinity for precipitation and eveparation
    nfig = nfig + 1;
    figure(nfig)
    iALK = find(par.alkraw(iwet)>0);
    %
    O = par.alkraw(iwet(iALK));
    M = ALK(iwet(iALK)); % already including anthropogenic CO2 
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
    R2 = sprintf('%3.2f',rsquare(O,M));
    fprintf('valid ALK data points % 3.0d \n', length(O))
    fprintf('R^2 for ALK is %3.3f \n',str2num(R2))
    text('string',['N = ' num2str(length(O))],'units','normalized','position',[0.25,0.95],'fontsize',12)
    text('string',['R^2 = ' R2],'units','normalized','position',[0.25,0.90],'fontsize',12)
    
    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    if expFigs == on 
        exportfig(gcf,'mkFigs/SeaWiFS/ALK_MvsO_1to100','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    end
    
    if isfield(data,'DOC') 
        DOC = data.DOC + data.DOCr ;%+ data.DOCl;
    end 
    
    iDOC_ATL = find(par.docraw(iwet)>0 & ATL(iwet)>0) ;
    iDOC_PAC = find(par.docraw(iwet)>0 & PAC(iwet)>0) ;
    iDOC_IND = find(par.docraw(iwet)>0 & IND(iwet)>0) ;
    iDOC_ARC = find(par.docraw(iwet)>0 & ARC(iwet)>0) ;
    iDOC_MED = find(par.docraw(iwet)>0 & MED(iwet)>0) ;
    
    nfig = nfig + 1;
    figure(nfig)
    plot(par.docraw(iwet(iDOC_ATL)), DOC(iwet(iDOC_ATL)),'ro')
    hold on
    plot(par.docraw(iwet(iDOC_PAC)), DOC(iwet(iDOC_PAC)),'ks')
    hold on
    plot(par.docraw(iwet(iDOC_IND)), DOC(iwet(iDOC_IND)),'b^')
    hold on
    % plot(par.docraw(iwet(iDOC_ARC)), DOC(iwet(iDOC_ARC)),'g*')
    % hold on
    % plot(par.docraw(iwet(iDOC_MED)), DOC(iwet(iDOC_MED)),'c>')
    % hold on 
    legend('ATL','PAC','IND','Location','northwest')
    hold on
    plot([30 90],[30 90],'r-','linewidth',3)

    xlim([30 90])
    ylim([30 90])
    hold off
    xlabel('Observed DOC (mmol/m^3)')
    ylabel('Model DOC (mmol/m^3)')
    iDOC = find(par.docraw(iwet)>0 & DOC(iwet)>0) ;
    O = par.docraw(iwet(iDOC)) ;
    M = DOC(iwet(iDOC)) ;
    R2 = sprintf('%3.2f',rsquare(O,M));
    fprintf('valid DOC data points %3.0d \n', length(O))
    fprintf('R^2 for DOC is %3.2f \n', str2num(R2))
    text('string',['N = ' num2str(length(O))],'units','normalized','position',[0.5,0.95],'fontsize',12)
    text('string',['R^2 = ' R2],'units','normalized','position',[0.5,0.90],'fontsize',12)

    if expFigs == on 
        exportfig(gcf,'mkFigs/SeaWiFS/tDOC_MvsO_1to100','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    end
    %%%%%%%% compare to sediment trap %%%%%%%%%%%%
    POC = data.POC ;
    par.kappa_p  = 1/(720*60^2) ;
    if isfield(xhat,'bC_T')
        par.bC   = xhat.bC   ;
        par.bC_T = xhat.bC_T ;
    else
        par.bC   = xhat.bC ;
        par.bC_T = 0 ;
    end         
    
    [PFdiv,Gout] = buildPFD(par,'POC');
    w = Gout.w(:,:,1:24) ;
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
    R2 = sprintf('%3.2f',rsquare(log10(POC_flux(ikeep)),log10(fPOC(ikeep))));
    fprintf('valid fPOC data points %3.0d \n', length(POC_flux(ikeep)))
    fprintf('R^2 for fPOC is %3.2f \n', str2num(R2))
    text('string',['N = ' num2str(length(fPOC(ikeep)))],'units','normalized','position',[0.5,0.95],'fontsize',12)
    text('string',['R^2 = ' R2],'units','normalized','position',[0.5,0.90],'fontsize',12)

    PIC = data.PIC ;
    if isfield(xhat,'d')
        par.d   = xhat.d  ;
    end         
    par.tauPIC = 30*spd   ;
    [~,Gout] = buildPFD(par,'PIC') ;
    w = Gout.w(:,:,1:24)  ;
    fPIC = -w.*PIC*spd*12 ;
    if expFigs == on 
        exportfig(gcf,'mkFigs/SeaWiFS/fPOC_MvsO_1to100','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    end
    
    rRatio = fPOC./fPIC ;
    load RainRatio.mat rR
    rR([1:15,75:end],:,:) = nan ;
    irRatio_ATL = find(rR(iwet)>0 & rR(iwet)>0 & ATL(iwet)>0) ;
    irRatio_PAC = find(rR(iwet)>0 & rR(iwet)>0 & PAC(iwet)>0) ;
    irRatio_IND = find(rR(iwet)>0 & rR(iwet)>0 & IND(iwet)>0) ;
    % irRatio_ARC = find(rR(iwet)>0 & rR(iwet)>0 & ARC(iwet)>0) ;
    % irRatio_MED = find(rR(iwet)>0 & rR(iwet)>0 & MED(iwet)>0) ;
    
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
    R2 = sprintf('%3.2f',rsquare(log10(rR(iwet(ikp))),log10(rRatio(iwet(ikp)))));
    fprintf('valid rRR data points %3.0d \n', length(rR(iwet(ikp))))
    fprintf('R^2 for rain ratio is %3.3f \n', str2num(R2))
    text('string',['N = ' num2str(length(rR(iwet(ikp))))],'units','normalized','position',[0.25,0.95],'fontsize',12)
    text('string',['R^2 = ' R2],'units','normalized','position',[0.75,0.90],'fontsize',12)

    if expFigs == on 
        exportfig(gcf,'mkFigs/SeaWiFS/rRatio_MvsO_1to100','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    end
end 

% ---------------------------------------------------
if (par.Omodel == on)
    if ~exist('O2')
        O2 = data.O2 ;
    end 
    nfig = nfig + 1;
    figure(nfig)
    o2obs = o2raw ;
    
    io2 = find(O2(iwet)>0 & o2obs(iwet)>0);
    %
    M = O2(iwet(io2));
    O = o2obs(iwet(io2));
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
    R2 = sprintf('%3.2f',rsquare(O,M)) ;
    fprintf('valid O2 data points %3.0d \n', length(O))
    fprintf('R^2 for O2 is %3.3f \n',str2num(R2)) 
    text('string',['N = ' num2str(length(O))],'units','normalized','position',[0.25,0.95],'fontsize',12)
    text('string',['R^2 = ' R2],'units','normalized','position',[0.25,0.90],'fontsize',12)
    
    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    if expFigs == on 
        exportfig(gcf,'mkFigs/SeaWiFS/O2_MvsO_1to100','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    end
end 

% -----------------------------------------------------
if (par.Simodel == on)
    if ~exist(DSi)
        DSi = data.DSi ;
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
    text('string',['R^2 = ' num2str(rsquare(O,M))],'units','normalized','position',[0.75,0.95],'fontsize',12)
    text('string',['N = ' num2str(length(O))],'units','normalized','position',[0.25,0.85],'fontsize',12)
    
    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    if expFigs == on 
        exportfig(gcf,'mkFigs/SeaWiFS/DSi_MvsO_1to100','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
    end 
end