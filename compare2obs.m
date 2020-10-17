clc; clear all; close all
if ismac 
    addpath('~/Dropbox/myfunc'     )
    addpath('~/Documents/DATA/'    )
    addpath('~/Documents/DATA/OCIM')
else 
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/');
    addpath('/DFS-L/DATA/primeau/weilewang/my_func');
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
end 
on = true; off = false;
TR_ver = 91 ;
%
Pmodel  = on ;
Cmodel  = on ;
Omodel  = off ;
Simodel = off ;
fscale  = 1.0 ; % factor to weigh DOC in the objective function
                %
GridVer   = 91 ;
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

if ismac
    input_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer); 
elseif isunix
    input_dir = sprintf('/DFS-L/DATA/primeau/weilewang/TempSensi/MSK%2d/',GridVer);
end
VER = strcat(input_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if (Cmodel == off & Omodel == off & Simodel == off)
    fname = strcat(VER,'_P');
elseif (Cmodel == on & Omodel == off & Simodel == off)
    base_name = strcat(VER,'_PCv2Zscore'); %minmax_PME4All');
    catDOC = sprintf('_DOC%2.0e',fscale);
    fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == on & Simodel == off)
    base_name = strcat(VER,'_PCOv3');
    catDOC = sprintf('_DOC%2.0e',fscale);
    fname = strcat(base_name,catDOC);
elseif (Cmodel == on & par.Omodel == off & Simodel == on)
    base_name = strcat(VER,'_PCSi');
    catDOC = sprintf('_DOC%2.0e',fscale);
    fname = strcat(base_name,catDOC);
elseif (Cmodel == on & Omodel == on & Simodel == on)
    base_name = strcat(VER,'_PCOSi');
    catDOC = sprintf('_DOC%2.0e',fscale);
    fname = strcat(base_name,catDOC);
end

if GridVer == 90
    load transport_v4.mat
    load M3d90x180x24v2.mat MSKS
    load Sobs_90x180x24.mat
    load GLODAPv2_90x180x24raw.mat
    load PME_TS_90x180x24.mat modT modS
    load DICant_90x180x24.mat DICant
    grd  = grid;
elseif GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    load M3d91x180x24.mat MSKS
    load Sobs_91x180x24.mat
    load PME_TS_91x180x24.mat modT modS
    load GLODAPv2_91x180x24raw.mat
    load DICant_91x180x24.mat DICant
    load DOMobs_91x180x24.mat
    grd = output.grid;
    M3d = output.M3d;
end
ATL = MSKS.ATL ;
PAC = MSKS.PAC ;
IND = MSKS.IND ;
ARC = MSKS.ARC ;
MED = MSKS.MED ;

iarc = find(ARC(:)) ;
imed = find(MED(:)) ;
% DOCobs(iarc)  = nan ;
% DOCobs(imed)  = nan ;
alkraw(iarc)  = nan ;
alkraw(imed)  = nan ; 
dicraw(iarc)  = nan ;
dicraw(imed)  = nan ;
% o2raw(iarc)   = nan ;
% o2raw(imed)   = nan ;
% po4raw(iarc)  = nan ;
% po4raw(imed)  = nan ; 
% sio4raw(iarc) = nan ;
% sio4raw(imed) = nan ;

load(fname)
%
rho = 1024.5     ; % seawater density;
permil = rho*1e-3; % from umol/kg to mmol/m3;

iwet = find(M3d(:));
nwet = length(iwet);
dVt  = grd.DXT3d.*grd.DYT3d.*grd.DZT3d;
par.M3d  =  M3d ;
par.grd  = grd  ;
par.iwet = iwet ;
par.modS = modS ;
par.MSKS = MSKS ;
par.Salt = Sobs ;
par.o2raw = o2raw ;
par.po4raw = po4raw  ;
par.sio4raw = sio4raw;

par.dicraw = dicraw*permil ;
par.alkraw = alkraw*permil ;
par.dicant = DICant*permil;

nfig = 0;
%%%%%%%%%%%%%%%%% compare DIP  %%%%%%%%%%%%%%w
if (Pmodel == on)
    if ~exist('DOP')
        DOP = data.DOP ;
    end 
    
    dopraw = DOPobs - 0.03 ; % less refractory DOP 
    idop   = find(dopraw(iwet) > 0.0) ;
    iDOP_ATL = find(dopraw(iwet)>0 & ATL(iwet)>0) ;
    iDOP_PAC = find(dopraw(iwet)>0 & PAC(iwet)>0) ;
    iDOP_IND = find(dopraw(iwet)>0 & IND(iwet)>0) ;
    iDOP_ARC = find(dopraw(iwet)>0 & ARC(iwet)>0) ;
    iDOP_MED = find(dopraw(iwet)>0 & MED(iwet)>0) ;
    fprintf('R^2 for DIP is %3.3f \n',rsquare(dopraw(idop),DOP(idop)))
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
    ipo4 = find(DIP(iwet)>0 & po4raw(iwet)>0);
    O = po4raw(iwet(ipo4));
    M = DIP(iwet(ipo4));
    fprintf('R^2 for DIP is %3.3f \n',rsquare(O,M))
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
end 

% -----------------------------------------------------
if (Cmodel == on)
    if ~exist('DIC')
        DIC = data.DIC ;
    end 
    nfig = nfig + 1;
    figure(nfig)
    
    DICobs = par.dicraw ;
    iDIC = find(DICobs(iwet)>0);
    %
    O = DICobs(iwet(iDIC));
    % already including anthropogenic CO2
    M = DIC(iwet(iDIC)) ;
    % not include anthropogenic CO2
    % M = DIC(iwet(iDIC))+par.dicant(iwet(iDIC)); 
    fprintf('R^2 for DIC is %3.3f \n',rsquare(O,M))
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
    
    subplot('position',[0.82 0.2 0.05 0.6]);
    contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
    contour([1 2],cr,[cr(:),cr(:)],cr);
    hold off
    %set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    % exportfig(gcf,'DIC_MvsO','fontmode','fixed','fontsize',12,...
              % 'color','rgb','renderer','painters')

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
    fprintf('R^2 for ALK is %3.3f \n',rsquare(O,M))
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

    if isfield(data,'DOC') 
        DOC = data.DOC  ;
    end 
    par.DOCobs = DOCobs ;
    DOCclean = RemoveRef(par) ;
    
    DOCclean(DOCclean(:)>0);
    
    ibad = find( DOCclean(iarc) > 50 ) ;
    DOCclean(iarc(ibad)) = nan ;
    
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
    plot([0 60],[0 60],'r-','linewidth',3)
    xlim([0 60])
    ylim([0 60])
    
    nfig = nfig + 1;
    figure(nfig)

    iDOC = find(DOCclean(iwet)>0 & DOC(iwet)>0) ;
    O = DOCclean(iwet(iDOC)) ;
    M = DOC(iwet(iDOC)) ; 
    fprintf('R^2 for DOC is %3.3f \n',rsquare(O,M))
    OvsM = [O, M] ;
    
    Woc  = dVt(iwet(iDOC))/sum(dVt(iwet(iDOC))) ;
    % mu   = sum(Woc*DOCclean(iwet(iDOC)))/sum(diag(Woc)) ;
    % var  = sum(Woc*(DOCclean(iwet(iDOC))-mu).^2)/sum(diag(Woc));
    % Woc  = full(fscale*Woc/var) ;

    [bandwidth,density,X,Y] = mykde2d(OvsM,500,[0 0],[50 50],Woc);
    cr = 5:5:95 ;
    dx = X(3,5)-X(3,4) ; 
    dy = Y(4,2)-Y(3,2) ;
    [q,ii] = sort(density(:)*dx*dy,'descend') ;
    D  = density ;
    D(ii) = cumsum(q) ;
    subplot('position',[0.2 0.2 0.6 0.6])
    contourf(X,Y,100*(1-D),cr) ; hold on
    contour(X,Y,100*(1-D),cr)  ;
    
    caxis([5 95])
    set(gca,'FontSize',16);
    grid on
    axis square
    xlabel('Observed DOC (mmol/m^3)') ;
    ylabel('Model DOC (mmol/m^3)') ;
    title('model V.S. observation')
    plot([0 50],[0 50],'r--','linewidth',2) ;
    
    subplot('position',[0.82 0.2 0.05 0.6]) ;
    contourf([1 2],cr,[cr(:),cr(:)],cr) ; hold on
    contour([1 2],cr,[cr(:),cr(:)],cr)  ;
    hold off
    set(gca,'FontSize',14);
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right');
    set(gca,'TickLength',[0 0])
    ylabel('(percentile)')
    % exportfig(gcf,'DOC_MvsO','fontmode','fixed','fontsize',12,...
              % 'color','rgb','renderer','painters')
end

% ---------------------------------------------------
if (Omodel == on)
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
if (Simodel == on)
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