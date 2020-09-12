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
    input_dir = sprintf('/DFS-L/DATA/primeau/weilewang/COP4WWF/MSK%2d/',GridVer);
end
VER = strcat(input_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if (Cmodel == off & Omodel == off & Simodel == off)
    fname = strcat(VER,'_P');
elseif (Cmodel == on & Omodel == off & Simodel == off)
    fname = strcat(VER,'_PC');
elseif (Cmodel == on & Omodel == on & Simodel == off)
    fname = strcat(VER,'_PCO');
elseif (Cmodel == on & par.Omodel == off & Simodel == on)
    fname = strcat(VER,'_PCSi');
elseif (Cmodel == on & Omodel == on & Simodel == on)
    fname = strcat(VER,'_PCOSi');
end

if GridVer == 90
    load transport_v4.mat
    load M3d90x180x24v2.mat MSKS 
    load GLODAPv2_90x180x24raw.mat
    load DICant_90x180x24.mat DICant
    grd  = grid;
elseif GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    load M3d91x180x24.mat MSKS 
    load GLODAPv2_91x180x24raw.mat
    load DICant_91x180x24.mat DICant
    grd = output.grid;
    M3d = output.M3d;
end
ARC  = MSKS.ARC;
iarc = find(ARC(:));
o2raw(iarc)   = nan;
dicraw(iarc)  = nan ;
po4raw(iarc)  = nan ;
sio4raw(iarc) = nan ; 
load(fname)
%
iwet = find(M3d(:));
nwet = length(iwet);
dVt  = grd.DXT3d.*grd.DYT3d.*grd.DZT3d;
nfig = 0;
%%%%%%%%%%%%%%%%% compare DIP  %%%%%%%%%%%%%%w
if (Pmodel == on)
    if ~exist('DIP')
        DIP = data.DIP ;
    end 
    
    nfig = nfig+1;
    figure(nfig)
    DIPobs = po4raw; 
    ipo4 = find(DIP(iwet)>0 & DIPobs(iwet)>0);
    O = DIPobs(iwet(ipo4));
    M = DIP(iwet(ipo4));
    rsquare(O,M)
    OvsM = [O,M];
    W = (dVt(iwet(ipo4))./sum(dVt(iwet(ipo4))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,200,[0 0],[4 4],W);
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
    % exportfig(gcf,'DIP_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
end 

%% ---------------------------------------------------
if (Omodel == on)
    if ~exist('O2')
        O2 = data.O2 ;
    end 
    nfig = nfig + 1;
    figure(nfig)
    o2obs = o2raw;
    % convert unit form [ml/l] to [umol/l].
    % o2obs = o2obs.*44.661;  
    % o2 correction based on Bianchi et al.(2012) [umol/l] .
    % o2obs = o2obs.*1.009-2.523;
    io2 = find(O2(iwet)>0 & o2obs(iwet)>0);
    %
    M = O2(iwet(io2));
    O = o2obs(iwet(io2));
    rsquare(O,M)
    OvsM = [O, M];
    %
    W = (dVt(iwet(io2))./sum(dVt(iwet(io2))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,200,[0 0],[300 300],W);
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

%% -----------------------------------------------------
if (Cmodel == on)
    if ~exist('DIC')
        DIC = data.DIC ;
    end 
    nfig = nfig + 1;
    figure(nfig)
    rho = 1024.5     ; % seawater density;
    permil = rho*1e-3; % from umol/kg to mmol/m3;
    DICobs = dicraw*permil; % GLODAP dic obs [mmol/m3];
    iDIC = find(DICobs(iwet)>0);
    %
    O = DICobs(iwet(iDIC));
    M = DIC(iwet(iDIC)); % already including anthropogenic CO2 
    rsquare(O,M)
    %
    OvsM = [O, M];
    W = (dVt(iwet(iDIC))./sum(dVt(iwet(iDIC))));
    [bandwidth,density,X,Y] = mykde2d(OvsM,200,[2000 2000],[2500 2500],W);
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
    % exportfig(gcf,'DIC_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
end

%% -----------------------------------------------------
if (Simodel == on)
    if ~exist(SIL)
        SIL = data.SIL ;
    end 
    
    nfig = nfig + 1;
    figure(nfig)
    iSIL = find(SIL(iwet)>0 & sio4raw(iwet)>0);
    
    cr = 5:5:95;
    M = SIL(iwet(iSIL));
    O = sio4raw(iwet(iSIL));
    OvsM = [O, M];
    rsquare(O,M)
    W = (dVt(iwet(iSIL))./sum(dVt(iwet(iSIL))));
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
    xlabel('Observed DIC (mmol/m^3)');
    ylabel('Model DIC (mmol/m^3)');
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
    % exportfig(gcf,'DIC_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
end 