clc; clear alul; close all
global iter
iter = 0 ;
on   = true  ;
off  = false ;
format long
%
GridVer  = 91  ;
operator = 'A' ;
% GridVer: choose from 90 and 91; Ver 90 is for a Transport
% operator without diapycnal mixing but optimized using DIP ;
% Ver 91 include a bunch of operators that include diapycnal
% mixing. These operators represent sensiviity tests on He
% constraint and on mixing parameterizations (DeVries et al, 2018).
% A -> CTL_He; B -> CTL_noHe; C -> KiHIGH_He; D -> KiHIGH_noHe;
% E -> KvHIGH_KiLOW_He; F -> KvHIGH_KiLOW_noHe; G -> KiLOW_He;
% H -> KiLOW_noHe; I -> KvHIGH_He; J -> KvHIGH_noHe; K -> KvHIGH_KiHIGH_noHe
Gtest = off ;
Htest = off ;
par.optim   = on ;
par.Cmodel  = on ;
par.Omodel  = on ;
par.Simodel = off ;
par.LoadOpt = off ; % if load optimial par.
par.pscale  = 0.0 ;
par.cscale  = 0.0 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigma = on ;
par.opt_kP_T  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ;
par.opt_bP    = on ;
par.opt_beta  = on ;
par.opt_alpha = on ;
% C model parameters
par.opt_bC_T  = on ;
par.opt_bC    = on ;
par.opt_d     = on ;
par.opt_kC_T  = on ;
par.opt_kdC   = on ;
par.opt_R_Si  = on ;
par.opt_rR    = on ;
par.opt_cc    = on ;
par.opt_dd    = on ;
% O model parameters
par.opt_O2C_T = on ;
par.opt_rO2C  = on ;
par.opt_O2P_T = on ;
par.opt_rO2P  = on ;
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;

% --------------load data and set up parameters -----------
SetUp ;

% save results
% ATTENTION: Change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer);
elseif isunix
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
    % 'MSK%2d/'],GridVer);
    output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
                        'MSK%2d/PME4DICALK/'],GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
    % 'MSK%2d/'],GridVer);
end
VER = strcat(output_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_P');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
        base_name = strcat(VER,'_PCv2');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
        base_name = strcat(VER,'_PCOv1');
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
end
pfname    = strcat(fname,'_pert.mat');
par.fname = fname ;
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ;
load(fname)
load(fxhat)

%---------------- inital guesses on C and O ---------------
DIC = data.DIC - par.dicant ;

GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); ...
       data.PIC(iwet); data.ALK(iwet)];
% GC  = GC + 1e-6*randn(5*nwet,1) ;
if par.Omodel == on
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end
%--------------------- prepare parameters ------------------
if par.optim == on
    % load optimal parameters from a file or set them to default values
    par = SetPar(par) ;
    % pack parameters into an array, assign them corresponding indices.
    par = PackPar(par) ;
end
p0 = par.p0 ;
% ---------------------perturbe temperature ---------------
vT0  = par.Temp(iwet) ;
vT1  = vT0 + 2        ;
Tz1  = (vT1-min(vT0))/(max(vT0)-min(vT0));
Tz0  = (vT0-min(vT0))/(max(vT0)-min(vT0));
Tz3d = M3d + nan      ;
Tz3d(iwet) = Tz1      ;
par.aveT   = nanmean(Tz3d(:,:,1:3),3) ;
par.Tz     = Tz1*1e-8 ;

DIP = data.DIP ;
POP = data.POP ;
DOP = data.DOP ;

if isfile(pfname)
    load(pfname)
    pDIP = pdata.DIP ;
    pDOP = pdata.DOP ;
    pPOP = pdata.POP ;
else
    [par, P ] = eqPcycle(p0, par) ;
    %
    pDIP = M3d+nan ;  pDIP(iwet) = P(1+0*nwet:1*nwet) ;
    pPOP = M3d+nan ;  pPOP(iwet) = P(1+1*nwet:2*nwet) ;
    pDOP = M3d+nan ;  pDOP(iwet) = P(1+2*nwet:3*nwet) ;
end
par.DIP  = pDIP(iwet) ;
% ---------------------- make plots --------------------
nfig = 0 ;
nfig = nfig + 1 ;
figure(nfig)
dDIP = pDIP - data.DIP ;
% make a zonal cross section of the age
contourf(grd.yt,-grd.zt,squeeze(dDIP(:,170,:))')
set(gca,'color','black')
colormap(darkb2r(-0.015,0.025)), colorbar
xlabel('latitude (deg)');
ylabel('depth (m)')
t = sprintf('DIP anomally x = %4.1f deg', 170);
title(t);

nfig = nfig + 1 ;
figure(nfig)
% make a zonal average of age for the Pacific basin
PAC = MSKS.PAC;
PZA = squeeze(nansum(PAC.*dDIP.*dVt,2)./sum(PAC.*dVt,2))';
subplot(2,1,1) ;
contourf(grd.yt,-grd.zt(1:9),PZA(1:9,:),[-0.015:0.005:0.025]);
set(gca,'color','black')
colormap(darkb2r(-0.015,0.025)), colorbar
ylabel('depth (m)');
title('Pacific zonal average DIP anomaly')
%
subplot(2,1,2) ;
contourf(grd.yt,-grd.zt(10:end),PZA(10:end,:),[-0.015:0.005:0.025]);
set(gca,'color','black')
colormap(darkb2r(-0.015,0.025)), colorbar
caxis([-0.015 0.025])
colorbar
xlabel('latitutde (deg)')
ylabel('depth (m)')
set(gcf, 'InvertHardcopy', 'off')
exportfig(gcf,'Figs91/pza_dip','fontmode','fixed','fontsize',12, ...
          'color','rgb','renderer','painters')
%
nfig = nfig + 1 ;
figure(nfig)
% make a zonal average of age for the Pacific basin
ATL = MSKS.ATL;
AZA = squeeze(nansum(ATL.*dDIP.*dVt,2)./sum(ATL.*dVt,2))';
subplot(2,1,1) ;
contourf(grd.yt,-grd.zt(1:9),AZA(1:9,:),[-0.015:0.005:0.025]);
set(gca,'color','black')
colormap(darkb2r(-0.015,0.025)), colorbar
ylabel('depth (m)')
title('Atlantic zonal average DIP anomaly')
%
subplot(2,1,2) ;
contourf(grd.yt,-grd.zt(10:end),AZA(10:end,:),[-0.015:0.005:0.025]);
set(gca,'color','black')
colormap(darkb2r(-0.015,0.025)), colorbar
xlabel('latitutde (deg)')
ylabel('depth (m)')
set(gcf, 'InvertHardcopy', 'off')
exportfig(gcf,'Figs91/aza_dip','fontmode','fixed','fontsize',12, ...
          'color','rgb','renderer','painters')

nfig = nfig + 1 ;
figure(nfig)
pcolor(nanmean(dDIP(:,:,1:3),3));colorbar;shading flat
set(gca,'color','black')
colormap(darkb2r(-0.01, 0.05)), colorbar
set(gcf, 'InvertHardcopy', 'off')
exportfig(gcf,'Figs91/surface_dip_anomaly','fontmode','fixed','fontsize',12, ...
          'color','rgb','renderer','painters')

if Cmodel == on
    if isfile(pfname)
        pDIC = pdata.DIC ;
        pDOC = pdata.DOC ;
        pPOC = pdata.POC ;
        pPIC = pdata.PIC ;
    else
        DIC = data.DIC ;
        POC = data.POC ;
        DOC = data.DOC ;
        PIC = data.PIC ;

        [par, C ] = eqCcycle(p0, par) ;
        pDIC = M3d+nan ;  pDIC(iwet) = C(0*nwet+1:1*nwet) ;
        pPOC = M3d+nan ;  pPOC(iwet) = C(1*nwet+1:2*nwet) ;
        pDOC = M3d+nan ;  pDOC(iwet) = C(2*nwet+1:3*nwet) ;
        pPIC = M3d+nan ;  pPIC(iwet) = C(3*nwet+1:4*nwet) ;
        pALK = M3d+nan ;  pALK(iwet) = C(4*nwet+1:5*nwet) ;
        par.DIC  = pDIC(iwet) ;
        par.DOC  = pDOC(iwet) ;
    end

    dDIC = pDIC + par.dicant - data.DIC ;
    nfig = nfig + 1 ;
    figure(nfig)
    % make a zonal cross section of the age
    contourf(grd.yt,-grd.zt,squeeze(dDIC(:,170,:))')
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    xlabel('latitude (deg)');
    ylabel('depth (m)')
    t = sprintf('DIC anomaly x = %4.1f deg', 170);
    title(t);

    nfig = nfig + 1 ;
    figure(nfig)
    % make a zonal average of age for the Pacific basin
    PAC = MSKS.PAC;
    PZA = squeeze(nansum(PAC.*dDIC.*dVt,2)./sum(PAC.*dVt,2))';
    subplot(2,1,1) ;
    contourf(grd.yt,-grd.zt(1:9),PZA(1:9,:),[-15:2:15]);
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    ylabel('depth (m)');
    title('Pacific zonal average DIP anomaly')
    %
    subplot(2,1,2) ;
    contourf(grd.yt,-grd.zt(10:end),PZA(10:end,:),[-15:2:15]);
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    xlabel('latitutde (deg)');
    ylabel('depth (m)');
    set(gcf, 'InvertHardcopy', 'off')
    exportfig(gcf,'Figs91/pza_dic','fontmode','fixed','fontsize',12, ...
              'color','rgb','renderer','painters')

    nfig = nfig + 1 ;
    figure(nfig)
    % make a zonal average of age for the Pacific basin
    ATL = MSKS.ATL;
    AZA = squeeze(nansum(ATL.*dDIC.*dVt,2)./sum(ATL.*dVt,2))';
    subplot(2,1,1) ;
    contourf(grd.yt,-grd.zt(1:9),AZA(1:9,:),[-15:2:15]);
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    ylabel('depth (m)');
    title('Atlantic zonal average DIC anomaly')
    %
    subplot(2,1,2) ;
    contourf(grd.yt,-grd.zt(10:end),AZA(10:end,:),[-15:2:15]);
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    xlabel('latitutde (deg)');
    ylabel('depth (m)');
    set(gcf, 'InvertHardcopy', 'off')
    exportfig(gcf,'Figs91/aza_dic','fontmode','fixed','fontsize',12, ...
              'color','rgb','renderer','painters')

    nfig = nfig + 1 ;
    figure(nfig)
    pcolor(nanmean(dDIC(:,:,1:3),3));colorbar;shading flat
    set(gca,'color','black')
    colormap(darkb2r(-5, 20)), colorbar
    set(gcf, 'InvertHardcopy', 'off')
    exportfig(gcf,'Figs91/surface_dic_anomaly','fontmode','fixed','fontsize',12, ...
              'color','rgb','renderer','painters')

end

if Omodel == on
    O2 = data.O2 ;
    GO = real(O2(iwet)) + 1e-9*randn(par.nwet,1);
    [par, O2] = eqOcycle(p0, par) ;
end
