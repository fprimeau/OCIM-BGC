clc; clear all; close all
global iter
iter = 0 ;
on   = true  ;
off  = false ;
format long
expfig = off ;
% 
GridVer  = 91  ;
operator = 'A' ;

par.optim   = off ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.LoadOpt = on ; % if load optimial par. 
par.pscale  = 0.0 ;
par.cscale  = 1.0 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigP  = on ; 
par.opt_Q10P  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ; 
par.opt_bP    = on ;
par.opt_alpha = on ;
par.opt_beta  = on ;
% C model parameter
par.opt_sigC  = on ; 
par.opt_kru   = on ;
par.opt_krd   = on ;
par.opt_etau  = on ;
par.opt_etad  = off ;
par.opt_bC_T  = on ;
par.opt_bC    = on ; 
par.opt_d     = on ;
par.opt_Q10C  = on ;
par.opt_kdC   = on ; 
par.opt_R_Si  = on ; 
par.opt_rR    = on ; 
par.opt_cc    = on ;
par.opt_dd    = on ;
% O model parameters
par.opt_O2C_T = on ;
par.opt_rO2C  = on ;
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;
%
%-------------load data and set up parameters---------------------
SetUp ;
load MLD_CESM_91x180.mat
R = par.R ;
% save results 
% ATTENTION: Change this direcrtory to where you wanna save your output files
if ismac
    output_dir = sprintf('~/OneDrive/rDOC/MSK%2d/',GridVer); 
elseif isunix
    output_dir = sprintf('~/rDOC-OP/MSK%2d/',GridVer) ;
end
VER = strcat(output_dir,TRdivVer);

% Creat output file names based on which model(s) is(are) optimized
if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
    fname = strcat(VER,'_Pv1');
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
    base_name = strcat(VER,'_PCv3');
    catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale, par.pscale);
    fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
    % base_name = strcat(VER,'_PCO_gamma_half_VGPM_3layers_GM15_aveT1000_diffSigma_C2O_DOP1to5');
    % base_name = strcat(VER,'_PCO_gamma_1third_POC2DIC_GM15_VGPM_aveT1000_diffSigma_O2C');
    % base_name = strcat(VER,'_PCO_Gamma1to2_POC2DIC_GM15_VGPM_aveTeu_diffSigma_O2C_grd90')
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSigma_O2C_noArcMed_grd90');
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta') ;
    % base_name = strcat(VER,'_PCO_gamma_1third_POC2DIC_GM15_VGPM_aveT1000_diffSigma_O2C_uniEta');
    % base_name = strcat(VER,'_PCO_gamma_1third_POC2DIC_GM15_VGPM_aveTeu_diffSigma_O2C_uniEta') ;
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta_lkd1d');
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta_kp60d');
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_MODIS_CbPM_aveTeu_diffSig_O2C_uniEta');
    % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta');
    base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta');
    % base_name = strcat(VER,'_PCO_Gamma1to100_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta');
    catDOC = sprintf('_DOC%2.0e_DOP%2.0e', par.cscale, par.pscale);
    fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
    base_name = strcat(VER,'_PCSi');
    catDOC = sprintf('_DOC%2.0e_DOP%2.0e', par.cscale, par.pscale);
    fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
    base_name = strcat(VER,'_PCOSi');
    catDOC = sprintf('_DOC%2.0e_DOP%2.0e', par.cscale, par.pscale);
    fname = strcat(base_name,catDOC);
end
par.fname = strcat(fname,'.mat') ; 
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');

par.fxhat = fxhat ; 
load(par.fname) ;
load(par.fxhat) ;
%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par) ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

%------------------ extract parameters ---------------------------
% production allocation to DOP/DOC
gamma = par.gamma ;
sigP  = par.sigP  ;
sigC  = par.sigC  ;
kru   = par.kru   ;
krd   = par.krd   ;
% linear parameter of npp to DIP assimilation function. 
alpha = par.alpha ; 
% exponential parameter of npp to DIN assimilation function.
beta  = par.beta ; 
%------------------ prepare NPP for the model --------------------
% DIP assimilation
LAM  = M3d + nan;
LAM(:,:,1:par.nl)  = par.npp.^beta.*par.Lambda(:,:,1:par.nl);
L   = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];

%------------------ extract data ---------------------------------
PO4 = po4obs(iwet)   ;  
DIP = data.DIP(iwet) ;
% DIP(DIP<0.02) = 0.02   ;
if par.Cmodel == on 
    DIC  = data.DIC(iwet)  ;
    POC  = data.POC(iwet)  ;
    DOC  = data.DOC(iwet)  ;
    PIC  = data.PIC(iwet)  ;
    ALK  = data.ALK(iwet)  ;
    DOCr = data.DOCr(iwet) ;
    DOCl = data.DOCl(iwet) ;
end 
TRdiv = par.TRdiv ;
I     = par.I     ;
% -------------- C:P uptake ratio --------------------------------
W = d0(dVt(iwet)) ;
C2P3D = M3d + nan ;
C2P3D(iwet) = 1./(par.cc*PO4 + par.dd) ;

%--------------- calculate primary production --------------------
G       = M3d+nan     ;
G(iwet) = alpha*L*DIP ; % primary production in P unit.
[NPP,Int_PNPP,Int_CNPP] = deal(M3d(:,:,1)*0);
Int_PNPP = nansum(G.*dVt*spa,3) ;
Int_CNPP = nansum(G.*dVt.*C2P3D*spa,3) ;
NPP = nansum(G.*C2P3D.*grd.DZT3d*spa,3) ;

% unit conversion from mg C/m^3/s to g C/m^2/year;
tmp_CNPP = Int_CNPP*12*1e-18 ;
Sum_CNPP = nansum(tmp_CNPP(:)) ;
fprintf('Model NPP is %3.3e \n',Sum_CNPP) ;

%---------------- calculate phosphorus export --------------------
PFD = buildPFD(par, 'POP') ; 

F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;

[nx,ny,nz] = size(M3d) ;
nppMSK = par.nppMSK ;
ieup = find(nppMSK(:) > 0) ;
idrk = find(nppMSK(:) == 0);
[euMSK,dkMSK] = deal(M3d * 0);
euMSK(ieup) = 1 ;
dkMSK(idrk) = 1 ;

% POP export 
% mmol/year
POPexp = nansum(par.kappa_p*data.POP.*dVt.*dkMSK*spa,3) ;
% mmol/m2/year, save to file 
POPexp = POPexp./dAt(:,:,3);

% DOP exp 
kappa_l = par.kappa_l ;
Omega = dkMSK(iwet) ;
Prod  = G(iwet)    ;
% adjoint method to calculate total P export.
tf    = (par.vT - 30)/10 ;
kP    = d0( par.kdP * par.Q10P .^ tf ) ;
T_diag   = inv(W)*TRdiv'*W ;
kP_diag  = inv(W)*kP'*W    ; 

Jex_P = sigP*d0(kP*Prod)*((T_diag + kP_diag)\Omega); 
Jex_Pl = gamma*kappa_l*d0(Prod)*((T_diag + kappa_l*I)\Omega);

P3d = M3d+nan;
P3d(iwet) = Jex_P + Jex_Pl;

%  mmol P/m^3/s 
DOPexp = P3d.*grd.DZT3d*spa ;
DOPexp = nansum(DOPexp, 3) ; % mmol/m2/year, save to file 
tmp_Pexp = P3d .* dVt*spa ;
%---------------- calculate carbon export -------------------------
[PFD, Gout] = buildPFD(par, 'POC') ; 
kC = d0( par.kdC * par.Q10C .^ tf ) ;
bC = xhat.bC_T*par.aveT + xhat.bC ;

% PIC export
[~,Gout] = buildPFD(par,'PIC') ;
w = -Gout.w ;
PICexp = M3d(:,:,1)*0 ;
for ji = 1 : par.nl
    PICexp = PICexp + data.PIC(:,:,ji).*w(:,:,ji).*12*spd ;
end
PICexp = PICexp/12*365 ;  % mmol/m2/year, save to file 
tmp_PICexp = PICexp.*dAt(:,:,3);
Sum_PICexp = nansum(tmp_PICexp(:))*12*1e-18;
fprintf('Model PIC export is %3.3e \n\n',Sum_PICexp);

% adjoint method to calculate total C export.
F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;
kC_diag  = inv(W)*kC'*W    ; 
Prod  = G(iwet).*C2P3D(iwet) ;
% DOCr export is 6 magnitude lower than the other two due to its long decay time.
% Here DOCr export is ignored.
UM = par.UM ;
DM = par.DM ;
WM = par.WM ;
M200 = M3d;
M200(:,:,1:5) = 0; 
kappa_r =  kru*UM +  krd*DM ;

Jex_Cr = krd*DM*DOCr ;
Jex_Cl_200 = kappa_l*DOCl;
Jex_C_200  = kC*DOC ;
Jex_C = sigC*d0(kC*Prod)*((T_diag + kC_diag)\Omega); 
Jex_Cl = gamma*kappa_l*d0(Prod)*((T_diag + kappa_l*I)\Omega);

[lDOCexp200,sDOCexp200] = deal(M3d + nan) ;
[lDOCexp3D,sDOCexp3D,rDOCexp3D] = deal(M3d + nan) ;
lDOCexp3D(iwet) = Jex_Cl ;
sDOCexp3D(iwet) = Jex_C  ;
rDOCexp3D(iwet) = Jex_Cr ;

lDOCexp200(iwet) = Jex_Cl_200;
sDOCexp200(iwet) = Jex_C_200 ;
% DOC export (mmol/year)
lDOCexp = nansum(lDOCexp3D.*dVt.*euMSK*spa,3) ;
sDOCexp = nansum(sDOCexp3D.*dVt.*euMSK*spa,3) ;
rDOCexp = nansum(rDOCexp3D.*dVt*spa,3) ;

int_lDOCexp200 = lDOCexp200.*dVt.*M200*spa; % 3D
int_sDOCexp200 = sDOCexp200.*dVt.*M200*spa; % 3D
% POC export (mmol/year)
POCexp = nansum(par.kappa_p*data.POC.*dVt.*dkMSK*spa,3) ;
[x,y] = size(POCexp) ;
for kk = 1 : x
    for tt = 1 : y
        tmp = squeeze(M3d(kk,tt,:)) ;
        itmp = find(tmp == 1) ;
        if length(itmp) > 12
            POCadj(kk,tt) = POCexp(kk,tt).*(100/grd.zw(4)).^(-bC(kk,tt)) ;
            sDOCadj(kk,tt) = sDOCexp(kk,tt).*(100/grd.zw(4)).^(-bC(kk,tt)) ;
            POC2000(kk,tt) = POCexp(kk,tt).*(2000/grd.zw(4)).^(-bC(kk,tt)) ;
            if MLD(kk,tt) > grd.zw(4);
                POCmld(kk,tt) = POCexp(kk,tt).*(MLD(kk,tt)/grd.zw(4)).^(-bC(kk,tt)) ;
                sDOCmld(kk,tt) = sDOCexp(kk,tt).*(MLD(kk,tt)/grd.zw(4)).^(-bC(kk,tt)) ;
                TOCmld(kk,tt) = POCmld(kk,tt) + sDOCmld(kk,tt) ;
            else
                POCmld(kk,tt) = POCexp(kk,tt);
                sDOCmld(kk,tt) = sDOCexp(kk,tt);
                TOCmld(kk,tt) = POCexp(kk,tt) + sDOCmld(kk,tt) + lDOCexp(kk,tt);
            end 
        else
            POCadj(kk,tt) = POCexp(kk,tt).*(100/grd.zw(3)).^(-bC(kk,tt)) ;
            sDOCadj(kk,tt) = sDOCexp(kk,tt).*(100/grd.zw(3)).^(-bC(kk,tt)) ;
            POC2000(kk,tt) = POCexp(kk,tt).*(2000/grd.zw(3)).^(-bC(kk,tt)) ;
            if MLD(kk,tt) > grd.zw(4);
                POCmld(kk,tt) = POCexp(kk,tt).*(MLD(kk,tt)/grd.zw(3)).^(-bC(kk,tt)) ;
                sDOCmld(kk,tt) = sDOCexp(kk,tt).*(MLD(kk,tt)/grd.zw(3)).^(-bC(kk,tt)) ;
                TOCmld(kk,tt) = POCmld(kk,tt) + sDOCmld(kk,tt) ;
            else
                POCmld(kk,tt) = POCexp(kk,tt);
                sDOCmld(kk,tt) = sDOCexp(kk,tt);
                TOCmld(kk,tt) = POCexp(kk,tt) + sDOCmld(kk,tt) + lDOCexp(kk,tt);
            end 
        end 
    end
end

% production at the third layer (mmol/year)
CNPP3l = G(:,:,3).*C2P3D(:,:,3).*dVt(:,:,3)*spa;

% production below 100 m (mmol/year)
CNPP3l = CNPP3l*(grd.zw(4)-100)/grd.dzt(3);
POCadj = POCadj - CNPP3l*(1-gamma-sigC);
sDOCadj = sDOCadj - CNPP3l*sigC ;
%
Sum_lDOCexp200 = sum(int_lDOCexp200(iwet)*12*1e-18);
fprintf('labile DOC export out of 200 m is %3.3e \n',Sum_lDOCexp200);
Sum_sDOCexp200 = sum(int_sDOCexp200(iwet)*12*1e-18);
fprintf('semilabile DOC export out of 200 m is %3.3e \n',Sum_sDOCexp200);
Sum_POCexp = nansum(POCexp(:))*12*1e-18;
fprintf('Model POC export out of EU is %3.3e \n',Sum_POCexp);
Sum_POCexp = nansum(POCadj(:))*12*1e-18;
fprintf('Model POC export out of 100 m is %3.3e \n',Sum_POCexp);
Sum_POCexp = nansum(POC2000(:))*12*1e-18;
fprintf('Model POC export out of 1000 m is %3.3e \n',Sum_POCexp);
Sum_POCexp = nansum(POCmld(:))*12*1e-18;
fprintf('Model POC export out of MLD is %3.3e \n',Sum_POCexp);
Sum_lDOCexp = nansum(lDOCexp(:))*12*1e-18;
fprintf('Model labile DOC export out of EU is %3.3e \n',Sum_lDOCexp);
Sum_sDOCexp = nansum(sDOCexp(:))*12*1e-18;
fprintf('Model semi-labile DOC export out of EU is %3.3e \n',Sum_sDOCexp);
Sum_rDOCexp = nansum(rDOCexp(:))*12*1e-18;
fprintf('Model semi-labile DOC export out of EU is %3.3e \n',Sum_rDOCexp);
Sum_sDOCadj = nansum(sDOCadj(:))*12*1e-18;
fprintf('Model semi-labile DOC export out of 100 m is %3.3e \n',Sum_sDOCadj);
Sum_TOCmld = nansum(TOCmld(:))*12*1e-18;
fprintf('Model semi-labile DOC export out of 100 m is %3.3e \n',Sum_TOCmld);

% mmol/m2/year, save to file 
lDOCexp = lDOCexp./dAt(:,:,3);
sDOCexp = sDOCexp./dAt(:,:,3);
sDOCmld = sDOCmld./dAt(:,:,3);
POCexp = POCexp./dAt(:,:,3);
POCmld = POCmld./dAt(:,:,3);
TOCmld = TOCmld./dAt(:,:,3);
POCadj = POCadj./dAt(:,:,3);
POC2000 = POC2000./dAt(:,:,12);

% total carbon export (mmol/m2/year)
TOCexp = POCexp + sDOCexp + lDOCexp ;
% export ratio
eRatio = TOCexp ./ NPP ;
% calculate export ratio based on Laws et al 2011 L&O
tp = NPP*12/365;
T = nanmean(par.Temp(:,:,1:2),3) ;
eRatio_Laws2011 = (0.5857-0.0165*T).*NPP./(51.7+NPP) ;
R = R(:,:,1);
iR1 = find(R(:) == 1);
iR2 = find(R(:) == 2);
iR3 = find(R(:) == 3);
iR4 = find(R(:) == 4);
iR5 = find(R(:) == 5);
iR6 = find(R(:) == 6);
iR7 = find(R(:) == 7);
iR8 = find(R(:) == 8);
iR9 = find(R(:) == 9);
iR10 = find(R(:) == 10);
iR11 = find(R(:) == 11);
iR12 = find(R(:) == 12);

ikp1 = find(~isnan(sum([eRatio(iR1), eRatio_Laws2011(iR1)],2)));
ikp2 = find(~isnan(sum([eRatio(iR2), eRatio_Laws2011(iR2)],2)));
ikp3 = find(~isnan(sum([eRatio(iR3), eRatio_Laws2011(iR3)],2)));
ikp4 = find(~isnan(sum([eRatio(iR4), eRatio_Laws2011(iR4)],2)));
ikp5 = find(~isnan(sum([eRatio(iR5), eRatio_Laws2011(iR5)],2)));
ikp6 = find(~isnan(sum([eRatio(iR6), eRatio_Laws2011(iR6)],2)));
ikp7 = find(~isnan(sum([eRatio(iR7), eRatio_Laws2011(iR7)],2)));
ikp8 = find(~isnan(sum([eRatio(iR8), eRatio_Laws2011(iR8)],2)));
ikp9 = find(~isnan(sum([eRatio(iR9), eRatio_Laws2011(iR9)],2)));
ikp10 = find(~isnan(sum([eRatio(iR10), eRatio_Laws2011(iR10)],2)));
ikp11 = find(~isnan(sum([eRatio(iR11), eRatio_Laws2011(iR11)],2)));
ikp12 = find(~isnan(sum([eRatio(iR12), eRatio_Laws2011(iR12)],2)));
figure
subplot(3,4,1)
plot(eRatio_Laws2011(iR1(ikp1)),eRatio(iR1(ikp1)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,2)
plot(eRatio_Laws2011(iR2(ikp2)),eRatio(iR2(ikp2)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,3) 
plot(eRatio_Laws2011(iR3(ikp3)),eRatio(iR3(ikp3)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,4) 
plot(eRatio_Laws2011(iR4(ikp4)),eRatio(iR4(ikp4)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,5) 
plot(eRatio_Laws2011(iR5(ikp5)),eRatio(iR5(ikp5)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,6) 
plot(eRatio_Laws2011(iR6(ikp6)),eRatio(iR6(ikp6)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,7) 
plot(eRatio_Laws2011(iR7(ikp7)),eRatio(iR7(ikp7)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,8) 
plot(eRatio_Laws2011(iR8(ikp8)),eRatio(iR8(ikp8)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,9) 
plot(eRatio_Laws2011(iR9(ikp9)),eRatio(iR9(ikp9)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,10) 
plot(eRatio_Laws2011(iR10(ikp10)),eRatio(iR10(ikp10)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,11) 
plot(eRatio_Laws2011(iR11(ikp11)),eRatio(iR11(ikp11)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
subplot(3,4,12) 
plot(eRatio_Laws2011(iR12(ikp12)),eRatio(iR12(ikp12)),'b.')
hold on
plot([0 0.8],[0 0.8],'r--','linewidth',2);
hold off
xlim([0 0.8]); ylim([0 0.8])
if expfig 
    exportfig(gcf,'Figs/eRatio_TOC_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
end 
%--------------- compare to ANCP -----------------------------
% TOCexp = smoothit(grd,M3d,TOCexp,3,1e5);
% POCexp = smoothit(grd,M3d,POCexp,3,1e5);
% DOCexp = smoothit(grd,M3d,DOCexp,3,1e5);

Lat_HOT  = 22+45/60 ; Lon_HOT = mod(-158,360);
Lat_BATS = 31+40/60 ; Lon_BATS = mod((-64-10/60),360);
Lat_OSP  = 50+1/60  ; Lon_OSP  = mod((-144-9/60),360);

indx_hot = length(find(grd.xt<Lon_HOT));
indy_hot = length(find(grd.yt<Lat_HOT));

indx_bats = length(find(grd.xt<Lon_BATS));
indy_bats = length(find(grd.yt<Lat_BATS));

indx_osp = length(find(grd.xt<Lon_OSP));
indy_osp = length(find(grd.yt<Lat_OSP));

% find ANCP at specific location and convert unit from
%  mol/m2/year;
TOCexp_HOT = nanmean([TOCexp(indy_hot-1,indx_hot-1), ...
                    TOCexp(indy_hot,    indx_hot-1), ...
                    TOCexp(indy_hot+1,  indx_hot-1), ...
                    TOCexp(indy_hot-1,  indx_hot),   ...
                    TOCexp(indy_hot,    indx_hot),   ...
                    TOCexp(indy_hot+1,  indx_hot),   ...
                    TOCexp(indy_hot-1,  indx_hot+1), ...
                    TOCexp(indy_hot,    indx_hot+1), ...
                    TOCexp(indy_hot+1,  indx_hot+1)])/1000 ;

TOCexp_HOTstd = nanstd([TOCexp(indy_hot-1,indx_hot-1), ...
                    TOCexp(indy_hot,    indx_hot-1), ...
                    TOCexp(indy_hot+1,  indx_hot-1), ...
                    TOCexp(indy_hot-1,  indx_hot),   ...
                    TOCexp(indy_hot,    indx_hot),   ...
                    TOCexp(indy_hot+1,  indx_hot),   ...
                    TOCexp(indy_hot-1,  indx_hot+1), ...
                    TOCexp(indy_hot,    indx_hot+1), ...
                    TOCexp(indy_hot+1,  indx_hot+1)])/1000 ;

TOCexp_BATS = nanmean([TOCexp(indy_bats-1,indx_bats-1), ...
                    TOCexp(indy_bats,    indx_bats-1), ...
                    TOCexp(indy_bats+1,  indx_bats-1), ...
                    TOCexp(indy_bats-1,  indx_bats),   ...
                    TOCexp(indy_bats,    indx_bats),   ...
                    TOCexp(indy_bats+1,  indx_bats),   ...
                    TOCexp(indy_bats-1,  indx_bats+1), ...
                    TOCexp(indy_bats,    indx_bats+1), ...
                    TOCexp(indy_bats+1,  indx_bats+1)])/1000 ;

TOCexp_BATSstd = nanstd([TOCexp(indy_bats-1,indx_bats-1), ...
                    TOCexp(indy_bats,    indx_bats-1), ...
                    TOCexp(indy_bats+1,  indx_bats-1), ...
                    TOCexp(indy_bats-1,  indx_bats),   ...
                    TOCexp(indy_bats,    indx_bats),   ...
                    TOCexp(indy_bats+1,  indx_bats),   ...
                    TOCexp(indy_bats-1,  indx_bats+1), ...
                    TOCexp(indy_bats,    indx_bats+1), ...
                    TOCexp(indy_bats+1,  indx_bats+1)])/1000 ;

TOCexp_OSP  = nanmean([TOCexp(indy_osp-1,indx_osp-1), ...
                    TOCexp(indy_osp,    indx_osp-1), ...
                    TOCexp(indy_osp+1,  indx_osp-1), ...
                    TOCexp(indy_osp-1,  indx_osp),   ...
                    TOCexp(indy_osp,    indx_osp),   ...
                    TOCexp(indy_osp+1,  indx_osp),   ...
                    TOCexp(indy_osp-1,  indx_osp+1), ...
                    TOCexp(indy_osp,    indx_osp+1), ...
                    TOCexp(indy_osp+1,  indx_osp+1)])/1000 ;

TOCexp_OSPstd  = nanstd([TOCexp(indy_osp-1,indx_osp-1), ...
                    TOCexp(indy_osp,    indx_osp-1), ...
                    TOCexp(indy_osp+1,  indx_osp-1), ...
                    TOCexp(indy_osp-1,  indx_osp),   ...
                    TOCexp(indy_osp,    indx_osp),   ...
                    TOCexp(indy_osp+1,  indx_osp),   ...
                    TOCexp(indy_osp-1,  indx_osp+1), ...
                    TOCexp(indy_osp,    indx_osp+1), ...
                    TOCexp(indy_osp+1,  indx_osp+1)])/1000 ;

% print the results to the screen mol/m2/year
fprintf('TOC export at HOT  is %2.2f pm %2.2f mol C m-2 day-1 \n', TOCexp_HOT, TOCexp_HOTstd)
fprintf('TOC export at BATS is %2.2f pm %2.2f mol C m-2 day-1 \n',TOCexp_BATS, TOCexp_BATSstd)
fprintf('TOC export at OSP  is %2.2f pm %2.2f mol C m-2 day-1 \n\n', TOCexp_OSP, TOCexp_OSPstd)

msk_tropical = M3d(:,:,3)*0;
msk_tropical(length(find(grd.yt<-15)):length(find(grd.yt<15)),:) = 1;
iTropical = find(msk_tropical(:)) ;

junk1 = M3d(:,:,3) * 0 ;
junk1(length(find(grd.yt<-30)):length(find(grd.yt<30)),:) = 1;
msk_subtro = junk1 - msk_tropical ;
iSubtro = find(msk_subtro(:)) ;

junk2  = M3d(:,:,3) * 0 ;
junk2(length(find(grd.yt<-45)):length(find(grd.yt<45)),:) = 1;
msk_subtro_subpo = junk2 - junk1 ;
iSubtro_subpo = find(msk_subtro_subpo(:)) ;

junk3 =  M3d(:,:,3) * 0 ;
junk3(length(find(grd.yt<-60)):length(find(grd.yt<60)),:) = 1;
msk_subpolar = junk3 - junk2 ;
iSubpolar = find(msk_subpolar(:)) ;

% Area of the third layer
Area3 = grd.DXT3d(:,:,3) .* grd.DYT3d(:,:,3) ;

% units mol/m2/year;
TOCexp_tropical = nansum(TOCexp(iTropical).*Area3(iTropical)) / nansum(Area3(iTropical)) ;
TOCexp_subtro = nansum(TOCexp(iSubtro) .* Area3(iSubtro)) / nansum(Area3(iSubtro)) ;
TOCexp_subtro_subpo = nansum(TOCexp(iSubtro_subpo).*Area3(iSubtro_subpo)) / nansum(Area3(iSubtro_subpo)) ;
TOCexp_subpolar = nansum(TOCexp(iSubpolar).*Area3(iSubpolar)) / nansum(Area3(iSubpolar)) ;

mean_TOC_tropical = TOCexp_tropical/1000;
mean_TOC_subtro = TOCexp_subtro/1000;
mean_TOC_subtro_subpo = TOCexp_subtro_subpo/1000;
mean_TOC_subpolar = TOCexp_subpolar/1000;

fprintf('TOC export at tropical is %2.2f mol/m2/yr\n', mean_TOC_tropical)
fprintf('TOC export at subtropical is %2.2f mol/m2/yr \n',mean_TOC_subtro)
fprintf('TOC export at subtropical-subpolar is %2.2f mol/m2/yr \n', mean_TOC_subtro_subpo)
fprintf('TOC export at subpolar is %2.2f mol/m2/year \n\n', mean_TOC_subpolar)

% calculate DOC to TOC export ratio for the four biomes.
DOCexp = lDOCexp + sDOCexp;
D2T = DOCexp./TOCexp ;
D2T(isinf(D2T)) = nan ;
D2T_tropical = D2T(iTropical) ;
D2T_subtro = D2T(iSubtro) ;
D2T_subtro_subpo = D2T(iSubtro_subpo) ;
D2T_subpolar = D2T(iSubpolar) ;

mean_D2T_tropical = nansum(D2T_tropical.*Area3(iTropical)) / nansum(Area3(iTropical)) ;
mean_D2T_subtro   = nansum(D2T_subtro.*Area3(iSubtro)) / nansum(Area3(iSubtro)) ;
mean_D2T_subtro_subpo = nansum(D2T_subtro_subpo.*Area3(iSubtro_subpo)) / nansum(Area3(iSubtro_subpo));
mean_D2T_subpolar = nansum(D2T_subpolar.*Area3(iSubpolar)) / nansum(Area3(iSubpolar));

fprintf('tropial zonal mean DOC to TOC export ratio is %2.2f percent\n', ...
        mean_D2T_tropical*100)
fprintf('subtropical zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subtro*100)
fprintf('subtropical subpolar zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subtro_subpo*100)
fprintf('subpolar zonal mean DOC to TOC export ratio is %2.2f percent\n\n', ...
        mean_D2T_subpolar*100)

save mkFigs/CbPM_Cexp_Gamma1to100 POC2000 POCadj POCexp POCmld TOCexp TOCmld lDOCexp sDOCexp sDOCmld