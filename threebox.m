% indices
clc; clear all; close all
on   = true    ; off = false    ;
spd  = 24*60^2 ; spa  = 365*spd ;
% addpath according to opterating system
if ismac 
    addpath('~/Dropbox/myfunc'     )
    addpath('~/Documents/DATA/'    )
    addpath('~/Documents/DATA/OCIM')
elseif isunix
    addpath('/DFS-L/DATA/primeau/weilewang/my_func/'  )
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/'     )
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
end
format long
%
Cmodel  = on ; 
Omodel  = off ; 
Simodel = off ;
LoadOpt = on  ;
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

fprintf('Transport version: % s \n', TRdivVer)
if Cmodel == on
    fprintf('---- C model is on ---- \n')
end
if Omodel == on
    fprintf('---- O model is on ---- \n')
end 
if Simodel == on
    fprintf('---- Si model is on ---- \n')
end 
fprintf('\n')
% P model parameters
par.optim     = off    ; 
par.Cmodel    = Cmodel ;
par.Omodel    = Omodel ;
par.Simodel   = Simodel;
par.LoadOpt   = LoadOpt;
%
par.opt_sigma = off ; 
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
par.opt_RR    = on ; 
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
%
% save results 
% ATTENTION: please change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer); 
elseif isunix
    output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
                        'MSK%2d/'],GridVer);
end
VER = strcat(output_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if (Cmodel == off & Omodel == off & Simodel == off)
    fname = strcat(VER,'_P');
elseif (Cmodel == on & Omodel == off & Simodel == off)
    fname = strcat(VER,'_PC');
elseif (Cmodel == on & Omodel == on & Simodel == off)
    fname = strcat(VER,'_PCO');
elseif (Cmodel == on & Omodel == off & Simodel == on)
    fname = strcat(VER,'_PCSi');
elseif (Cmodel == on & Omodel == on & Simodel == on)
    fname = strcat(VER,'_PCOSi');
end
pfname    = strcat(fname,'_pert.mat');
par.fname = fname ; 
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ; 
%
if GridVer == 90
    load transport_v4.mat grid M3d TR
    load M3d90x180x24v2.mat MSKS 
    load Sobs_90x180x24.mat
    load tempobs_90x180x24.mat
    load po4obs_90x180x24.mat       % WOA PO4 observation
    load Siobs_90x180x24.mat Siobs
    %
    load DICant_90x180x24.mat
    load GLODAPv2_90x180x24raw.mat
    load splco2_mod_monthly.mat     % monthly CO2 data
    load co2syspar90.mat co2syspar
    load cbpm_npp_annual_90x180.mat
    load kw660_90x180.mat
    %
    grd = grid ;

elseif GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    load M3d91x180x24.mat MSKS 
    load Sobs_91x180x24.mat
    load po4obs_91x180x24.mat % WOA PO4 observation
    load tempobs_91x180x24.mat
    load Siobs_91x180x24.mat Siobs
    %
    load DICant_91x180x24.mat
    load GLODAPv2_91x180x24raw.mat
    load splco2_mod_monthly % monthly CO2 data
    load co2syspar91.mat co2syspar
    load cbpm_npp_annual_91x180.mat
    load kw660_91x180.mat
    %
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
end

load(fname)
load(fxhat) 
% get rid of arctice o2 observations
ARC  = MSKS.ARC;
iarc = find(ARC(:)) ;
o2raw(iarc)   = nan ;
dicraw(iarc)  = nan ;
po4raw(iarc)  = nan ;
sio4raw(iarc) = nan ; 
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
dAt  = grd.DXT3d.*grd.DYT3d ;
dVt  = dAt.*grd.DZT3d ;
%
[par.kw,par.P] = kw(M3d,grd);
par.Salt  = Sobs    ;
par.Temp  = tempobs ;
par.dVt   = dVt     ;
par.Kw660 = Kw660   ;
par.p4    = p4      ;
par.c2p   = 110     ;
par.M3d   = M3d     ;
par.iwet  = iwet    ;
par.nwet  = nwet    ;
par.TRdiv = -TR     ;
par.grd   = grd     ;
par.I     = speye(nwet)  ;
par.rho   = 1024.5       ; % seawater density;
permil    = par.rho*1e-3 ; % from umol/kg to mmol/m3;

par.po4obs    = po4obs  ;
par.human_co2 = DICant*permil;

% transiant CO2 concentraion;
par.year      = splco2_mod(:,1) ;
par.pco2_air  = splco2_mod(:,2) ;
par.co2syspar = co2syspar       ;

% load optimal parameters from a file
% or set them to default values 
par = SetPara(par) ;
%
% pack adjustable parameters in an array and
% assign them corresponding indices.
[p0, par] = PackPar(par) ;
%
PrintPara(p0, par)     ;
% PME part;
[modT,modS] = PME(par) ;
par.modS    = modS     ;
par.modT    = modT     ;
Tz0   = (modT(iwet)-mean(modT(iwet)))/std(modT(iwet)) ;
modT1 = modT + 2       ;
Tz1   = (modT1(iwet)-mean(modT(iwet)))/std(modT(iwet)) ;
Tz3d  = M3d + nan      ;
Tz3d(iwet)  = Tz1      ;
par.aveT    = nanmean(Tz3d(:,:,1:3),3) ;
par.Tz      = Tz1*1e-8 ;


%%%%%%% prepare NPP for the model %%%%%%%%
par.nzo   = 2 ;
par.p2c   = 0.006+0.0069*po4obs ;
inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan) = 0 ;

par.npp    = npp/(12*spd) ;
par.npp1   = (0.5*par.npp./grd.dzt(1)).*par.p2c(:,:,1) ; 
par.npp2   = (0.5*par.npp./grd.dzt(2)).*par.p2c(:,:,2) ; 
par.Lambda = M3d*0 ;
par.Lambda(:,:,1) = 1./(1e-6+po4obs(:,:,1)) ;
par.Lambda(:,:,2) = 1./(1e-6+po4obs(:,:,2)) ;

par.Lambda(:,:,3:end) = 0 ;

% build part of the biological DIP uptake operator
Lambda     = par.Lambda;
LAM        = 0*M3d;
LAM(:,:,1) = (par.npp1.^par.beta).*Lambda(:,:,1);
LAM(:,:,2) = (par.npp2.^par.beta).*Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
par.L      = L;

DIP = data.DIP ;
POP = data.POP ;
DOP = data.DOP ;

par.DIPbar = nansum(po4obs(iwet).*dVt(iwet))/nansum(dVt(iwet)) ;

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

%%% Solve for steady-state carbon distribution
options.iprint = 1 ; 
options.atol = 1e-10 ;
options.rtol = 1e-10 ;
fprintf('Solving C model ...\n') ;

% initial guess 
DIC = data.DIC(iwet) ;  POC = data.POC(iwet) ;
DOC = data.DOC(iwet) ;  PIC = data.CaC(iwet) ;
pco2atm = par.pco2_air(1) ;  % uatm

C  = [DIC; POC; DOC; PIC; pco2atm] ;
% C  = [DIC; POC; DOC; PIC] ;

X0  = C ;
[C,ierr] = nsnew(X0,@(X) C_eqn(X, par),options) ;
if (ierr ~= 0)
    fprintf('eqCcycle did not converge.\n') ;
end

%%%%
function [F,FD] = C_eqn(X, par)
    grd   = par.grd   ;
    M3d   = par.M3d   ;
    TRdiv = par.TRdiv ;
    iwet  = par.iwet  ;
    nwet  = par.nwet  ;
    I     = par.I     ;
    % volumes and areas from T99
    Na  = 1.773e20/1000 ; % molar volume of atmosphere

    Tz  = par.Tz ;
    DIC = X(0*nwet+1:1*nwet) ; 
    POC = X(1*nwet+1:2*nwet) ;
    DOC = X(2*nwet+1:3*nwet) ;
    PIC = X(3*nwet+1:4*nwet) ;
    par.pco2atm = X(end) ;
    %
    PO4 = par.po4obs(iwet);
    %
    % fixed parameters
    kappa_p = par.kappa_p ;
    % parameters need to be optimized
    sigma = par.sigma ;
    d     = par.d     ;
    RR    = par.RR    ;
    bC_T  = par.bC_T  ;
    bC    = par.bC    ;
    kC_T  = par.kC_T  ;
    kdC   = par.kdC   ;
    alpha = par.alpha ;
    beta  = par.beta  ;
    cc    = par.cc    ;
    dd    = par.dd    ;
    dVt   = par.dVt   ;
    vw    = dVt(iwet)./Na ;
    
    kC    = d0(kC_T * Tz + kdC) ;
    C2P   = 1./(cc*PO4 + dd) ;
    par.C2P = C2P ;
    % particle flux div_rergence [s^-1];
    PFDa = buildPFD(par,'PIC') ;
    PFDc = buildPFD(par,'POC') ;
    par.PFDa = PFDa ;
    par.PFDc = PFDc ;
    par.DIC  = DIC  ;
    % Air-Sea gas exchange
    vout  = Fsea2air(par, 'CO2');
    KG    = vout.KG    ; % flux_dic 
    KA    = vout.KA    ; % flux_co2atm
    JgDIC = vout.JgDIC ; % flux mmol/m3/s 
    % biological DIC uptake operator
    G = uptake_C(par); par.G = G;

    eq1 =     (1+RR)*G*C2P + TRdiv*DIC - kC*DOC - kappa_p*PIC - JgDIC;
    eq2 = -(1-sigma)*G*C2P + (PFDc+kappa_p*I)*POC;
    eq3 =     -sigma*G*C2P + (TRdiv+kC)*DOC - kappa_p*POC;
    eq4 =        -RR*G*C2P + (PFDa+kappa_p*I)*PIC;
    eqa =         JgDIC'*vw ; % umol/mol/s 

    F   = [eq1; eq2; eq3; eq4; eqa];
    
    if nargout > 1
        % % column 1 dFdDIC
        Jc{1,1} = TRdiv - KG ; 
        Jc{2,1} = 0*I ;
        Jc{3,1} = 0*I ; 
        Jc{4,1} = 0*I ;
        
        % % column 2 dFdPOC
        Jc{1,2} =  0*I ;
        Jc{2,2} =  PFDc + kappa_p*I ;
        Jc{3,2} = -kappa_p*I ;
        Jc{4,2} =  0*I ;
        
        % % column 3 dFdDOC
        Jc{1,3} = -kC ;
        Jc{2,3} = 0*I ;
        Jc{3,3} = TRdiv + kC ;
        Jc{4,3} = 0*I ;
        
        % % column 4 dFdPIC
        Jc{1,4} = -kappa_p*I ;
        Jc{2,4} = 0*I ;
        Jc{3,4} = 0*I ;
        Jc{4,4} = PFDa + kappa_p*I ;

        % Jacobian without Atmosphere box
        PIII = cell2mat(Jc) ;

        % % column 1 dFdDIC
        Jc{5,1} = KG*vw ;
        % % column 5 dFdatmco2
        Jc{1,5} = -KA ;
        Jc{5,5} =  KA'*vw ;
        
        AIII = [[PIII, sparse(4*nwet,1)]; sparse(1,4*nwet+1)];
        AIII(5,1:nwet) = Jc{5,1} ;
        AIII(1:nwet,5) = Jc{1,5} ;
        AIII(5,5)      = Jc{5,5} ;
        
        FD  = mfactor(AIII) ;
        
    end 
end 



% step forward in time with step dt (yr)
% dt = 1*spa ;

% total carbon inventory
% CT = 2.87805e18 ;

% DIC = data.DIC(iwet) ;  POC = data.POC(iwet) ;
% DOC = data.DOC(iwet) ;  PIC = data.CaC(iwet) ;
% pco2atm = par.pco2_air(1) ;  % uatm
% C(:,1)  = [DIC; POC; DOC; PIC; pco2atm] ;

% biological DIC uptake operator
% G = uptake_C(par); par.G = G;
% kappa_p  = par.kappa_p ;
% sigma    = par.sigma ;
% kC_T     = par.kC_T  ;
% kdC      = par.kdC   ;
% cc       = par.cc    ;
% dd       = par.dd    ;
% RR       = par.RR    ;
% PO4      = po4obs(iwet) ;
% kC    = d0(kC_T * par.Tz + kdC) ;
% C2P   = 1./(cc*PO4 + dd) ;
% par.C2P = C2P ;
% particle flux div_rergence [s^-1];
% PFDa = buildPFD(par,'PIC') ;
% PFDc = buildPFD(par,'POC') ;
% par.PFDa = PFDa ;
% par.PFDc = PFDc ;
% par.DIC  = DIC  ;
% par.pco2atm = pco2atm ;
% TRdiv = par.TRdiv  ;
% I     = speye(nwet);
% vout  = Fsea2air(par, 'CO2')   ;
% JgDIC = vout.JgDIC ;


% F   = [eq1; eq2; eq3; eq4; eqa];
% construct the LHS matrix for the offline model
% disp('Preparing LHS and RHS matrix:')
% colum 1 dFdDIC
% number of time-steps
% l  = 1e5 ;

% for t  = 1 : l
% vout  = Fsea2air(par, 'CO2')   ;
% JgDIC = vout.JgDIC ;
% KG    = vout.KG    ;
% KA    = vout.KA    ;

% column 1 dFdDIC
% Jc{1,1} = TRdiv - KG ; 
% Jc{2,1} = 0*I ;
% Jc{3,1} = 0*I ; 
% Jc{4,1} = 0*I ;
% Jc{5,1} = diag(KG).*dVt(iwet)*1000./Na ;

% % column 2 dFdPOC
% Jc{1,2} =  0*I ;
% Jc{2,2} =  PFDc + kappa_p*I ;
% Jc{3,2} = -kappa_p*I ;
% Jc{4,2} =  0*I ;
% Jc{5,2} =    0 ;

% % column 3 dFdDOC
% Jc{1,3} = -kC ;
% Jc{2,3} = 0*I ;
% Jc{3,3} = TRdiv + kC ;
% Jc{4,3} = 0*I ;
% Jc{5,3} =   0 ;

% % column 4 dFdPIC
% Jc{1,4} = -kappa_p*I ;
% Jc{2,4} = 0*I ;
% Jc{3,4} = 0*I ;
% Jc{4,4} =  PFD + kappa_p*I ;
% Jc{5,4} =   0 ;

% % column 5 dFdatmco2
% Jc{1,5} = -KA ;
% Jc{2,5} =   0 ;
% Jc{3,5} =   0 ;
% Jc{4,5} =   0 ;
% Jc{5,5} =  sum(KA.*dVt(iwet)*1000)/Na ;

% PI = [[  I, 0*I, 0*I, 0*I]; ...
% [0*I,   I, 0*I, 0*I]; ...
% [0*I, 0*I,   I, 0*I]; ...
% [0*I, 0*I, 0*I,   I]];

% AI = [[PI, sparse(4*nwet,1)]; sparse(1,4*nwet+1)];
% AI(end,end) = 1;

% PII = [[TRdiv,  0*I,   0*I,  0*I]; ...
% [  0*I, PFDc,   0*I,  0*I]; ...
% [  0*I,  0*I, TRdiv,  0*I]; ...
% [  0*I,  0*I,   0*I, PFDa]];

% AII = [[PII, sparse(4*nwet,1)]; sparse(1,4*nwet+1)];

% PIII = [[Jc{1,1},Jc{1,2},Jc{1,3},Jc{1,4}]; ...
% [Jc{2,1},Jc{2,2},Jc{2,3},Jc{2,4}]; ...
% [Jc{3,1},Jc{3,2},Jc{3,3},Jc{3,4}]; ...
% [Jc{4,1},Jc{4,2},Jc{4,3},Jc{4,4}]];

% AIII = [[PIII, sparse(4*nwet,1)]; sparse(1,4*nwet+1)];
% AIII(5,1:nwet) = Jc{5,1} ;
% AIII(1:nwet,5) = Jc{1,5} ;
% AIII(5,5) = Jc{5,5} ;

% A = AI - (dt/2)*(AII + AIII) ;

% B = AI + (dt/2)*(AII + AIII) ;
% keyboard
% tic
% fprintf('factoring a huge matrix...\n');
% FA = mfactor(A);
% fprintf('done! \nstart time stepping....\n');
% toc

% C  = mfactor(FA, B*C) ;
% Jac = AII + AIII ;
% FD  = mfactor(PIII) ;

% eq1 =     (1+RR)*G*C2P + TRdiv*DIC - kC*DOC - kappa_p*PIC - JgDIC;
% eq2 = -(1-sigma)*G*C2P + (PFDc+kappa_p*I)*POC;
% eq3 =     -sigma*G*C2P + (TRdiv+kC)*DOC - kappa_p*POC;
% eq4 =        -RR*G*C2P + (PFDa+kappa_p*I)*PIC;
% eqa =   sum(JgDIC.*dVt(iwet)*1000)/Na; % umol/mol 

% F   = [eq1; eq2; eq3; eq4; eqa] ;

% dC = -mfactor(FD, F) ;
% C = C + dC ;

% DIC = C(0*nwet+1:1*nwet) ;
% POC = C(1*nwet+1:2*nwet) ;
% DOC = C(2*nwet+1:3*nwet) ;
% PIC = C(3*nwet+1:4*nwet) ;
% par.pco2atm = C(end) ;

% fprintf('current error %3.3e \n', norm(F))
% plot(t, norm(F), 'ro'); hold on; drawnow 

% if(norm(F) < 1.5e-11)
% break
% end
% end
% C = C(:,1:t) ;

% figure(1)
% pco2atm = C(1,:)*1e6;
% time = (1:t)*dt;
% plot(time,pco2atm); 
% title('Atmospheric pCO2 (ppm)')
