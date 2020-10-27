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

par.optim   = on ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.LoadOpt = off ; % if load optimial par. 
par.pscale  = 0.0 ;
par.cscale  = 0.75 ; % factor to weigh DOC in the objective function

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
%
%-------------load data and set up parameters---------------------
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
        base_name = strcat(VER,'_PCOv2');
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
par.fname = strcat(fname,'.mat') ; 
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ; 

% -------------------update initial guesses --------------
if isfile(par.fname)
    load(par.fname)
end 

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

%----------------- solve the P model ----------------------
p0 = par.p0 ;
[par, P ] = eqPcycle(p0, par)  ;
DIP = M3d+nan ;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
POP = M3d+nan ;  POP(iwet) = P(1+1*nwet:2*nwet) ;
DOP = M3d+nan ;  DOP(iwet) = P(1+2*nwet:3*nwet) ;
pdata.DIP = DIP ; pdata.DOP = DOP ; pdata.POP = POP ;
par.DIP   = DIP(iwet) ;

%  Solve for steady-state carbon distribution
% step forward in time with step dt (yr)
TRdiv = par.TRdiv ;
I     = par.I     ;
Tz    = par.Tz    ;
Na    = 1.773e20  ; %  molar volume of atmosphere
dVt   = par.dVt   ;
vw    = dVt(iwet)./Na ;
% fixed parameters
kappa_p = par.kappa_p ;
kPIC    = par.kappa_p ;
% parameters need to be optimized
sigma = par.sigma ;
d     = par.d     ;
R_Si  = par.R_Si  ;
rR    = par.rR    ;
bC_T  = par.bC_T  ;
bC    = par.bC    ;
kC_T  = par.kC_T  ;
kdC   = par.kdC   ;
alpha = par.alpha ;
beta  = par.beta  ;
cc    = par.cc    ;
dd    = par.dd    ;

kC    = d0(kC_T * Tz + kdC) ;
PO4   = po4obs(iwet)        ;
C2P   = 1./(cc * PO4 + dd)  ;
par.C2P = C2P ;
% particle flux div_rergence [s^-1];
PFDa = buildPFD(par,'PIC') ;
PFDc = buildPFD(par,'POC') ;
par.PFDa = PFDa ;
par.PFDc = PFDc ;
% biological DIC uptake operator
G = uptake_C(par); 
% rain ratio of PIC to POC
vout = mkPIC2P(par) ;
RR   = vout.RR   ;

DIC = data.DIC(iwet) - par.dicant(iwet) ;
POC = data.POC(iwet) ;
DOC = data.DOC(iwet) ;
PIC = data.PIC(iwet) ;
ALK = data.ALK(iwet) ;
pco2atm = par.pco2_air(1) ;  % uatm
C  = [DIC; POC; DOC; PIC; ALK; pco2atm] ;
% air sea gas exchange
par.DIC  = DIC  ;
par.ALK  = ALK  ;
par.pco2atm = pco2atm ;
vout  = Fsea2air(par, 'CO2')   ;
JgDIC = vout.JgDIC ;
G_dic = vout.G_dic ;
G_atm = vout.G_atm ;

PI = [[  I, 0*I, 0*I, 0*I, 0*I]; ...
      [0*I,   I, 0*I, 0*I, 0*I]; ...
      [0*I, 0*I,   I, 0*I, 0*I]; ...
      [0*I, 0*I, 0*I,   I, 0*I]; ...
      [0*I, 0*I, 0*I, 0*I,   I]];

AI = [[PI, sparse(5*nwet,1)]; sparse(1,5*nwet+1)];
AI(end,end) = 1;

PII = [[TRdiv,  0*I,   0*I,  0*I,   0*I]; ...
       [  0*I, PFDc,   0*I,  0*I,   0*I]; ...
       [  0*I,  0*I, TRdiv,  0*I,   0*I]; ...
       [  0*I,  0*I,   0*I, PFDa,   0*I]; ...
       [  0*I,  0*I,   0*I,  0*I, TRdiv]];

AII = [[PII, sparse(5*nwet,1)]; sparse(1,5*nwet+1)];

%%%%%%%%%%%%
N2C   = 16/117;
kappa_g = par.kappa_g ;
ftime = 0.05 ;
dt = ftime*spa ;
l  = 1e6 ;

A = AI + (dt/2)*AII ;
fprintf('factoring a huge matrix...\n');
FA = mfactor(A);
fprintf('Done, start stepping...\n');
B = AI - (dt/2)*AII ;
kk = 1 ; % indices for saving data vector;
sDICbar = par.sDICbar ;
sALKbar = par.sALKbar ;
for t  = 1 : l
    dDICdt = (I+(1-sigma)*RR)*(G*C2P) - kC*DOC - kPIC*PIC - JgDIC + pme*sDICbar; 
    dPOCdt = -(1-sigma)*G*C2P + kappa_p*POC      ; 
    dDOCdt = -sigma*G*C2P + kC*DOC - kappa_p*POC ; 
    dPICdt = -(1-sigma)*RR*(G*C2P) + kPIC*PIC    ; 
    dALKdt = 2*(1-sigma)*RR*(G*C2P) - N2C*G*C2P + N2C*kC*DOC ...
             - 2*kPIC*PIC - kappa_g*(ALK - par.ALKbar) + pme*sALKbar ;
    dATMdt = JgDIC'*vw*1000 ;
    
    dCdt = [dDICdt; dPOCdt; dDOCdt; dPICdt; dALKdt; dATMdt] ;
    C    = mfactor(FA, (B*C - dCdt*dt)) ;
    
    DIC = C(0*nwet+1:1*nwet) ;
    POC = C(1*nwet+1:2*nwet) ;
    DOC = C(2*nwet+1:3*nwet) ;
    PIC = C(3*nwet+1:4*nwet) ;
    ALK = C(4*nwet+1:5*nwet) ;
    pco2atm = C(end) ;
    
    sDICbar = sum(DIC(iwet(isrf)).*dVt(iwet(isrf)))./sum(dVt(iwet(isrf))) ;
    sALKbar = sum(ALK(iwet(isrf)).*dVt(iwet(isrf)))./sum(dVt(iwet(isrf))) ;
    
    par.DIC = DIC ;
    par.ALK = ALK ;
    par.pco2atm = pco2atm ;
    % Air-Sea gas exchange
    vout  = Fsea2air(par, 'CO2');
    G_dic = vout.G_dic ; % flux_dic 
    G_atm = vout.G_atm ; % flux_co2atm
    JgDIC = vout.JgDIC ; % flux mmol/m3/s 
    
    eq1 = (I+(1-sigma)*RR)*(G*C2P) + TRdiv*DIC - kC*DOC - kPIC*PIC ...
          - JgDIC + pme*sDICbar ;
    eq2 = -(1-sigma)*G*C2P + (PFDc+kappa_p*I)*POC        ; 
    eq3 = -sigma*G*C2P + (TRdiv+kC)*DOC - kappa_p*POC    ; 
    eq4 = -(1-sigma)*RR*(G*C2P) + (PFDa+kPIC*I)*PIC      ; 
    eq5 = 2*(1-sigma)*RR*(G*C2P) + TRdiv*ALK - N2C*G*C2P + N2C*kC*DOC ...
          - 2*kPIC*PIC - kappa_g*(ALK - par.ALKbar) + pme*sALKbar ;
    eqa = JgDIC'*vw*1000 ; % umol/mol/s 
    
    F   = [eq1; eq2; eq3; eq4; eq5; eqa] ;
    
    fprintf('.')
    if mod(t,100) == 0
        fprintf('\n')
    end
    
    if mod(t, 500) == 0
        kk = kk + 1 ;
        hist{kk} = [ftime*t; pco2atm] ;
        
        fprintf('current iteration % 3d \n', t)
        fprintf('current error %3.3e \n', norm(F))
        fprintf('Atm CO2 concentration %3.2f \n', pco2atm)
        % yyaxis left
        % plot(t*ftime, norm(F), 'co'); hold on; drawnow 
        % ylabel('Error')
        
        % yyaxis right
        % pco2atm = C(end) ;
        % plot(t*ftime, pco2atm, 'r*'); hold on; drawnow 
        % ylabel('Atm CO_2 (ppm)')
        % xlabel('Elapsed time (yr)')
    end 
    if(norm(F) < 1.0e-9)
        break
    end
end

if par.Omodel == on 
    par.DIC = C(0*nwet+1:1*nwet) ; 
    par.DOC = C(2*nwet+1:3*nwet) ;
    % O2C_T
    if (par.opt_O2C_T == on)
        par.O2C_T = x(par.pindx.O2C_T) ;
    end
    
    % rO2C
    if (par.opt_rO2C == on)
        lrO2C    = x(par.pindx.lrO2C) ;
        par.rO2C = exp(lrO2C)     ;
    end
    % O2P_T
    if (par.opt_O2P_T == on)
        par.O2P_T = x(par.pindx.O2P_T) ;
    end
    
    % rO2P
    if (par.opt_rO2P == on)
        lrO2P    = x(par.pindx.lrO2P) ;
        par.rO2P = exp(lrO2P)     ;
    end
    options.iprint = 1   ; 
    options.atol = 1e-10 ;
    options.rtol = 1e-10 ;
    fprintf('Solving C model ...\n') ;
    
    fprintf('Solving O model ...\n') ;
    X0  = GO ;
    [O,ierr] = nsnew(X0,@(X) O_eqn(X, par),options) ;
end 

function [F, FD] = O_eqn(O2, par)
    on = true; off = false;
    %
    % fixed parameters
    iwet  = par.iwet   ;
    nwet  = par.nwet   ;
    TRdiv = par.TRdiv ;
    I     = speye(nwet)   ;
    PO4   = par.po4obs(iwet) ;
    % variables from C model
    DOC   = par.DOC    ;
    Tz    = par.Tz     ;
    TZ    = par.Tz*1e8 ;
    %
    % tunable parameters;
    O2C_T = par.O2C_T   ; 
    rO2C  = par.rO2C    ;
    kC_T  = par.kC_T    ;
    kdC   = par.kdC     ;
    O2P_T = par.O2P_T   ;
    rO2P  = par.rO2P    ;
    kC    = d0(kC_T * Tz + kdC) ; 
    %
    vout  = mkO2P(par) ;
    O2P   = vout.O2P   ;
    %
    % O2 saturation concentration
    vout  = Fsea2air(par,'O2') ;
    KO2   = vout.KO2   ;
    o2sat = vout.o2sat ;
    % rate of o2 production
    G   = uptake_C(par) ;
    PO2 = G*O2P        ;
    
    % parobolic function for o2 consumption
    R    = 0.5 + 0.5*tanh(O2-10)    ;
    dRdO = 0.5 - 0.5*tanh(O2-10).^2 ;
    
    O2C  = O2C_T*TZ + rO2C ; 
    % rate of o2 utilization
    LO2  = kC*DOC.*O2C.*R  ;
    dLdO = d0(kC*DOC.*O2C.*dRdO) ;
    
    % O2 function
    F = TRdiv*O2 - PO2 + LO2 - KO2*(o2sat-O2) ;
    %
    if (nargout > 1)
        FD = mfactor(TRdiv + dLdO + KO2) ;
    end
end

output_dir = sprintf('/DFS-L/DATA/primeau/weilewang/COP4WWF/PerturbNPP/');
fname = strcat(output_dir,file1) ;
save(fname,'hist','P','C','JgDIC')
fprintf('-------------- end! ---------------\n');
