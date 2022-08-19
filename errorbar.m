clc; clear alul; close all
global iter
iter = 0 ;
on   = true  ;
off  = false ;
format long
% 
GridVer  = 91  ;
operator = 'A' ;

Gtest = off ;
Htest = off ;
par.optim   = on ; 
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
par.opt_beta  = on ;
par.opt_alpha = on ;
% C model parameters
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

% save results 
% ATTENTION: Change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/OneDrive/rDOC/MSK%2d/',GridVer); 
elseif isunix
    output_dir = sprintf('~/rDOC-OP/MSK%2d/',GridVer) ;
end
VER = strcat(output_dir,TRdivVer);

% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_Pv1');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
        base_name = strcat(VER,'_PCv2');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
        base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta');
        % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_MODIS_CbPM_aveTeu_diffSig_O2C_uniEta');
        % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta');
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
par.fxhat = strcat(fname,'_xhat.mat');
load(par.fxhat) ;
load(par.fname) ;
%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par) ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

idip = find(par.po4raw(iwet)>0)  ;
isil = find(par.sio4raw(iwet)>0) ;
idic = find(par.dicraw(iwet)>0)  ;
ialk = find(par.alkraw(iwet)>0)  ;
idoc = find(par.docraw(iwet)>0)  ;
io2  = find(par.o2raw(iwet)>0)   ;

ndip = length(idip) ;
nsil = length(isil) ;
ndic = length(idic) ;
nalk = length(ialk) ;
ndoc = length(idoc) ;
no2  = length(io2)  ;

if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
    sig = 2*xhat.f/ndip ;
    HH  = xhat.fxx/sig  ;
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
    if par.cscale ~= 0
        sig = (2*xhat.f)/(ndip+ndic+ndoc);
        HH  = xhat.fxx/sig  ;
    else
        sig = (2*xhat.f)/(ndip+ndic+nalk);
        HH  = xhat.fxx/sig  ;
    end 
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
    if par.cscale ~= 0
        sig = (2*xhat.f)/(ndip+ndic+nalk+ndoc+no2);
        HH  = xhat.fxx/sig  ;
    else
        sig = (2*xhat.f)/(ndip+ndic+no2);
        HH  = xhat.fxx/sig  ;
    end
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
    if par.cscale ~= 0
        sig = (2*xhat.f)/(ndip+ndic+nalk+ndoc+no2+nsil) ;
        HH  = xhat.fxx/sig  ;
    else
        sig = (2*xhat.f)/(ndip+ndic+no2+nsil) ;
        HH  = xhat.fxx/sig  ;
    end
elseif (par.Cmodel == off & par.Omodel == off & par.Simodel == on)
    sig = (2*xhat.f)/(ndip+nsil);
    HH  = xhat.fxx/sig  ; 
end

pindex = par.pindx ;
error  = sqrt(diag(inv(HH))) ;
nx = 0 ;
n  = numel(fieldnames(par.pindx)) ;
results.xhat   = zeros(n, 1) ;
results.upbar  = zeros(n, 1) ;
results.lowbar = zeros(n, 1) ;

name = cell(n, 1) ;
if exist('xhat') & isfield(xhat,'sigP')
    nx = nx + 1 ; 
    sigP = xhat.sigP ;
    sigP_up = exp(log(sigP)+error(pindex.lsigP)) - sigP;
    sigP_lo = sigP - exp(log(sigP)-error(pindex.lsigP));
    name{nx} = 'sigP' ;
    results.xhat(nx)   = sigP    ;
    results.upbar(nx)  = sigP_up ;
    results.lowbar(nx) = sigP_lo ;
end

if exist('xhat') & isfield(xhat,'Q10P')
    nx = nx + 1 ; 
    Q10P = xhat.Q10P ;
    Q10P_up  = exp(log(Q10P)+error(pindex.lQ10P)) - Q10P;
    Q10P_lo  = Q10P - exp(log(Q10P)-error(pindex.lQ10P));
    name{nx} = "Q10P" ;
    results.xhat(nx)   = Q10P    ;
    results.upbar(nx)  = Q10P_up ;
    results.lowbar(nx) = Q10P_lo ;
end 

if exist('xhat') & isfield(xhat,'kdP')
    nx = nx + 1 ; 
    kdP = 1/xhat.kdP/spa ;
    kdP_up = 1/exp(log(xhat.kdP)+error(pindex.lkdP))/spa - kdP ;
    kdP_lo = kdP - 1/exp(log(xhat.kdP)-error(pindex.lkdP))/spa ;
    name{nx} = "kdP" ;
    results.xhat(nx)   = kdP    ;
    results.upbar(nx)  = kdP_up ;
    results.lowbar(nx) = kdP_lo ;
end 

if exist('xhat') & isfield(xhat,'bP_T')
    nx = nx + 1 ; 
    bP_T = xhat.bP_T ;
    bP_T_up = (bP_T+error(pindex.bP_T)) - bP_T;
    bP_T_lo = bP_T - (bP_T-error(pindex.bP_T));
    name{nx} = "bP_T" ;
    results.xhat(nx)   = bP_T    ;
    results.upbar(nx)  = bP_T_up ;
    results.lowbar(nx) = bP_T_lo ;
end 

if exist('xhat') & isfield(xhat,'bP')
    nx = nx + 1 ; 
    bP  = xhat.bP ;
    bP_up = exp(log(bP)+error(pindex.lbP)) - bP;
    bP_lo = bP - exp(log(bP)-error(pindex.lbP));
    name{nx}     = "bP"  ;
    results.xhat(nx)   = bP    ;
    results.upbar(nx)  = bP_up ;
    results.lowbar(nx) = bP_lo ;
end 

if exist('xhat') & isfield(xhat,'alpha')
    nx = nx + 1 ; 
    alpha = xhat.alpha ;
    alpha_up = exp(log(alpha)+error(pindex.lalpha)) - alpha;
    alpha_lo = alpha - exp(log(alpha)-error(pindex.lalpha));
    name{nx}     = "alpha"  ;
    results.xhat(nx)   = (1/alpha)/spa    ;
    results.upbar(nx)  = (1/alpha_up)/spa ;
    results.lowbar(nx) = (1/alpha_lo)/spa ;
end 

if exist('xhat') & isfield(xhat,'beta')
    nx = nx + 1 ; 
    beta = xhat.beta ;
    beta_up = exp(log(beta)+error(pindex.lbeta)) - beta;
    beta_lo = beta - exp(log(beta)-error(pindex.lbeta));
    name{nx}     = "beta"  ;
    results.xhat(nx)   = beta    ;
    results.upbar(nx)  = beta_up ;
    results.lowbar(nx) = beta_lo ;
end 

% C model parameters
if exist('xhat') & isfield(xhat,'sigC')
    nx = nx + 1 ; 
    sigC = xhat.sigC ;
    sigC_up = exp(log(sigC)+error(pindex.lsigC)) - sigC;
    sigC_lo = sigC - exp(log(sigC)-error(pindex.lsigC));
    name{nx} = 'sigC' ;
    results.xhat(nx)   = sigC    ;
    results.upbar(nx)  = sigC_up ;
    results.lowbar(nx) = sigC_lo ;
end

if exist('xhat') & isfield(xhat,'kru')
    nx = nx + 1 ; 
    kru = 1/xhat.kru/spa ;
    kru_up = 1/exp(log(xhat.kru)+error(pindex.lkru))/spa - kru;
    kru_lo = kru - 1/exp(log(xhat.kru)-error(pindex.lkru))/spa;
    name{nx}     = "kru"  ;
    results.xhat(nx)   = kru    ;
    results.upbar(nx)  = kru_up ;
    results.lowbar(nx) = kru_lo ;
end 
keyboard
if exist('xhat') & isfield(xhat,'krd')
    nx = nx + 1 ; 
    krd = xhat.krd ;
    krd_up = exp(log(krd)+error(pindex.lkrd)) - krd;
    krd_lo = krd - exp(log(krd)-error(pindex.lkrd));
    name{nx}     = "krd"  ;
    results.xhat(nx)   = (1/krd)/spa    ;
    results.upbar(nx)  = (1/krd_up)/spa ;
    results.lowbar(nx) = (1/krd_lo)/spa ;
end 

if exist('xhat') & isfield(xhat,'etau')
    nx = nx + 1 ; 
    etau = xhat.etau ;
    etau_up = exp(log(etau)+error(pindex.letau)) - etau;
    etau_lo = etau - exp(log(etau)-error(pindex.letau));
    name{nx}     = "etau"  ;
    results.xhat(nx)   = etau    ;
    results.upbar(nx)  = etau_up ;
    results.lowbar(nx) = etau_lo ;
end 

if exist('xhat') & isfield(xhat,'etad')
    nx = nx + 1 ; 
    etad = xhat.etad ;
    etad_up = (etad+error(pindex.letad)) - etad;
    etad_lo = etad - (etad-error(pindex.letad));
    name{nx}     = "etad"  ;
    results.xhat(nx)   = etad    ;
    results.upbar(nx)  = etad_up ;
    results.lowbar(nx) = etad_lo ;
end 

if exist('xhat') & isfield(xhat,'bC_T')
    nx = nx + 1 ; 
    bC_T = xhat.bC_T ;
    bC_T_up = (bC_T+error(pindex.bC_T)) - bC_T;
    bC_T_lo = bC_T - (bC_T-error(pindex.bC_T));
    name{nx}     = "bC_T"  ;
    results.xhat(nx)   = bC_T    ;
    results.upbar(nx)  = bC_T_up ;
    results.lowbar(nx) = bC_T_lo ;
end 

if exist('xhat') & isfield(xhat,'bC')
    nx = nx + 1 ; 
    bC = xhat.bC ;
    bC_up = exp(log(bC)+error(pindex.lbC)) - bC;
    bC_lo = bC - exp(log(bC)-error(pindex.lbC));
    name{nx}     = "bC"  ;
    results.xhat(nx)   = bC    ;
    results.upbar(nx)  = bC_up ;
    results.lowbar(nx) = bC_lo ;
end 

if exist('xhat') & isfield(xhat,'d')
    nx = nx + 1 ; 
    d = xhat.d   ;
    d_up = exp(log(d)+error(pindex.ld)) - d;
    d_lo = d - exp(log(d)-error(pindex.ld));
    name{nx}     = "d"  ;
    results.xhat(nx)   = d    ;
    results.upbar(nx)  = d_up ;
    results.lowbar(nx) = d_lo ;
end 

if exist('xhat') & isfield(xhat,'Q10C')
    nx = nx + 1 ; 
    Q10C = xhat.Q10C;
    Q10C_up = exp(log(Q10C)+error(pindex.lQ10C)) - Q10C;
    Q10C_lo = Q10C - exp(log(Q10C)-error(pindex.lQ10C));
    name{nx}     = "Q10C"  ;
    results.xhat(nx)   = Q10C    ;
    results.upbar(nx)  = Q10C_up ;
    results.lowbar(nx) = Q10C_lo ;
end 

if exist('xhat') & isfield(xhat,'kdC')
    nx  = nx + 1 ; 
    kdC = xhat.kdC ;
    kdC_up = exp(log(kdC)+error(pindex.lkdC)) - kdC;
    kdC_lo = kdC - exp(log(kdC)-error(pindex.lkdC));
    name{nx}     = "kdC"  ;
    results.xhat(nx)   = 1/kdC/spa    ;
    results.upbar(nx)  = 1/kdC_up/spa ;
    results.lowbar(nx) = 1/kdC_lo/spa ;
end 

if exist('xhat') & isfield(xhat,'R_Si')
    nx = nx + 1 ; 
    R_Si = xhat.R_Si  ;
    R_Si_up = error(pindex.R_Si) ;
    R_Si_lo = error(pindex.R_Si) ;
    name{nx}     = "R_Si"  ;
    results.xhat(nx)   = R_Si    ;
    results.upbar(nx)  = R_Si_up ;
    results.lowbar(nx) = R_Si_lo ;
end

if exist('xhat') & isfield(xhat,'rR')
    nx = nx + 1 ; 
    rR = xhat.rR  ;
    rR_up = exp(log(rR)+error(pindex.lrR)) - rR;
    rR_lo = rR - exp(log(rR)-error(pindex.lrR));
    name{nx}     = "rR"  ;
    results.xhat(nx)   = rR    ;
    results.upbar(nx)  = rR_up ;
    results.lowbar(nx) = rR_lo ;
end

if exist('xhat') & isfield(xhat,'cc')
    nx = nx + 1 ; 
    cc = xhat.cc  ;
    cc_up = exp(log(cc)+error(pindex.lcc)) - cc;
    cc_lo = cc - exp(log(cc)-error(pindex.lcc));
    name{nx}     = "cc"  ;
    results.xhat(nx)   = cc    ;
    results.upbar(nx)  = cc_up ;
    results.lowbar(nx) = cc_lo ;
end 
if exist('xhat') & isfield(xhat,'dd')
    nx = nx + 1 ; 
    dd = xhat.dd  ;
    dd_up = exp(log(dd)+error(pindex.ldd)) - dd;
    dd_lo = dd - exp(log(dd)-error(pindex.ldd));
    name{nx}     = "dd"  ;
    results.xhat(nx)   = dd    ;
    results.upbar(nx)  = dd_up ;
    results.lowbar(nx) = dd_lo ;
end 
%
% O model parameters
if exist('xhat') & isfield(xhat,'O2C_T')
    nx = nx + 1 ; 
    O2C_T = xhat.O2C_T ;
    O2C_T_up = (O2C_T+error(pindex.O2C_T)) - O2C_T;
    O2C_T_lo = O2C_T - (O2C_T-error(pindex.O2C_T));
    name{nx} = "O2C_T" ;
    results.xhat(nx)   = O2C_T    ;
    results.upbar(nx)  = O2C_T_up ;
    results.lowbar(nx) = O2C_T_lo ;
end

if exist('xhat') & isfield(xhat,'rO2C')
    nx = nx + 1 ; 
    rO2C = xhat.rO2C ;
    rO2C_up = exp(log(rO2C)+error(pindex.lrO2C)) - rO2C;
    rO2C_lo = rO2C - exp(log(rO2C)-error(pindex.lrO2C));
    name{nx}     = "rO2C"  ;
    results.xhat(nx)   = rO2C    ;
    results.upbar(nx)  = rO2C_up ;
    results.lowbar(nx) = rO2C_lo ;
end 

%
% Si model parameters
if exist('xhat') & isfield(xhat,'dsi')
    nx = nx + 1 ; 
    dsi = xhat.dsi ;
    dsi_up = exp(log(dsi)+error(pindex.ldsi)) - dsi;
    dsi_lo = dsi - exp(log(dsi)-error(pindex.ldsi));
    name{nx}     = "dsi"  ;
    results.xhat(nx)   = dsi    ;
    results.upbar(nx)  = dsi_up ;
    results.lowbar(nx) = dsi_lo ;
end

if exist('xhat') & isfield(xhat,'at')
    nx = nx + 1 ; 
    at = xhat.at   ;
    at_up = exp(log(at)+error(pindex.lat)) - at;
    at_lo = at - exp(log(at)-error(pindex.lat));
    name{nx}     = "at"  ;
    results.xhat(nx)   = at    ;
    results.upbar(nx)  = at_up ;
    results.lowbar(nx) = at_lo ;
end 

if exist('xhat') & isfield(xhat,'bt')
    nx = nx + 1 ; 
    bt = xhat.bt   ;
    bt_up = exp(log(bt)+error(pindex.lbt)) - bt;
    bt_lo = bt - exp(log(bt)-error(pindex.lbt));
    name{nx}     = "bt"  ;
    results.xhat(nx)   = bt    ;
    results.upbar(nx)  = bt_up ;
    results.lowbar(nx) = bt_lo ;
end

if exist('xhat') & isfield(xhat,'aa')
    nx = nx + 1 ; 
    aa = xhat.aa   ;
    aa_up = exp(log(aa)+error(pindex.laa)) - aa;
    aa_lo = aa - exp(log(aa)-error(pindex.laa));
    name{nx}     = "aa"  ;
    results.xhat(nx)   = aa    ;
    results.upbar(nx)  = aa_up ;
    results.lowbar(nx) = aa_lo ;
end

if exist('xhat') & isfield(xhat,'bb')
    nx = nx + 1 ; 
    bb = xhat.bb   ;
    bb_up = exp(log(bb)+error(pindex.lbb)) - bb;
    bb_lo = bb - exp(log(bb)-error(pindex.lbb));
    name{nx}     = "bb"  ;
    results.xhat(nx)   = bb    ;
    results.upbar(nx)  = bb_up ;
    results.lowbar(nx) = bb_lo ;
end 

xhat   = results.xhat       ;
upbar  = results.upbar      ;
lowbar = results.lowbar     ;
name   = string(name) ;
T = table(xhat,upbar,lowbar,'RowNames',name)
