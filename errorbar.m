clc; clear all; close all
global iter
iter = 0 ;
on   = true  ;
off  = false ;
format long
%
GridVer  = 90  ;
operator = 'A' ;

Gtest = off ;
Htest = off ;
par.optim   = on ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.LoadOpt = on ; % if load optimial par.
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigma = on ;
par.opt_kP_T  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ;
par.opt_bP    = on ;
par.opt_alpha = on ;
par.opt_beta  = on ;
% C model parameters
par.opt_bC_T  = on ;
par.opt_bC    = on ;
par.opt_d     = on ;
par.opt_kC_T  = on ;
par.opt_kdC   = on ;
par.opt_R_Si  = on ;
par.opt_rR    = on ;
par.opt_cc    = off ; %always off if cell model on
par.opt_dd    = off ; %always off if cell model on
% O model parameters
par.opt_O2C_T = off ;
par.opt_rO2C  = on ;
par.opt_O2P_T = off ;
par.opt_rO2P  = on ;
% Si model parameters
par.opt_dsi   = off  ;
par.opt_at    = off ;
par.opt_bt    = off  ;
par.opt_aa    = off  ;
par.opt_bb    = off  ;
%Trait Model parameters
par.opt_Q10Photo     = on ;
par.opt_fStorage     = on;
par.opt_PLip_PCutoff = on;
par.opt_PLip_scale   = off;
par.opt_PStor_rCutoff = on;
par.opt_PStor_scale  = off;
par.opt_alphaS       = on;
par.opt_fRibE 	     = on;
par.opt_kST0 	     = off;
%
%-------------load data and set up parameters---------------------
SetUp ;

% save results
% ATTENTION: Change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer);
elseif isunix
	output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/', GridVer);
	% output_dir = sprintf('/DFS-L/DATA/primeau/weilewang/Cexp/');
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
    % 'MSK%2d/PME4DICALK/'],GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
                        % 'MSK%2d/'],GridVer);
end
VER = strcat(output_dir,TRdivVer);
catDOC = sprintf('_DOC%0.2g_DOP%0.2g',par.cscale,par.pscale); % used to add scale factors to file names
% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off & par.Cellmodel == off)
        fname = strcat(VER,'_P');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off & par.Cellmodel == off)
        base_name = strcat(VER,'_PC');
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == off)
        base_name = strcat(VER,'_PCO');
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on & par.Cellmodel == off)
        base_name = strcat(VER,'_PCSi');
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == off)
        base_name = strcat(VER,'_PCOSi');
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == off & par.Omodel == off & par.Simodel == off & par.Cellmodel == on) % cell model does nothing if C model is not on, so this case =Ponly
        base_name = strcat(VER,'_PCell');
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off & par.Cellmodel == on)
        base_name = strcat(VER,'_PCCellv3b');
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == on)
		base_name = strcat(VER,'_PCOCell');
		fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == on)
        base_name = strcat(VER,'_PCOSiCell');
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
idoc = find(par.docraw(iwet)>0)  ;
io2  = find(par.o2raw(iwet)>0)   ;

ndip = length(idip) ;
nsil = length(isil) ;
ndic = length(idic) ;
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
        sig = (2*xhat.f)/(ndip+ndic);
        HH  = xhat.fxx/sig  ;
    end
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
    sig = (2*xhat.f)/(ndip+ndic+no2);
    HH  = xhat.fxx/sig  ;
elseif (par.Cmodel == off & par.Omodel == off & par.Simodel == on)
    sig = (2*xhat.f)/(ndip+nsil);
    HH  = xhat.fxx/sig  ;
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
    sig = (2*xhat.f)/(ndip+ndic+no2+nsil) ;
    HH  = xhat.fxx/sig  ;
end

pindex = par.pindx ;
error  = sqrt(diag(inv(HH))) ;
nx = 0 ;
n  = numel(fieldnames(par.pindx)) ;
%%%
n= 14; %temporary value. change back once cell model parameters are added
%%%%
R.xhat   = zeros(n, 1) ;
R.upbar  = zeros(n, 1) ;
R.lowbar = zeros(n, 1) ;
name = cell(n, 1) ;
if exist('xhat') & isfield(xhat,'sigma')
    nx = nx + 1 ;
    sigma = xhat.sigma ;
    sigma_up = exp(log(sigma)+error(pindex.lsigma)) - sigma;
    sigma_lo = sigma - exp(log(sigma)-error(pindex.lsigma));
    name{nx} = 'sigma' ;
    R.xhat(nx)   = sigma    ;
    R.upbar(nx)  = sigma_up ;
    R.lowbar(nx) = sigma_lo ;
end

if exist('xhat') & isfield(xhat,'kP_T')
    nx = nx + 1 ;
    kP_T = xhat.kP_T ;
    kP_T_up  = (kP_T+error(pindex.kP_T)) - kP_T;
    kP_T_lo  = kP_T - (kP_T-error(pindex.kP_T));
    name{nx} = 'kP_T' ;
    R.xhat(nx)   = kP_T    ;
    R.upbar(nx)  = kP_T_up ;
    R.lowbar(nx) = kP_T_lo ;
end

if exist('xhat') & isfield(xhat,'kdP')
    nx = nx + 1 ;
    kdP = xhat.kdP ;
    kdP_up = exp(log(kdP)+error(pindex.lkdP)) - kdP ;
    kdP_lo = kdP - exp(log(kdP)-error(pindex.lkdP)) ;
    name{nx} = 'kdP' ;
    R.xhat(nx)   = kdP    ;
    R.upbar(nx)  = kdP_up ;
    R.lowbar(nx) = kdP_lo ;
end

if exist('xhat') & isfield(xhat,'bP_T')
    nx = nx + 1 ;
    bP_T = xhat.bP_T ;
    bP_T_up = (bP_T+error(pindex.bP_T)) - bP_T;
    bP_T_lo = bP_T - (bP_T-error(pindex.bP_T));
    name{nx} = 'bP_T' ;
    R.xhat(nx)   = bP_T    ;
    R.upbar(nx)  = bP_T_up ;
    R.lowbar(nx) = bP_T_lo ;
end

if exist('xhat') & isfield(xhat,'bP')
    nx = nx + 1 ;
    bP  = xhat.bP ;
    bP_up = exp(log(bP)+error(pindex.lbP)) - bP;
    bP_lo = bP - exp(log(bP)-error(pindex.lbP));
    name{nx}     = 'bP'  ;
    R.xhat(nx)   = bP    ;
    R.upbar(nx)  = bP_up ;
    R.lowbar(nx) = bP_lo ;
end

if exist('xhat') & isfield(xhat,'alpha')
    nx = nx + 1 ;
    alpha = xhat.alpha ;
    alpha_up = exp(log(alpha)+error(pindex.lalpha)) - alpha;
    alpha_lo = alpha - exp(log(alpha)-error(pindex.lalpha));
    name{nx}     = 'alpha'  ;
    R.xhat(nx)   = alpha    ;
    R.upbar(nx)  = alpha_up ;
    R.lowbar(nx) = alpha_lo ;
end

if exist('xhat') & isfield(xhat,'beta')
    nx = nx + 1 ;
    beta = xhat.beta ;
    beta_up = exp(log(beta)+error(pindex.lbeta)) - beta;
    beta_lo = beta - exp(log(beta)-error(pindex.lbeta));
    name{nx}     = 'beta'  ;
    R.xhat(nx)   = beta    ;
    R.upbar(nx)  = beta_up ;
    R.lowbar(nx) = beta_lo ;
end

% C model parameters
if exist('xhat') & isfield(xhat,'bC_T')
    nx = nx + 1 ;
    bC_T = xhat.bC_T ;
    bC_T_up = (bC_T+error(pindex.bC_T)) - bC_T;
    bC_T_lo = bC_T - (bC_T-error(pindex.bC_T));
    name{nx}     = 'bC_T'  ;
    R.xhat(nx)   = bC_T    ;
    R.upbar(nx)  = bC_T_up ;
    R.lowbar(nx) = bC_T_lo ;
end

if exist('xhat') & isfield(xhat,'bC')
    nx = nx + 1 ;
    bC = xhat.bC ;
    bC_up = exp(log(bC)+error(pindex.lbC)) - bC;
    bC_lo = bC - exp(log(bC)-error(pindex.lbC));
    name{nx}     = 'bC'  ;
    R.xhat(nx)   = bC    ;
    R.upbar(nx)  = bC_up ;
    R.lowbar(nx) = bC_lo ;
end

if exist('xhat') & isfield(xhat,'d')
    nx = nx + 1 ;
    d = xhat.d   ;
    d_up = exp(log(d)+error(pindex.ld)) - d;
    d_lo = d - exp(log(d)-error(pindex.ld));
    name{nx}     = 'd'  ;
    R.xhat(nx)   = d    ;
    R.upbar(nx)  = d_up ;
    R.lowbar(nx) = d_lo ;
end

if exist('xhat') & isfield(xhat,'kC_T')
    nx = nx + 1 ;
    kC_T = xhat.kC_T;
    kC_T_up = (kC_T+error(pindex.kC_T)) - kC_T;
    kC_T_lo = kC_T - (kC_T-error(pindex.kC_T));
    name{nx}     = 'kC_T'  ;
    R.xhat(nx)   = kC_T    ;
    R.upbar(nx)  = kC_T_up ;
    R.lowbar(nx) = kC_T_lo ;
end

if exist('xhat') & isfield(xhat,'kdC')
    nx  = nx + 1 ;
    kdC = xhat.kdC ;
    kdC_up = exp(log(kdC)+error(pindex.lkdC)) - kdC;
    kdC_lo = kdC - exp(log(kdC)-error(pindex.lkdC));
    name{nx}     = 'kdC'  ;
    R.xhat(nx)   = kdC    ;
    R.upbar(nx)  = kdC_up ;
    R.lowbar(nx) = kdC_lo ;
end

if exist('xhat') & isfield(xhat,'R_Si')
    nx = nx + 1 ;
    R_Si = xhat.R_Si  ;
    R_Si_up = error(pindex.R_Si) ;
    R_Si_lo = error(pindex.R_Si) ;
    name{nx}     = 'R_Si'  ;
    R.xhat(nx)   = R_Si    ;
    R.upbar(nx)  = R_Si_up ;
    R.lowbar(nx) = R_Si_lo ;
end

if exist('xhat') & isfield(xhat,'rR')
    nx = nx + 1 ;
    rR = xhat.rR  ;
    rR_up = exp(log(rR)+error(pindex.lrR)) - rR;
    rR_lo = rR - exp(log(rR)-error(pindex.lrR));
    name{nx}     = 'rR'  ;
    R.xhat(nx)   = rR    ;
    R.upbar(nx)  = rR_up ;
    R.lowbar(nx) = rR_lo ;
end

if exist('xhat') & isfield(xhat,'cc')
    nx = nx + 1 ;
    cc = xhat.cc  ;
    cc_up = exp(log(cc)+error(pindex.lcc)) - cc;
    cc_lo = cc - exp(log(cc)-error(pindex.lcc));
    name{nx}     = "cc"  ;
    R.xhat(nx)   = cc    ;
    R.upbar(nx)  = cc_up ;
    R.lowbar(nx) = cc_lo ;
end
if exist('xhat') & isfield(xhat,'dd')
    nx = nx + 1 ;
    dd = xhat.dd  ;
    dd_up = exp(log(dd)+error(pindex.ldd)) - dd;
    dd_lo = dd - exp(log(dd)-error(pindex.ldd));
    name{nx}     = "dd"  ;
    R.xhat(nx)   = dd    ;
    R.upbar(nx)  = dd_up ;
    R.lowbar(nx) = dd_lo ;
end
%
% O model parameters
if exist('xhat') & isfield(xhat,'O2C_T')
    nx = nx + 1 ;
    O2C_T = xhat.O2C_T ;
    O2C_T_up = (O2C_T+error(pindex.O2C_T)) - O2C_T;
    O2C_T_lo = O2C_T - (O2C_T-error(pindex.O2C_T));
    name{nx} = "O2C_T" ;
    R.xhat(nx)   = O2C_T    ;
    R.upbar(nx)  = O2C_T_up ;
    R.lowbar(nx) = O2C_T_lo ;
end

if exist('xhat') & isfield(xhat,'rO2C')
    nx = nx + 1 ;
    rO2C = xhat.rO2C ;
    rO2C_up = exp(log(rO2C)+error(pindex.lrO2C)) - rO2C;
    rO2C_lo = rO2C - exp(log(rO2C)-error(pindex.lrO2C));
    name{nx}     = "rO2C"  ;
    R.xhat(nx)   = rO2C    ;
    R.upbar(nx)  = rO2C_up ;
    R.lowbar(nx) = rO2C_lo ;
end

if exist('xhat') & isfield(xhat,'O2P_T')
    nx = nx + 1 ;
    O2P_T = xhat.O2P_T ;
    O2P_T_up = (O2P_T+error(pindex.O2P_T)) - O2P_T;
    O2P_T_lo = O2P_T - (O2P_T-error(pindex.O2P_T));
    name{nx} = "O2P_T" ;
    R.xhat(nx)   = O2P_T    ;
    R.upbar(nx)  = O2P_T_up ;
    R.lowbar(nx) = O2P_T_lo ;
end

if exist('xhat') & isfield(xhat,'rO2P')
    nx = nx + 1 ;
    rO2P = xhat.rO2P ;
    rO2P_up = exp(log(rO2P)+error(pindex.lrO2P)) - rO2P;
    rO2P_lo = rO2P - exp(log(rO2P)-error(pindex.lrO2P));
    name{nx}     = "rO2P"  ;
    R.xhat(nx)   = rO2P    ;
    R.upbar(nx)  = rO2P_up ;
    R.lowbar(nx) = rO2P_lo ;
end
%
% Si model parameters
if exist('xhat') & isfield(xhat,'dsi')
    nx = nx + 1 ;
    dsi = xhat.dsi ;
    dsi_up = exp(log(dsi)+error(pindex.ldsi)) - dsi;
    dsi_lo = dsi - exp(log(dsi)-error(pindex.ldsi));
    name{nx}     = "dsi"  ;
    R.xhat(nx)   = dsi    ;
    R.upbar(nx)  = dsi_up ;
    R.lowbar(nx) = dsi_lo ;
end

if exist('xhat') & isfield(xhat,'at')
    nx = nx + 1 ;
    at = xhat.at   ;
    at_up = exp(log(at)+error(pindex.lat)) - at;
    at_lo = at - exp(log(at)-error(pindex.lat));
    name{nx}     = "at"  ;
    R.xhat(nx)   = at    ;
    R.upbar(nx)  = at_up ;
    R.lowbar(nx) = at_lo ;
end

if exist('xhat') & isfield(xhat,'bt')
    nx = nx + 1 ;
    bt = xhat.bt   ;
    bt_up = exp(log(bt)+error(pindex.lbt)) - bt;
    bt_lo = bt - exp(log(bt)-error(pindex.lbt));
    name{nx}     = "bt"  ;
    R.xhat(nx)   = bt    ;
    R.upbar(nx)  = bt_up ;
    R.lowbar(nx) = bt_lo ;
end

if exist('xhat') & isfield(xhat,'aa')
    nx = nx + 1 ;
    aa = xhat.aa   ;
    aa_up = exp(log(aa)+error(pindex.laa)) - aa;
    aa_lo = aa - exp(log(aa)-error(pindex.laa));
    name{nx}     = "aa"  ;
    R.xhat(nx)   = aa    ;
    R.upbar(nx)  = aa_up ;
    R.lowbar(nx) = aa_lo ;
end

if exist('xhat') & isfield(xhat,'bb')
    nx = nx + 1 ;
    bb = xhat.bb   ;
    bb_up = exp(log(bb)+error(pindex.lbb)) - bb;
    bb_lo = bb - exp(log(bb)-error(pindex.lbb));
    name{nx}     = "bb"  ;
    R.xhat(nx)   = bb    ;
    R.upbar(nx)  = bb_up ;
    R.lowbar(nx) = bb_lo ;
end
% -------Cell Parameters-------------------
if exist('xhat') & isfield(xhat,'Q10Photo')
    nx = nx + 1 ;
    Q10Photo = xhat.Q10Photo  ;
    Q10Photo_up = exp(log(Q10Photo)+error(pindex.lQ10Photo)) - Q10Photo;
    Q10Photo_lo = Q10Photo - exp(log(Q10Photo)-error(pindex.lQ10Photo));
    name{nx}     = 'Q10Photo'  ;
    R.xhat(nx)   = Q10Photo    ;
    R.upbar(nx)  = Q10Photo_up ;
    R.lowbar(nx) = Q10Photo_lo ;
end

x_hat   = R.xhat       ;
upbar  = R.upbar      ;
lowbar = R.lowbar     ;
name   = string(name) ;
T = table(x_hat,upbar,lowbar,'RowNames',name)
