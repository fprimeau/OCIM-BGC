%% notes: imaginary values in the covariance matrix indicate:
% need to run cell model once more after the second deriv so the imaginary part isn't carried through?

clc; clear all; close all
on   = true  ;
off  = false ;

%RunVer = 'a0q1f2k1r0g1s2_CTL_He_PCCell_DOC0.25_DOP0'
%RunVer = 'optCell_CTL_He_PCCell01b_DOC0.25_DOP0'
%RunVer = 'optC_GM15_CTL_He_PC_DOC0.25_DOP0'
RunVer = 'optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0'

%model output directory
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/';
%figDir = strcat(outputDir,'FIGS_optCell01/v01b_');
figDir = strcat(outputDir,'FIGS_optC_Cell/');
%outPath = figDir;

% load model output fields
fname = strcat(outputDir, RunVer, '.mat');
load(fname);
model = data;
clear data

% load optimal parameter values
fxhat = strcat(outputDir, RunVer,'_xhat.mat');
load(fxhat);

GridVer  = 91  ;
operator = 'A' ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function
par.dynamicP = off;

%-------------load data and set up parameters---------------------

SetUp ;
xhat;
iwet = par.iwet;
iprod = find(par.M3d(:,:,1:2));

% need to turn on the opt_parmetername fields to use PackPar
% but we will need pindx later in this code. saved pindx in xhat instead

%----------- display all parameters -----------------------
fprintf('All Parameter Values: \n')
params = xhat.allparams

% fnames= fieldnames(xhat.allparams);
% for ii=1:length(fnames)
%     paramname=fnames{ii};
%     fprintf('%12s : %.5e \n', paramname, xhat.allparams.(paramname))
% end

fprintf('\n')
fprintf('Optimized Parameter Values: \n')
pnames = fieldnames(xhat);
for ii=1:length(xhat.fx)
    fprintf('%16s : %.5e \n', pnames{ii}, xhat.(pnames{ii}));
end

%% ------- Calculate Error Bars--------------
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

pindex = xhat.pindx ;
error  = sqrt(diag(inv(HH))) ;
nx = 0 ;
n  = numel(fieldnames(xhat.pindx)) ;

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
	alpha_opt = alpha;
	clear alpha
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

if exist('xhat') & isfield(xhat,'fStorage')
    nx = nx + 1 ;
    fStorage = xhat.fStorage  ;
    fStor_up = exp(log(fStorage)+error(pindex.lfStorage)) - fStorage;
    fStor_lo = fStorage - exp(log(fStorage)-error(pindex.lfStorage));
    name{nx}     = 'fStorage'  ;
    R.xhat(nx)   = fStorage    ;
    R.upbar(nx)  = fStor_up ;
    R.lowbar(nx) = fStor_lo ;
end

% up and lo units in error same as transformed units
if exist('xhat') & isfield(xhat,'fRibE')
    nx = nx + 1 ;
	%fRibE = 0.5*(1+tanh(tfRibE))
	%tfRibE = atanh(2*fRibE-1);
    fRibE = xhat.fRibE  ;
    fRibE_up = 0.5*(1+tanh( atanh(2*fRibE-1)+error(pindex.tfRibE))) - fRibE;
    fRibE_lo = fRibE - 0.5*(1+tanh( atanh(2*fRibE-1) - error(pindex.tfRibE)));
    name{nx}     = 'fRibE'  ;
    R.xhat(nx)   = fRibE    ;
    R.upbar(nx)  = fRibE_up ;
    R.lowbar(nx) = fRibE_lo ;
end

if exist('xhat') & isfield(xhat,'kST0')
    nx = nx + 1 ;
    kST0 = xhat.kST0   ;
    kST0_up = exp(log(kST0)+error(pindex.lkST0)) - kST0 ;
    kST0_lo = kST0 - exp(log(kST0)-error(pindex.lkST0));
    name{nx}     = 'kST0'  ;
    R.xhat(nx)   = kST0     ;
    R.upbar(nx)  = kST0_up ;
    R.lowbar(nx) = kST0_lo ;
end

if exist('xhat') & isfield(xhat,'PLip_PCutoff')
    nx = nx + 1 ;
    PCutoff = xhat.PLip_PCutoff   ;
    PCutoff_up = exp(log(PCutoff)+error(pindex.lPLip_PCutoff)) - PCutoff ;
    PCutoff_lo = PCutoff - exp(log(PCutoff)-error(pindex.lPLip_PCutoff));
    name{nx}     = 'PLip_PCutoff'  ;
    R.xhat(nx)   = PCutoff    ;
    R.upbar(nx)  = PCutoff_up ;
    R.lowbar(nx) = PCutoff_lo ;
end

if exist('xhat') & isfield(xhat,'PLip_scale')
    nx = nx + 1 ;
    PLipscale = xhat.PLip_scale  ;
    PLipscale_up = exp(log(PLipscale)+error(pindex.lPLip_scale)) - PLipscale ;
    PLipscale_lo = PLipscale - exp(log(PLipscale)-error(pindex.lPLip_scale));
    name{nx}     = 'PLip_scale'  ;
    R.xhat(nx)   = PLipscale    ;
    R.upbar(nx)  = PLipscale_up ;
    R.lowbar(nx) = PLipscale_lo ;
end

if exist('xhat') & isfield(xhat,'PStor_rCutoff')
    nx = nx + 1 ;
    rCutoff = xhat.PStor_rCutoff   ;
    rCutoff_up = exp(log(rCutoff)+error(pindex.lPStor_rCutoff)) - rCutoff ;
    rCutoff_lo = rCutoff - exp(log(rCutoff)-error(pindex.lPStor_rCutoff));
    name{nx}     = 'PStor_rCutoff'  ;
    R.xhat(nx)   = rCutoff    ;
    R.upbar(nx)  = rCutoff_up ;
    R.lowbar(nx) = rCutoff_lo ;
end

if exist('xhat') & isfield(xhat,'PStor_scale')
    nx = nx + 1 ;
    PStorscale = xhat.PStor_scale  ;
    PStorscale_up = exp(log(PStorscale)+error(pindex.lPStor_scale)) - PStorscale ;
    PStorscale_lo = PStorscale - exp(log(PStorscale)-error(pindex.lPStor_scale));
    name{nx}     = 'PStor_scale'  ;
    R.xhat(nx)   = PStorscale    ;
    R.upbar(nx)  = PStorscale_up ;
    R.lowbar(nx) = PStorscale_lo ;
end

if exist('xhat') & isfield(xhat,'alphaS')
    nx = nx + 1 ;
    alphaS = xhat.alphaS   ;
    alphaS_up = exp(log(alphaS)+error(pindex.lalphaS)) - alphaS ;
    alphaS_lo = alphaS - exp(log(alphaS)-error(pindex.lalphaS));
    name{nx}     = 'alphaS'  ;
    R.xhat(nx)   = alphaS     ;
    R.upbar(nx)  = alphaS_up ;
    R.lowbar(nx) = alphaS_lo ;
end

if exist('xhat') & isfield(xhat,'gammaDNA')
    nx = nx + 1 ;
    gammaDNA = xhat.gammaDNA   ;
	gammaDNA_up = 0.5*(1+tanh( atanh(2*gammaDNA-1)+error(pindex.tgammaDNA))) - gammaDNA;
    gammaDNA_lo = gammaDNA - 0.5*(1+tanh( atanh(2*gammaDNA-1) - error(pindex.tgammaDNA)));
    name{nx}     = 'gammaDNA'  ;
    R.xhat(nx)   = gammaDNA     ;
    R.upbar(nx)  = gammaDNA_up ;
    R.lowbar(nx) = gammaDNA_lo ;
end

%-----------Make Table-----------------------
x_hat   = R.xhat       ;
upbar  = R.upbar      ;
lowbar = R.lowbar     ;
name   = string(name) ;
T = table(x_hat,upbar,lowbar,'RowNames',name)

%------- Cell Output Stats -----------

%%------------- Calculate C export ------------------
% ----------------------------------------------
nn = par.nzo ; %number fo verticle boxes in euphotic zone / export depth

% -------------- C:P uptake ratio --------------------------------
W = d0(dVt(iwet)) ;
dVtwet = M3d*nan;
dVtwet(iwet) = dVt(iwet);
Wiprod = dVtwet(:,:,1:nn)/nansum(dVtwet(:,:,1:nn),'all');

C2P3D = M3d + nan ;
if par.Cellmodel==on
	C2P3D(iwet) = model.CellOut.C2P(iwet); %zero beneath the surface  layers
elseif par.Cmodel ==on
	C2P3D(iwet) = 1./(params.cc*par.po4obs(iwet) + params.dd) ;  % DIP or PO4?
	%C2P3D(iwet) = 1./(7e-4*PO4 + 5e-3); % WL
	%C2P3D(iwet) = 1000./(6.6*PO4 + 5.3); %Qian
else
	C2P3D = M3d +nan;
	C2P3D(iwet) = 106; 		% redfield C:P
end

% DIP assimilation
LAM        = 0*M3d;
LAM(:,:,1) = (par.npp1.^params.beta).*par.Lambda(:,:,1);
LAM(:,:,2) = (par.npp2.^params.beta).*par.Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];

%--------------- calculate primary production --------------------
fprintf('----------- Cexp ------------- \n')
G        = M3d*0        ;
G(iwet)  = params.alpha*L*model.DIP(iwet)  ; % primary production [unit: mmol P/m^3/s]

% inegG = find(G<0); % should negative production values be removed?
% G(inegG)=nan;

Int_PNPP = G(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31;
Int_CNPP = G(:,:,1:nn).*grd.DZT3d(:,:,1:nn).*C2P3D(:,:,1:nn)*12;

PNPP = Int_PNPP*spa*1e-3 ; % convert to g P/m^2/yr
CNPP = Int_CNPP*spa*1e-3 ; % convert production from mg C/m^3/s to gC/m^2/year;
tem_PNPP = PNPP.*dAt(:,:,1:nn)*1e-15 ;
tem_CNPP = CNPP.*dAt(:,:,1:nn)*1e-15 ;
Sum_CNPP = nansum(tem_CNPP(:))    ;
fprintf('Model NPP (P) is %3.3e Pg P/yr \n',nansum(tem_PNPP(:))) ; %Pg/yr
fprintf('Model NPP is %3.3e Pg C/yr \n\n',Sum_CNPP) ; %Pg/yr

clear tem_PNPP tem_CNPP

%----- C2P --------------
C2Pavg = nansum(C2P3D(:,:,1:nn).*Wiprod,'all');
fprintf('Average C:P uptake in Euphotic Layers is %4.2f \n',C2Pavg);

% prod weighted C:P
WeightNPP = CNPP.*dAt(:,:,1:nn)/nansum(CNPP.*dAt(:,:,1:nn),'all'); %gC/yr
C2Pavg = nansum(C2P3D(:,:,1:nn).*WeightNPP,'all');
fprintf('Average C:P uptake in Euphotic Layers (NPP weighted) is %4.2f \n\n',C2Pavg);


% ------- P export -----
%POP export: integrte POP beneath the euphotic zone
POPexp = nansum(params.kappa_p*model.POP(:,:,3:end).*dVt(:,:,3:end),3)*31*spd;
Sum_POPexp =nansum(POPexp(:))*365*1e-18;
fprintf('Model POP export is %3.3e Pg P /yr  (beneath EZ) \n',Sum_POPexp);

%integrated DOP remineralization below the euphotic zone. (should equal the TOP export calculated by the adjoint method)
DOPexpint = params.kdP*model.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
tem_DOPexpint = nansum(DOPexpint.*dAt(:,:,3:end),3);
Sum_DOPexpint =nansum(tem_DOPexpint(:))*365*1e-18;
fprintf('Model TOP export (integrated DOP below the Euphotic zone) is %3.3f Pg P /yr \n\n',Sum_DOPexpint);
DOPexpint = sum(DOPexpint,3,'omitnan');


%----- C export ---------------
%POC export: integrate POC beneath the euphotic zone ----------
%POCexp = nansum(par.kappa_p*model.POC(:,:,3:end).*dVt(:,:,3:end),3)*12*spd;  % [1/s]*[mmol/m^3]*[m^3]*[mg/mmol]*[s/day] =  [mg/day]
POCexp = params.kappa_p*model.POC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
tem_POCexp = nansum(POCexp.*dAt(:,:,3:end),3);
Sum_POCexp =nansum(tem_POCexp(:))*365*1e-18;
fprintf('Model POC export is %3.3e Pg C /yr  (beneath EZ) \n',Sum_POCexp);
POCexp = sum(POCexp,3,'omitnan');

%integrated DOC remineralization below the euphotic zone. (should equal the TOC export calculated by the adjoint method)
DOCexpint = params.kdC*model.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
tem_DOCexpint = nansum(DOCexpint.*dAt(:,:,3:end),3);
Sum_DOCexpint =nansum(tem_DOCexpint(:))*365*1e-18;
fprintf('Model TOC export (integrated DOC below the Euphotic zone) is %7.5f Pg C /yr \n\n',Sum_DOCexpint);
DOCexpint = sum(DOCexpint,3,'omitnan');


%%---------- compare to Obs ------------------
fprintf('------- Compare to Obs ----------- \n')
idop   = find(par.dopraw(iwet) > 0 & model.DOP(iwet)>0) ;
fprintf('R^2 for DOP is %3.4f \n',rsquare(par.dopraw(iwet(idop)),model.DOP(iwet(idop))))

ipo4 = find(model.DIP(iwet) > 0 & par.po4raw(iwet) > 0.02);
O = par.po4raw(iwet(ipo4));
M = model.DIP(iwet(ipo4));
fprintf('R^2 for DIP is %3.4f \n',rsquare(O,M))

% surface only DIP
%ipo4 = find(DIP(iwetsurf) > 0 & po4raw(iwetsurf) > 0.02);
ipo4 = find(model.DIP(iprod) & ~isnan(par.po4raw(iprod)) );
O = par.po4raw(iprod(ipo4));
M = model.DIP(iprod(ipo4));
fprintf('R^2 for surface DIP is %3.4f \n',rsquare(O,M))

DICmod = model.DIC - par.dicant;
DICobs = par.dicraw ;
iDIC = find(DICobs(iwet)>0);
%
O = DICobs(iwet(iDIC));
% already including anthropogenic CO2
% M = DIC(iwet(iDIC)) ;
% not include anthropogenic CO2
M = DICmod(iwet(iDIC))+par.dicant(iwet(iDIC));
fprintf('R^2 for DIC is %3.4f \n',rsquare(O,M))

% ---- surface only DIC -----
iDIC = find(DICobs(iprod)>0);
O = DICobs(iprod(iDIC));
% already including anthropogenic CO2
% M = DIC(iwet(iDIC)) ;
% not include anthropogenic CO2
M = DICmod(iprod(iDIC))+par.dicant(iprod(iDIC));
fprintf('R^2 for surface DIC is %3.4f \n',rsquare(O,M))

% -- ALK ----
iALK = find(par.alkraw(iwet)>0);
%
O = par.alkraw(iwet(iALK));
M = model.ALK(iwet(iALK)); % already including anthropogenic CO2
fprintf('R^2 for ALK is %3.4f \n',rsquare(O,M))

%---- DOC ---------
iDOC = find(DOCclean(iwet)>0 & model.DOC(iwet)>0) ;
O = DOCclean(iwet(iDOC)) ;
M = model.DOC(iwet(iDOC)) ;
fprintf('R^2 for DOC is %3.4f \n',rsquare(O,M))
