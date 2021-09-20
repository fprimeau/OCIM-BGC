%not currently used: does not produce any figures
%testing how to set up start of figure files
clc; clear all; close all
on   = true  ;
off  = false ;

%ver = datestr(now,'mmmdd');
RunVer = 'Tv4_PCCellv3b_DOC0.25_DOP0';

%model output directory
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/';
figDir = strcat(outputDir,'FIGS_PCCellv3b_DOC0.25_DOP0/');
outPath = figDir;

% load model output fields
fname = strcat(outputDir, RunVer, '.mat');
load(fname);
model = data;

% load optimal parameter values
fxhat = strcat(outputDir, RunVer,'_xhat.mat');
load(fxhat);

GridVer  = 90  ;
operator = 'A' ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

%-------------load data and set up parameters---------------------
SetUp ;

%------- read in which parameters were optimized ----------
%{
% default all parameters = off; then parameters in xhat are turned on
% P model parameters
par.opt_sigma = off ;
par.opt_kP_T  = off ;
par.opt_kdP   = off ;
par.opt_bP_T  = off ;
par.opt_bP    = off ;
par.opt_alpha = off ;
par.opt_beta  = off ;
% C model parameters
par.opt_bC_T  = off ;
par.opt_bC    = off ;
par.opt_d     = off ;
par.opt_kC_T  = off ;
par.opt_kdC   = off ;
par.opt_R_Si  = off ;
par.opt_rR    = off ;
par.opt_cc    = off ; %always off if cell model on
par.opt_dd    = off ; %always off if cell model on
% O model parameters
par.opt_O2C_T = off ;
par.opt_rO2C  = off ;
par.opt_O2P_T = off ;
par.opt_rO2P  = off ;
% Si model parameters
par.opt_dsi   = off  ;
par.opt_at    = off ;
par.opt_bt    = off  ;
par.opt_aa    = off  ;
par.opt_bb    = off  ;
%Trait Model parameters
par.opt_Q10Photo     = off ;
par.opt_fStorage     = off ;
par.opt_PLip_PCutoff = off ;
par.opt_PLip_scale   = off ;
par.opt_PStor_rCutoff = off;
par.opt_PStor_scale  = off ;
par.opt_alphaS       = off ;
par.opt_fRibE 	     = off ;
par.opt_kST0 	     = off;
% par.BIO.opt_gammaDNA = off;
% par.BIO.opt_gammaLipid = off;
% par.BIO.opt_DNT0 = off;
% par.BIO.opt_DPT0 = off;
% par.BIO.opt_Q10Diffusivity = off;
% par.BIO.opt_AMin = off;
% par.BIO.opt_PhiS = off;

%turn on parameters that were optimized in the run
parnames = fieldnames(xhat);

for ii=2:length(parnames)
    theField = ['opt_' parnames{ii}];
	par.(theField) = on;
end
%}


%--------------------- prepare parameters ------------------
xhat
% if par.LoadOpt == on
%     % load optimal parameters from a file or set them to default values
%     par = SetPar(par) ;
%     % pack parameters into an array, assign them corresponding indices.
%     par = PackPar(par) ;
% 	x0    = par.p0 ;
% 	iter = 0;
% 	PrintPara(x0, par) ;
% end

%{
%----------------surface indices---------------------------
iprod = find(par.M3d(:,:,1:2));
iprod1 = find(par.M3d(:,:,1)); % surface layer only
iprod2 = find(par.M3d(:,:,2)); % second EZ depth layer
%}

fprintf('-------------- end! ---------------\n');
