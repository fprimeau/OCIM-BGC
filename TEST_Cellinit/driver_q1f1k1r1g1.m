clc; clear all; close all
global iter
%global fmin_history

fprintf(['\n' 'NOTE: test cell model initial values \n' ...
	'run model with Cell model C:P \n' ...
	'BGC parameters fixed to output of GM15 model run \n'...
	'PLip_PCutoff fixed at 1e-6 \n'...
	'sigma fixed at zero, \n'...
	'beta fixed at 0.5 \n'...
	'alpha = 3.78e-5'
	'Q10Photo on initial = 1 \n'...
	'fStorage initialized at 1e-3 \n'...
	'dynamicP = off \n'])

cd ..

% set up parallel pool for cell model
if isempty(gcp('nocreate'))
	poolobj = parpool(2); % if running in parallel; comment out to run serial
	%to close the parallel pool: delete(poolobj)
	poolobj.IdleTimeout = 120;
end

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

VerName = 'a0q1f1k1r1g1_'; 		% optional version name. leave as an empty character array
					% or add a name ending with an underscore
VerNum = '';		% optional version number

Gtest = off ;		% do gradient test
Htest = off ;		% do hessian test
par.optim   = on ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.dynamicP = off ; % if on, cell model uses modeled DIP. if off, cell model uses WOA observed DIP field.

par.LoadOpt = on ; % if load optimial par.
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

% to load parameter values from a run with a different name. need to delete or comment out for loadOpt to work normally
par.fxhatload = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK91/testNPP_CTL_He_PCfixb_DOC0.25_DOP0_xhat.mat' ;
%par.fxhatload = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/ESS225_xhat.mat' ;

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
par.opt_bPC	  = off; %new variable to optimize bP and bC with the same value.
					 % no temperature dependence. must have opt_bP on and opt_bP_T, opt_bC, opt_bC_T off
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
par.opt_Q10Photo     = on ;
par.opt_fStorage     = on ; %
par.opt_fRibE 	     = off ;
par.opt_kST0 	     = on ; %
par.opt_PLip_PCutoff = off ; %
par.opt_PLip_scale   = off ;
par.opt_PStor_rCutoff = on; %
par.opt_PStor_scale  = off ;
par.opt_alphaS       = off ;
par.opt_gammaDNA	 = on ; %
% par.BIO.opt_gammaLipid = off;
% par.BIO.opt_DNT0 = off;
% par.BIO.opt_DPT0 = off;
% par.BIO.opt_Q10Diffusivity = off;
% par.BIO.opt_AMin = off;
% par.BIO.opt_PhiS = off;

if par.opt_bPC && (par.opt_bC || par.opt_bC_T)
	fprintf('Error: cannot have opt_bPC turned on if opt_bC or opt_bC_T are on');
	return
end

%-------------load data and set up parameters---------------------
SetUp ;

% save results
% ATTENTION: Change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer);
elseif isunix
    % output_dir = sprintf('/DFS-L/DATA/primeau/weilewang/Cexp/');
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
                        % 'MSK%2d/PME4DICALK/'],GridVer);

	%output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/test%s/', datestr(now,'ddmmmyy'));
	output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/testCellinit/', GridVer);
	%output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper/');
end
VER = strcat(output_dir,VerName,TRdivVer);
catDOC = sprintf('_DOC%0.2g_DOP%0.2g',par.cscale,par.pscale); % used to add scale factors to file names
% Creat output file names based on which model(s) is(are) optimized
%if Gtest == on
%    fname = strcat(VER,'_GHtest');
%elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off & par.Cellmodel == off)
        fname = strcat(VER,'_P',VerNum);
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off & par.Cellmodel == off)
        base_name = strcat(VER,'_PC',VerNum);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == off)
        base_name = strcat(VER,'_PCO',VerNum);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on & par.Cellmodel == off)
        base_name = strcat(VER,'_PCSi',VerNum);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == off)
        base_name = strcat(VER,'_PCOSi',VerNum);
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == off & par.Omodel == off & par.Simodel == off & par.Cellmodel == on) % cell model does nothing if C model is not on, so this case =Ponly
        base_name = strcat(VER,'_PCell',VerNum);
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off & par.Cellmodel == on)
        base_name = strcat(VER,'_PCCell',VerNum);
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == on)
		base_name = strcat(VER,'_PCOCell',VerNum);
		fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == on)
        base_name = strcat(VER,'_PCOSiCell',VerNum);
        fname = strcat(base_name,catDOC);
    end
%end
% -------------------- Set up output files ---------------
par.fname = strcat(fname,'.mat') ;
fxhat     = strcat(fname,'_xhat.mat') %;
par.fxhat = fxhat ;
if Htest ==on
	fGHtest = strcat(fname,'_GHtest.mat')  ;
end
%if sensitivityTest ==on
%	fSensiTest = strcat(fname,'_pert.mat')  ;
%end

%par.fhistory = strcat(fname,'_history.mat');

% -------------------update initial guesses --------------
if isfile(par.fname)
    load(par.fname)

	% make sure these fields are real values only
	fnames = fieldnames(data);
	for ii = 1:length(fnames)
		if isnumeric(data.(fnames{ii}))
			data.(fnames{ii}) = real(data.(fnames{ii}));
		end
	end
end

%---------------- inital guesses on C and O ---------------
if par.Cmodel == on
    DIC = data.DIC - par.dicant ;

    GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); ...
           data.PIC(iwet); data.ALK(iwet)];

	% temprarily use following line to create ALK initial field for 90x180 grid
	% GC = [DIC(iwet); data.POC(iwet); data.DOC(iwet); data.PIC(iwet); data.DIC(iwet)];

	% GC is used inside eqCcycle. Should GC be set to a global variable within this script (driver)?
	% global GC is called at  start of eqCcycle.
    GC  = GC + 1e-6*randn(5*nwet,1) ;
end
if par.Omodel == on
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end

%--------------------- prepare parameters ------------------
% if (par.optim == on) | (par.LoadOpt == on)
    % load optimal parameters from a file or set them to default values
	%cd TEST_Cellinit
    %par = SetPar_Cellv1(par) ;
	%cd ..

	spd  = 24*60^2 ;  	% seconds per day
	spa  = 365*spd ;  	% seconds per year
    % fixed parameters
    par.kappa_g  = 1/(1e6*spa)  ; % geological restoring time [1/s];
    par.taup     = 30*spd       ; % [s] pic dissolution time-scale
    % par.tau_TA   = 1./par.taup  ;
    par.kappa_p  = 1/(30*spd) ;   % [1/s] POM solubilization rate constant (fixed) (for both POP and POC)
								  % the parameters kappa_p and b, which affect
								  % the sinking flux attenuation profile, are not
								  % independently identified by the data. We therefore
								  % prescribe the value kappa_p of 1/(30 days)
    par.tauPIC = 30*spd ;   	% [s] The CaCO3 e-folding dissolution time is 30 days. converted to seconds.
    par.kPIC   = 1/par.tauPIC ;
								% PIC dissolution constant 0.38 day^-1 based on
								% first-order reaction kinetics according to
								% Sarmiento and Gruber book (p.271);

    % load optimal parameters if they exist
    % if isfile(par.fxhat) & par.LoadOpt == on
    %     load(par.fxhat)
    % end
	if par.LoadOpt ==on
		if isfield(par,'fxhatload') & isfile(par.fxhatload)
			load(par.fxhatload)
			% keep all parameters the same, including those that were not optimized
			xhat = xhat.allparams;
		elseif isfile(par.fxhat)
			load(par.fxhat)
		end
	end

    % P model parameters
    if exist('xhat') & isfield(xhat,'sigma')
        par.sigma = xhat.sigma ;
    else
        %par.sigma = 0.30 ;  			% Fraction of organic P production allocated directly to the dissolved pool
										% default was 1.0e-01 , however 0.30 is used by (Wang, 2019: Nitrogen fixation)
										% SONTrap value = 1/3
		%par.sigma = 0.10 ;
		par.sigma = 0.0 ; 				% all DOM comes from particles. no direct mechanism creates DOM during production.
    end
    if exist('xhat') & isfield(xhat,'kP_T')
        par.kP_T = xhat.kP_T ;
    else
        par.kP_T = 0.00 ;				% linear temperature dependence of DOP remineralization
    end
    if exist('xhat') & isfield(xhat,'kdP')
        par.kdP = xhat.kdP ;
    else
        %par.kdP = 2.42e-08 ;			% DOP remineralization constant [s^-1] (SONTrap = 1.57e-7)
		par.kdP = 1.23e-08 ;
    end
    if exist('xhat') & isfield(xhat,'bP_T')
        par.bP_T = xhat.bP_T ;
    else
        par.bP_T = 0.00e+00 ;			% Martin exponent of POP solubilization (slope)
    end
    if exist('xhat') & isfield(xhat,'bP')
        par.bP  = xhat.bP ;
    else
        %par.bP  = 9.60e-01 ;			% Martin exponent of POP solubilization (const term) (SONTrap paper =0.7)
		par.bP = 1.24;
    end
    if exist('xhat') & isfield(xhat,'alpha')
        par.alpha = xhat.alpha ;
    else
        %par.alpha = 3.37e-08   ;		% WeiLei's NPP scaling factor for DIP uptake rate (used until v9b)
		%par.alpha = 1.0   ;		% NPP scaling factor for DIP uptake rate
		%par.alpha = 1.0e-06 ;
		par.alpha = 3.78e-05;
    end
    if exist('xhat') & isfield(xhat,'beta')
        par.beta = xhat.beta ;
    else
        %par.beta = 1.65e-02  ;			% Weilei's NPP scaling exponent for DIP uptake rate (used until v9b)
		%par.beta = 1.0 ;			% NPP scaling exponent for DIP uptake rate
		%par.beta = 1.0e-03 ;
		%par.beta = 5.93e-01 ; 	% inbetween value to test. result of PCa1e-8b1e-3
		par.beta = 0.5 ;
    end
	% adding a new parameter for b, same for both P and C; incomplete
	if exist('xhat') & isfield(xhat,'bPC')
        par.bPC = xhat.bPC ;
    else
        par.bPC = 9.60e-1    ;			% Martin exponent of POC solubilization set equal to that of POP
    end

    % C model parameters
    if exist('xhat') & isfield(xhat,'bC_T')
        par.bC_T = xhat.bC_T ;
    else
        par.bC_T =  0.00e+00 ;			% Martin exponent of POP solubilization (linear T dependence)
    end
    if exist('xhat') & isfield(xhat,'bC')
        par.bC = xhat.bC ;
    else
        %par.bC = 9.38e-01    ;			% Martin exponent of POC solubilization (const term)
		%par.bC = 9.60e-01    ;			% set to match bP
		par.bC = 9.71e-01    ;
    end
    if exist('xhat') & isfield(xhat,'d')
        par.d = xhat.d   ;
    else
        %par.d = 4.54e+03;    			% CaCO3 dissolution length scale [m]
		par.d = 4.80e+03;
    end
    if exist('xhat') & isfield(xhat,'kC_T')
        par.kC_T = xhat.kC_T ;
    else
        par.kC_T = 0.00e+00 ;			%  DOC remineralization constant linear T dependence
    end
    if exist('xhat') & isfield(xhat,'kdC')
        par.kdC = xhat.kdC ;
    else
        %par.kdC = 7.35e-08 ;			% DOC remineralization constant (const. term) [s^-1]
		par.kdC = 3.25e-08 ;
    end
    if exist('xhat') & isfield(xhat,'R_Si')
        par.R_Si = xhat.R_Si ;
    else
        par.R_Si = 0.00 ;				% CaCO3 to POC production ratio linear silica dependence
    end
    if exist('xhat') & isfield(xhat,'rR')
        par.rR = xhat.rR  ;
    else
        %par.rR = 2.64e-02 ;				% CaCO3 to POC production ratio (const. term)
		par.rR = 3.22e-02 ;
    end
    if exist('xhat') & isfield(xhat,'cc')
        par.cc = xhat.cc  ;
    else
        %par.cc = 7.51e-4 ;				% slope for P:C as a linear function of DIP (WeiLei's value)
		%par.cc = 1.02e-7 ; 				% WeiLei's new value
		par.cc = 6.9e-3 ; 				% value from Galbraith & Martiny 2015
		%par.cc = 0.0 ;					% slope for P:C as a linear function of DIP. if cc is off, P:C is a constant
										% with cc initially set to 0.0, optimization doesn't work for this param. instead make it a very small value
    end
    if exist('xhat') & isfield(xhat,'dd')
        par.dd = xhat.dd  ;
    else
        %par.dd = 5.56e-03 ;			% intercept for P:C as a linear function of DIP (WeiLei's Value)
		par.dd = 6.0e-03 ;				% intercept for P:C as a linear function of DIP (value from Galbraith & MArtiny 2015)
		%par.dd = 1/106 ;				% intercept for P:C as a linear function of DIP (Redfield)
    end

    % O model parameters
    if exist('xhat') & isfield(xhat,'O2C_T')
        par.O2C_T = xhat.O2C_T ;
    else
        par.O2C_T = 0.00 ;
    end
    if exist('xhat') & isfield(xhat,'rO2C')
        par.rO2C = xhat.rO2C ;
    else
        par.rO2C = 1.10e+00 ;
    end
    if exist('xhat') & isfield(xhat,'O2P_T')
        par.O2P_T = xhat.O2P_T ;
    else
        par.O2P_T = 0.00 ;
    end
    if exist('xhat') & isfield(xhat,'rO2P')
        par.rO2P = xhat.rO2P ;
    else
        par.rO2P = 1.70e+02 ;
    end
    %
    % Si model parameters
    if exist('xhat') & isfield(xhat,'dsi')
        par.dsi = xhat.dsi ;
    else
        par.dsi = 3300     ;
    end
    if exist('xhat') & isfield(xhat,'at')
        par.at = xhat.at   ;
    else
        par.at = 1.32e16/spd;
    end
    if exist('xhat') & isfield(xhat,'bt')
        par.bt = xhat.bt   ;
    else
        par.bt = 11481     ;
    end
    if exist('xhat') & isfield(xhat,'aa')
        par.aa = xhat.aa   ;
    else
        par.aa = 1         ;
    end
    if exist('xhat') & isfield(xhat,'bb')
        par.bb = xhat.bb   ;
    else
        par.bb = 0.968     ;
    end
	%
	% Cell model parameters
	if exist('xhat') & isfield(xhat,'Q10Photo')
		par.BIO.Q10Photo = real(xhat.Q10Photo);
	else
		%v1
		par.BIO.Q10Photo = 1.0;
		%v2
		%par.BIO.Q10Photo = 2.0;		% Q10 of photosynthesis. v1 (default = 1.983)
		%par.BIO.Q10Photo = 1.88; 	% Eppley value (measured in natural communities)
		%par.BIO.Q10Photo = 1.46; 	% Anderson et al. 2021 (median for all phytoplankton; compilation of culture studies)
	end
	if exist('xhat') & isfield(xhat,'fStorage')
		par.BIO.fStorage = real(xhat.fStorage);
	else
		%v1
		par.BIO.fStorage = 0.001; %0.7 %1; %exp(-.358);  % strength of luxury P storage [L/molC]
		%v2
		%par.BIO.fStorage = 1; %0.7 %1; %exp(-.358);  % strength of luxury P storage [L/molC]
	end
	if exist('xhat') & isfield(xhat,'fRibE')
		par.BIO.fRibE = real(xhat.fRibE);
	else
		par.BIO.fRibE = .618;           % ribosome fraction of biosynthetic apparatus
	end
	if exist('xhat') & isfield(xhat,'kST0')
		par.BIO.kST0 = real(xhat.kST0);
	else
		%v1
		par.BIO.kST0 =0.185;            % specific synthesis rate of synthetic apparatus at 25degC [1/hr]
										% empirically derived specific efficiency: 0.168 hr^-1 (shuter 1979)
		%v2
		%par.BIO.kST0 =0.115;            % specific synthesis rate of synthetic apparatus at 25degC [1/hr]
	end
	if exist('xhat') & isfield(xhat,'PLip_PCutoff')
		par.BIO.PLip_PCutoff = real(xhat.PLip_PCutoff);
	else
		%par.BIO.PLip_PCutoff = exp(-14.408);  % log of [P (mol/L)] below which more PLipids are substituted with Slipids
											  % original default = exp(-14.408) = 5.5295e-07
		par.BIO.PLip_PCutoff = 1.0e-06;
	end
	if exist('xhat') & isfield(xhat,'PLip_scale')
		par.BIO.PLip_scale = real(xhat.PLip_scale);
	else
		par.BIO.PLip_scale = 3.0e6;  % scale factor for logistic function controlling phospholipid quota (changed from 1.0 b.c. not using log(P). changed from 1e6 to 3e6 to make transition sharper)
	end
	if exist('xhat') & isfield(xhat,'PStor_rCutoff')
		par.BIO.PStor_rCutoff = real(xhat.PStor_rCutoff);
	else
		%par.BIO.PStor_rCutoff = 2.25;  % radius [um] above which cell stores luxury phosphorus?
		%v1
		par.BIO.PStor_rCutoff = 1;
		%v2
		%par.BIO.PStor_rCutoff = 10;
	end
	if exist('xhat') & isfield(xhat,'PStor_scale')
		par.BIO.PStor_scale = real(xhat.PStor_scale);
	else
		par.BIO.PStor_scale = 3.00;  % scale factor for logistic function controlling luxury phosphorus storage (changed default from 1 to 3 to give sharper transition)
	end
	if exist('xhat') & isfield(xhat,'alphaS')
		par.BIO.alphaS = real(xhat.alphaS);
	else
		par.BIO.alphaS = 0.225;          % radius at which cell is all periplasm and membrane [um]
	end
	if exist('xhat') & isfield(xhat,'gammaDNA')
		par.BIO.gammaDNA = real(xhat.gammaDNA);
	else
		%par.BIO.gammaDNA = 0.016;          % DNA Fraction of cell
		%v1
		par.BIO.gammaDNA = 0.005;
		%v2
		%par.BIO.gammaDNA = 0.030;
	end

% cell model parameters that don't change
	if (par.Cellmodel==on)
		%if exist('xhat') & isfield(xhat,'gammaLipid')
		%par.BIO.gammaDNA 	= .016;		% DNA fraction of cell
		par.BIO.gammaLipid 	= .173;		% structural Lipid (non-membrane or periplasm) fraction of cell
		%par.BIO.lPCutoff 	= -7.252;	% log of max [P] for which Plipids will be substituted with Slipids
		%par.BIO.r0Cutoff 	= 2.25;		% % NEED TO REDEFINE: r0Cutoff =  rFullA; ASK GEORGE % now PStor_rCutoff
		par.BIO.DNT0 		= 1e-12*3.6e2*3600;    % Diffusivity of Nitrate at 25degC [m^2/hr]
		par.BIO.DPT0 		= 1e-12*3.6e2*3600;    % Diffusivity of Phosphate at 25degC [m^2/hr]
		par.BIO.Q10Diffusivity = 1.5;
		par.BIO.AMin 		= .05;		% minimal fraction of cell dry mass that is nutrient uptake proteins
		%par.BIO.CStor 		= 1.00;		% replaced by PStor_scale
		par.BIO.PhiS 		= .67;		% specific carbon cost of synthesis [gC/gC]
		%%% BIO parameters below should remain fixed
		par.BIO.pDry 		= .47;		% Dry mass fraction of the cell
		par.BIO.rho 		= 1e-12;	% cell density [g/um^3]
		par.BIO.fProtM 		= 0.25;		% protein fraction of cell membranes
		par.BIO.fProtL 		= .7;		% protein fraction of light harvesting apparatus
		par.BIO.PDNA 		= .095;		% phosphorus mass fraction in DNA [gP/g]
		par.BIO.PRib 		= 0.047;	% phosphorus mass fraction in ribosomes [gP/g]
		par.BIO.PPhospholipid = 0.042;	% phosphorus mass fraction in phospholipids [gP/g]
		par.BIO.NProt 		= .16;		% nitrogen mass fraction in proteins [gN/g]
		par.BIO.NDNA 		= .16;		% nitrogen mass fraction in DNA [gN/g]
		par.BIO.NRib 		= .16;		% nitrogen mass fraction in Ribosomes [gN/g]
		par.BIO.CProt 		= .53;		% carbon mass fraction in proteins [gC/g]
		par.BIO.CDNA 		= .36;		% carbon mass fraction in DNA [gC/g]
		par.BIO.CPhospholipid = .65;	% carbon mass fraction in phospholipids [gC/g] - why seperate form lipids?
		par.BIO.CLipid 		= .76;		% carbon mass fraction in other lipids (that are not phospholipids) [gC/g]
		par.BIO.CRib 		= .419;		% carbon mass fraction in ribosomes [gC/g] (technically, only correct for eukaryotes)
		par.BIO.alphaPLip 	= 0.12;		% phospholipid fraction of cell membrane
	end

	clear xhat


%--------------------- pack parameters ------------------
    % pack parameters into an array, assign them corresponding indices.
    par = PackPar(par) ;
% end

%-------------------set up fminunc -------------------------
x0    = par.p0 ;
myfun = @(x) neglogpost(x, par);

objfuntolerance = 5e-12;
options = optimoptions(@fminunc                  , ...
                       'Algorithm','trust-region', ...
                       'GradObj','on'            , ...
                       'Hessian','on'            , ...
                       'Display','iter'          , ...
                       'MaxFunEvals',2000        , ...
                       'MaxIter',2000            , ...
					   'StepTolerance',objfuntolerance      , ...
                       'OptimalityTolerance',objfuntolerance, ...
                       'DerivativeCheck','off'   , ...
                       'FinDiffType','central'   , ...
                       'PrecondBandWidth',Inf    , ...
					   'SubproblemAlgorithm','factorization') ; %testing this
					   %'OutputFcn',@fminoutfun ) 	 ; % maybe use 'SubproblemAlgorithm','factorization' -changing subproblem did not seem to effect anything
					   %'TolX',5e-7               , ...
                       %'TolFun',5e-7             , ...
%
nip = length(x0);
if (Gtest);
	fprintf('derivative test \n')
	GHtest.fx_cstep = NaN([nip,1]);
	GHtest.fxx_cstep = NaN(nip);
	GHtest.pindx = par.pindx;
    dx = sqrt(-1)*eps.^3*eye(nip);
    % for ii = nip : -1 : 13
    for ii = 1 : nip
        x  = real(x0)+dx(:,ii);
		iter = 5; %bypasses the BadStep check and extra save commands
        if Htest == on
            [f,fx,fxx] = neglogpost(x, par) ;
			GHtest.fx_cstep(ii) = imag(f)/eps.^3     ;
			GHtest.fxx_cstep(:,ii) = imag(fx)/eps.^3 ;
			%save complex step gradient and hessian
			fprintf('saving GHtest to file: %s \n',fGHtest)
			save(fGHtest,'GHtest')

            % print relative errors
            diff = (real(fx(ii)) - imag(f)/eps.^3)/(imag(f)/eps.^3);
            fprintf('gradient relative error (%i): % .3e  \n',ii,diff);
            diffx = (real(fxx(:,ii))-imag(fx)/eps.^3)./(imag(fx)/eps.^3+eps.^3);
            for jj = 1:length(fx)
                fprintf('% .3e  ', diffx(jj));
            end
        else
            [f,fx] = neglogpost(x, par) ;
            diff = (real(fx(ii)) - imag(f)/eps.^3)/(imag(f)/eps.^3) ;
            fprintf('%i % .3e  \n',ii,diff);
            fprintf('\n');
        end
        fprintf('\n');
    end
	format shortE
	real(fx)
	real(full(fxx))
    %keyboard
else
	if (par.optim)
	    [xsol,fval,exitflag] = fminunc(myfun,x0,options);
		fprintf('objective function tolerance = %5.1e \n',objfuntolerance);
		fprintf('----fminunc complete----\n')
	    [f,fx,fxx,data,xhat] = neglogpost(xsol,par);
		fprintf('----neglogpost solved for final parameter values----\n')
	    %load(fxhat,'xhat')
		xhat.pindx = par.pindx;
	    xhat.f   = f   ;
	    xhat.fx  = fx  ;
	    xhat.fxx = fxx ;
	    % save results
		fprintf('saving optimized parameters to file: %s \n',fxhat)
		fprintf('saving model solution to file: %s \n',par.fname)
	    save(fxhat, 'xhat')
	    save(par.fname, 'data')
		%save(par.fhistory,'fmin_history')
	else
		xsol = x0;
		iter = 6;
		[f,fx,fxx,data] = neglogpost(xsol,par);
		fprintf('----neglogpost complete----\n')
		%fprintf('saving model solution to file: %s \n',par.fname)
	    %save(par.fname, 'data')
	end
end
if ~isempty(gcp('nocreate'))
	delete(gcp('nocreate')) % shut down parallel pool
end
fprintf('-------------- end! ---------------\n');
