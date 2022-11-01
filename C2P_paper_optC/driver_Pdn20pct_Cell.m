clc; clear all; close all
global iter
%global fmin_history

fprintf([' NOTE: Phosphate decrease by 20 percent; \n'...
	'in Pdn20pct, only decrease po4obs; no decline in modelled DIP; \n'...
	'production does not change. \n'])

cd ..

% set up parallel pool
if isempty(gcp('nocreate'))
	poolobj = parpool(8); % if running in parallel; comment out to run serial
	%to close the parallel pool: delete(poolobj)
	poolobj.IdleTimeout = 120;
end

iter = 0 ;
on   = true  ;
off  = false ;
format long
%

% filename to load parameter values from (output of steady state optimization)(different from future projection run name)
%Cell model
%par.fxhatload = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper/optCell_CTL_He_PCCell01b_DOC0.25_DOP0_xhat.mat';
%par.fmodelload = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper/optCell_CTL_He_PCCell01b_DOC0.25_DOP0.mat';

%GM15
par.fxhatload = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0_xhat.mat';
par.fmodelload = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0.mat';

fprintf('load parameters from run: %s \n',par.fxhatload)

VerName = 'Pdn20pct_Cell_'; 		% optional version name. leave as an empty character array
					% or add a name ending with an underscore
VerNum = '';		% optional version number

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

par.optim   = off ;
par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.dynamicP = off ; % if on, cell model uses modeled DIP. if off, cell model uses WOA observed DIP field.

par.LoadOpt = on ; % if load optimial par.
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

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
%Trait Model parameters
par.opt_Q10Photo     = off ;
par.opt_fStorage     = off ; %
par.opt_fRibE 	     = off ;
par.opt_kST0 	     = off ; %
par.opt_PLip_PCutoff = off ; %
par.opt_PLip_scale   = off ;
par.opt_PStor_rCutoff = off; %
par.opt_PStor_scale  = off ;
par.opt_alphaS       = off ;
par.opt_gammaDNA	 = off ; %



%-------------load data and set up parameters---------------------
SetUp ;

% ----overwrite temperature obs from SetUp -----
if GridVer == 90
	%load tempobs_90x180x24.mat
    load po4obs_90x180x24.mat       % WOA PO4 observation [units: umol/kg]
	%load no3obs_90x180x24.mat		% WOA NO3 obs [units: umol/kg]
elseif GridVer == 91
	%load tempobs_91x180x24.mat
    load po4obs_91x180x24.mat % WOA PO4 observation
	%load no3obs_91x180x24.mat % WOA NO3 obs
end
iprod = find(M3d(:,:,1:2));
po4obs(iprod) = po4obs(iprod).*0.8;

% Uniform  20 percent Phosphate decrease in euphotic zone
par.po4obs  = po4obs  ;

%load model solution for DIP



% save results
% ATTENTION: Change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer);
elseif isunix
	output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/futureproj/');

end
VER = strcat(output_dir,VerName,TRdivVer);
%catDOC = sprintf('_DOC%0.2g_DOP%0.2g',par.cscale,par.pscale); % used to add scale factors to file names
catDOC = '';
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
    par = SetPar(par) ;
    % pack parameters into an array, assign them corresponding indices.
    par = PackPar(par) ;
% end

%------------------- run model -------------------------
%
x0    = par.p0 ;
nip = length(x0);

	if (par.optim)

	else
		iter = 0;
		xhat = PrintPara(x0, par) ;
		% run P cycle
		% tic
		% [par, P] = eqPcycle(x0, par) ;
	    % DIP = M3d+nan  ;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
	    % POP = M3d+nan  ;  POP(iwet) = P(1+1*nwet:2*nwet) ;
	    % DOP = M3d+nan  ;  DOP(iwet) = P(1+2*nwet:3*nwet) ;
		%
	    % par.DIP  = DIP(iwet) ;
	    % model.DIP = DIP ; model.POP = POP ; model.DOP = DOP ;
		% fprintf('P model solved \n')
		% toc

		if (par.opt_beta == on)
	        lbeta = x(pindx.lbeta);
	        beta  = exp(lbeta);
	    else
	        beta  = par.beta;
	    end
	    par.beta  = beta; % pass parameter to C/Si/O models
		% build part of the biological DIP uptake operator
	    Lambda     = par.Lambda;
	    LAM        = 0*M3d;
	    LAM(:,:,1) = (par.npp1.^beta).*Lambda(:,:,1);
	    LAM(:,:,2) = (par.npp2.^beta).*Lambda(:,:,2);
	    L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
	    par.L      = L;

		% load P model solution from optimal run and then decrease DIP by 20 percent in surface
		d = load(par.fmodelload);
		DIP = d.data.DIP;
		DOP = d.data.DOP;
		POP = d.data.POP;
		clear d;

		iprod = find(M3d(:,:,1:2));

		%DIP(iprod) = DIP(iprod).*0.80;
		par.DIP  = DIP(iwet) ;
	    model.DIP = DIP ; model.POP = POP ; model.DOP = DOP ;


		if (par.Cellmodel == on)
			tic
			[par, C2P] = eqC2Puptake(x0, par, model); % this line replaces the rest of this section

			%model.CellOut = par.CellOut;
			model.CellOut.C2P =  par.CellOut.C2P;
			model.CellOut.N2P = par.CellOut.N2P;
			model.CellOut.C2N = par.CellOut.C2N;
			model.CellOut.LimType = par.CellOut.LimType;
			model.CellOut.r = par.CellOut.r;
			model.CellOut.mu = par.CellOut.mu;
			model.CellOut.E = par.CellOut.E;
			model.CellOut.L = par.CellOut.L;
			model.CellOut.A = par.CellOut.A; % M=A
			model.CellOut.PLip = par.CellOut.PLip;
			model.CellOut.PStor = par.CellOut.PStor;
			model.CellOut.QP = par.CellOut.QP;
			model.CellOut.QC = par.CellOut.QC;
			model.CellOut.dC2P_dDIP = par.CellOut.dC2P_dDIP;

			fprintf('All of Cell model solved \n ')
			toc
		end

		if (par.Cmodel == on)
			tic
	        [par, C] = eqCcycle(x0, par) ;
	        DIC = M3d+nan ;  DIC(iwet) = C(0*nwet+1:1*nwet) ;
	        POC = M3d+nan ;  POC(iwet) = C(1*nwet+1:2*nwet) ;
	        DOC = M3d+nan ;  DOC(iwet) = C(2*nwet+1:3*nwet) ;
	        PIC = M3d+nan ;  PIC(iwet) = C(3*nwet+1:4*nwet) ;
	        ALK = M3d+nan ;  ALK(iwet) = C(4*nwet+1:5*nwet) ;

	        par.DIC = DIC(iwet) ;
	        par.DOC = DOC(iwet) ;
	        DIC = DIC + par.dicant  ;
	        model.DIC = DIC ;  model.POC = POC ;
	        model.DOC = DOC ;  model.PIC = PIC ;
	        model.ALK = ALK ;
			fprintf('C model solved \n')
			toc
	    end

		% --------- Save all parameter values inlcuding fixed values in xhat ------
		% this is to make sure all parameter values are consistent when rerunning model
		allparams.kappa_p 	= par.kappa_p;
		allparams.kappa_g 	= par.kappa_g   ;
		allparams.sigma 	= par.sigma     ;
		allparams.kP_T		= par.kP_T      ;
		allparams.kdP 		= par.kdP       ;
		allparams.bP_T 		= par.bP_T      ;
		allparams.bP 		= par.bP        ;
		allparams.alpha 	= par.alpha    ;
		allparams.beta 		= par.beta      ;
		if (par.Cmodel == on)
			allparams.kPIC = par.kPIC  ;
			allparams.bC_T = par.bC_T  ;
			allparams.bC = par.bC    ;
			allparams.d = par.d   ;
			allparams.kC_T = par.kC_T  ;
			allparams.kdC = par.kdC   ;
			allparams.R_Si = par.R_Si  ;
			allparams.rR = par.rR    ;
			allparams.cc = par.cc    ;
			allparams.dd = par.dd    ;
		end
		if (par.Omodel == on)
			allparams.O2C_T = par.O2C_T ;
			allparams.rO2C = par.rO2C  ;
			allparams.O2P_T = par.O2P_T ;
			allparams.rO2P = par.rO2P  ;
		end
		if (par.Simodel==on)
			allparams.dsi =  par.dsi   ;
			allparams.at = par.at    ;
			allparams.bt = par.bt    ;
			allparams.iaa = par.iaa   ;
			allparams.bb = par.bb  ;
		end
		if (par.Cellmodel==on)
			%pull params from BIO substructure
			%allparams.BIO = par.BIO;
			BIOfnames = fieldnames(par.BIO);
			for ii = 1:length(BIOfnames)
				allparams.(BIOfnames{ii}) = par.BIO.(BIOfnames{ii});
			end
			clear BIOfnames
		end
		model.allparams = allparams;

		% ----- SAVE SOLUTION -------
		fprintf('saving model solution to file: %s \n',par.fname)
	    save(par.fname, 'model')
		%save all parameter values used in a substructure in model

	end

delete(gcp('nocreate')) % shut down parallel pool
fprintf('-------------- end! ---------------\n');
