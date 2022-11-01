clc; clear all; close all
global iter
%global fmin_history

fprintf(['\n' 'NOTE: optimize C cycle with Cell model for C:P \n' ...
	'hold P cycle parameters fixed to values in optPonly \n'...
	'sigma fixed at zero, \n'...
	'dynamicP = off \n'])

cd ..

% set up parallel pool for cell model
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

VerName = 'optC_Cellv2_'; 		% optional version name. leave as an empty character array
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
par.fxhatload = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/optPonly_CTL_He_P_xhat.mat';

% P model parameters
par.opt_sigma = off ;
par.opt_kP_T  = off ;
par.opt_kdP   = off ; %
par.opt_bP_T  = off ;
par.opt_bP    = off ; %
par.opt_alpha = off ; %
par.opt_beta  = off ; %
% C model parameters
par.opt_bC_T  = off ;
par.opt_bC    = on ; %
par.opt_bPC	  = off; %new variable to optimize bP and bC with the same value.
					 % no temperature dependence. must have opt_bP on and opt_bP_T, opt_bC, opt_bC_T off
par.opt_d     = on ; %
par.opt_kC_T  = off ;
par.opt_kdC   = on ; %
par.opt_R_Si  = off ;
par.opt_rR    = on ; %
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
par.opt_Q10Photo     = on ; %
par.opt_fStorage     = on ; %
par.opt_fRibE 	     = off ;
par.opt_kST0 	     = on ; %
par.opt_PLip_PCutoff = off ;
par.opt_PLip_scale   = off ;
par.opt_PStor_rCutoff = on; %
par.opt_PStor_scale  = off ;
par.opt_alphaS       = on ; %
par.opt_gammaDNA	 = off ;
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
	%output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/', GridVer);
	output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/');
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
    par = SetPar(par) ;

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
