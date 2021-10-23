clc; clear all; close all
global iter
%global fmin_history

% set up parallel pool
if isempty(gcp('nocreate'))
	poolobj = parpool(16) % if running in parallel; comment out to run serial
	%to close the parallel pool: delete(poolobj)
	poolobj.IdleTimeout = 120;
end

iter = 0 ;
on   = true  ;
off  = false ;
format long
%
GridVer  = 90  ;
operator = 'A' ;
% GridVer: choose from 90 and 91; Ver 90 is for a Transport
% operator without diapycnal mixing but optimized using DIP ;
% Ver 91 include a bunch of operators that include diapycnal
% mixing. These operators represent sensiviity tests on He
% constraint and on mixing parameterizations (DeVries et al, 2018).
% A -> CTL_He; B -> CTL_noHe; C -> KiHIGH_He; D -> KiHIGH_noHe;
% E -> KvHIGH_KiLOW_He; F -> KvHIGH_KiLOW_noHe; G -> KiLOW_He;
% H -> KiLOW_noHe; I -> KvHIGH_He; J -> KvHIGH_noHe; K -> KvHIGH_KiHIGH_noHe

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
par.opt_sigma = off ;
par.opt_kP_T  = off ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ;
par.opt_bP    = on ;
par.opt_alpha = on ;
par.opt_beta  = on ;
% C model parameters
par.opt_bC_T  = on ;
par.opt_bC    = on ;
par.opt_d     = on ;
par.opt_kC_T  = off ;
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
par.opt_fStorage     = on ;
par.opt_PLip_PCutoff = on ;
par.opt_PLip_scale   = off ;
par.opt_PStor_rCutoff = on;
par.opt_PStor_scale  = off ;
par.opt_alphaS       = off ;
par.opt_fRibE 	     = off ;
par.opt_kST0 	     = on ;
% par.BIO.opt_gammaDNA = off;
% par.BIO.opt_gammaLipid = off;
% par.BIO.opt_DNT0 = off;
% par.BIO.opt_DPT0 = off;
% par.BIO.opt_Q10Diffusivity = off;
% par.BIO.opt_AMin = off;
% par.BIO.opt_PhiS = off;


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
	output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/', GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
                        % 'MSK%2d/'],GridVer);
end
VER = strcat(output_dir,TRdivVer);
catDOC = sprintf('_DOC%0.2g_DOP%0.2g',par.cscale,par.pscale); % used to add scale factors to file names
% Creat output file names based on which model(s) is(are) optimized
%if Gtest == on
%    fname = strcat(VER,'_GHtest');
%elseif Gtest == off
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
        base_name = strcat(VER,'_PCCellv6b');
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == on)
		base_name = strcat(VER,'_PCOCell');
		fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == on)
        base_name = strcat(VER,'_PCOSiCell');
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

%par.fhistory = strcat(fname,'_history.mat');

% -------------------update initial guesses --------------
if isfile(par.fname)
    load(par.fname)
end

%---------------- inital guesses on C and O ---------------
if par.Cmodel == on
    DIC = data.DIC - par.dicant ;

    GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); ...
           data.PIC(iwet); data.ALK(iwet)];

	% temprarily use following line to create ALK initial field for 90x180 grid
	% GC = [DIC(iwet); data.POC(iwet); data.DOC(iwet); data.PIC(iwet); data.DIC(iwet)];

    GC  = GC + 1e-6*randn(5*nwet,1) ;
end
if par.Omodel == on
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end

%--------------------- prepare parameters ------------------
if (par.optim == on) | (par.LoadOpt == on)
    % load optimal parameters from a file or set them to default values
    par = SetPar(par) ;
    % pack parameters into an array, assign them corresponding indices.
    par = PackPar(par) ;
end

%-------------------set up fminunc -------------------------
x0    = par.p0 ;
myfun = @(x) neglogpost(x, par);

%fmin_history = struct;
% function stop = fminoutfun(x, optimValues, state)
% 	%global iter
% 	%global fmin_history
% 	stop = false;
% 	switch state
%   		case 'init'
% 			fmin_history = struct();
% 		    % fmin_history.x = x;
% 		    % fmin_history.optimValues = optimValues;
% 		    % fmin_history.timerVal = tic;
% 			fmin_history.x = zeros(100,length(x));
% 			fmin_history.gradient = zeros(100,length(x));
% 			fmin_history.fval = zeros(100,1);
% 			fmin_history.timerVal = zeros(100,1);
% 			fmin_history.timerStart = tic;
% 			fprintf('iter value at init state = %i',iter)
% 		case 'iter'
% 			% ind = length(fmin_history)+1;
% 		    % fmin_history(iter).x = x;
% 		    % fmin_history(iter).optimValues = optimValues;
% 		    % fmin_history(iter).timerVal = toc(fmin_history(1).timerVal);
% 			fmin_history.x(iter,:) = x;
% 			fmin_history.fval(iter) = optimValues.fval;
% 			fmin_history.gradient(iter,:) = optimValues.gradient;
% 			fmin_history.timerVal = toc(fmin_history.timerStart);
% 		case 'done'
% 			%
%   	end
% 	%history.x = [history.x; x];
% 	%history.fval = [history.fval; optimValues.fval];
% 	%history.gradient = [history.gradient; optimValues.gradient'];
% 	%save(par.fhistory,'history') %try making history an output of fminunc after exitflag, then only need to save once
% end
options = optimoptions(@fminunc                  , ...
                       'Algorithm','trust-region', ...
                       'GradObj','on'            , ...
                       'Hessian','on'            , ...
                       'Display','iter'          , ...
                       'MaxFunEvals',2000        , ...
                       'MaxIter',2000            , ...
                       'TolX',5e-7               , ...
                       'TolFun',5e-7             , ...
                       'DerivativeCheck','off'   , ...
                       'FinDiffType','central'   , ...
                       'PrecondBandWidth',Inf    , ...
					   'SubproblemAlgorithm','factorization') ; %testing this
					   %'OutputFcn',@fminoutfun ) 	 ; % maybe use 'SubproblemAlgorithm','factorization' -changing subproblem did not seem to effect anything
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
            fprintf('%i % .3e  \n',ii,diff);
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
    [xsol,fval,exitflag] = fminunc(myfun,x0,options);
	fprintf('----fminunc complete----\n')
    [f,fx,fxx,data] = neglogpost(xsol,par);
	fprintf('----neglogpost solved for final parameter values----\n')
    load(fxhat,'xhat')
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
end
delete(gcp('nocreate')) % shut down parallel pool
fprintf('-------------- end! ---------------\n');
