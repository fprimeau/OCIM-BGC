clc; clear all; close all
global iter
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
par.opt_kdP   = off ;
par.opt_bP_T  = off ;
par.opt_bP    = off ;
par.opt_alpha = off ;
par.opt_beta  = off ;
% C model parameters
par.opt_bC_T  = off ;
par.opt_bC    = off ;
par.opt_d     = on ;
par.opt_kC_T  = off ;
par.opt_kdC   = off ;
par.opt_R_Si  = off ;
par.opt_rR    = off ;
par.opt_cc    = off ;
par.opt_dd    = off ;
% O model parameters
par.opt_O2C_T = off ;
par.opt_rO2C  = on ;
par.opt_O2P_T = off ;
par.opt_rO2P  = on ;
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;
%Trait Model parameters
par.opt_Q10Photo = on ;
par.opt_fRibE 	 = off;
par.opt_fStorage = off;
par.opt_kST0 	 = off;
% par.BIO.opt_alphaS = off;
% par.BIO.opt_gammaDNA = off;
% par.BIO.opt_gammaLipid = off;
% par.BIO.opt_lPCutoff = off;
% par.BIO.opt_DNT0 = off;
% par.BIO.opt_DPT0 = off;
% par.BIO.opt_Q10Diffusivity = off;
% par.BIO.opt_AMin = off;
% par.BIO.opt_CStor = off;
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
	output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/test%s/', datestr(now,'ddmmmyy'));

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
        base_name = strcat(VER,'_PCCell');
        fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == on)
		base_name = strcat(VER,'_PCOCell');
		fname = strcat(base_name,catDOC);
	elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == on)
        base_name = strcat(VER,'_PCOSiCell');
        fname = strcat(base_name,catDOC);
    end
end
% -------------------- Set up output files ---------------
par.fname = strcat(fname,'.mat') ;
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ;

% -------------------update initial guesses --------------
if isfile(par.fname)
    load(par.fname)
end

%---------------- inital guesses on C and O ---------------
if par.Cmodel == on
    DIC = data.DIC - par.dicant ;

    GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); ...
           data.PIC(iwet); data.ALK(iwet)];
    GC  = GC + 1e-6*randn(5*nwet,1) ;
end
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

%-------------------set up fminunc -------------------------
x0    = par.p0 ;
myfun = @(x) neglogpost(x, par);
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
                       'PrecondBandWidth',Inf)   ;
%
nip = length(x0);
if(Gtest);
    dx = sqrt(-1)*eps.^3*eye(nip);
    % for ii = nip : -1 : 13
    for ii = 1 : nip
        x  = real(x0)+dx(:,ii);
        if Htest == on
            [f,fx,fxx] = neglogpost(x, par) ;
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
    keyboard
else
    [xhat,fval,exitflag] = fminunc(myfun,x0,options);
    [f,fx,fxx,data] = neglogpost(xhat,par);
    load(fxhat)
    xhat.f   = f   ;
    xhat.fx  = fx  ;
    xhat.fxx = fxx ;
    % save results
    save(fxhat, 'xhat')
    save(par.fname, 'data')
end

fprintf('-------------- end! ---------------\n');
