% driver_optC_Cell.m
%   Driver script to run biogeochemical model.
%   The carbon cycle uses a trait-based celllular growth model to compute C:P.
%   This run optimizes phosphorus, carbon, and oxygen cycle parameters
% ------------------------------------------------------------------------
clc; clear all; close all
global iter
iter = 0 ;
on   = true  ; off  = false ;
format short

% --- addpath to model code -----
addpath('../src/')

VerName = 'optPCO_Tz_v2_'; 		% optional version name. leave as an empty character array
					% or add a name ending with an underscore
VerNum = '';		% optional version number for testing

% Choose C2P function
par.C2Pfunctiontype = 'T';
% 'P' -> PO4 function ; 'C' -> Cell model; 'T' -> Temperature function; 'R' -> constant value (Redfield)
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
% L -> CTL_He_48layer

par.nl = 2; % number of layers in the model euphotic zone (doesn't change)

Gtest = off ; 
Htest = off ;
par.optim   = on ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.Cisotope  = off  ;
par.LoadOpt = on ; % if load optimial parameters. 
% to load parameter values from a run with a different name.
%par.fxhatload = '../../output/optPonly_CTL_He_P_xhat.mat';
%par.fxhatload = '../output/optPCO_GM15_CTL_He_PCOv1_DOC1_DOP0_xhat.mat';
% par.fxhatload = '/DFS-L/DATA/primeau/hojons1/Nature2023_BGC_reoptimized/src_Nature_parameter_Megan/MSK91/CTL_He_PCO_Gamma0_kl12h_O5_POC2DIC_GM15_Nowicki_npp1_aveTeu_diffSig_O2C_uniEta_DICrmAnthro_2L_Pnormal_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_xhat.mat' 
par.dynamicP = off ; % if on, cell model uses modeled DIP. if off, cell model uses WOA observed DIP field.

par.dopscale = 0.0 ;
par.dipscale = 1.0 ;
par.dicscale = 1.0 ;
par.docscale = 1.0 ; % factor to weigh DOC in the objective function
par.alkscale = 1.0 ;
par.o2scale  = 1.0 ;
% P model parameters
par.opt_sigP  = on ; 
par.opt_Q10P  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ; 
par.opt_bP    = on ;
par.opt_alpha = on ;
par.opt_beta  = on ;
% C model parameter
par.opt_sigC  = on ; 
par.opt_kru   = on ;
par.opt_krd   = on ;
par.opt_etau  = on ;
par.opt_etad  = off ; %keep off
par.opt_bC_T  = on ;
par.opt_bC    = on ; 
par.opt_d     = on ;
par.opt_Q10C  = on ;
par.opt_kdC   = on ; 
par.opt_R_Si  = on ; 
par.opt_rR    = on ; 
% --- C:P function parameters -----
% phosphate-dependent function parameters
par.opt_cc    = off ;
par.opt_dd    = off ; 
% temperature-dependent function parameters
par.opt_ccT   = on ; 
par.opt_ddT   = on ;
% Trait-based Cellular Growth Model parameters
par.opt_Q10Photo     = on ; % opt
par.opt_fStorage     = on ; % opt
par.opt_fRibE 	     = off ; 
par.opt_kST0 	     = on ; % opt
par.opt_PLip_PCutoff = off ;
par.opt_PLip_scale   = off ;
par.opt_PStor_rCutoff = on; % opt
par.opt_PStor_scale  = off ;
par.opt_alphaS       = on ; % opt
par.opt_gammaDNA	 = off ;
% O model parameters
par.opt_O2C_T = off ;
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
% ATTENTION: Change this directory to where you want to
% save your output files
output_dir = sprintf('../output/'); 

if ~isdir(output_dir)
    command = strcat('mkdir', " ", output_dir) ;
    system(command) ;
end

VER = strcat(output_dir,VerName,TRdivVer);
catDOC = sprintf('_DOC%0.2g_DOP%0.2g',par.docscale,par.dopscale); % used to add scale factors to file names
% Create output file names based on which model(s) is(are) optimized
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
% -------------------- Set up output files ---------------
par.fname = strcat(fname,'.mat') ;
fxhat     = strcat(fname,'_xhat.mat') 
fxpar     = strcat(fname,'_par.mat');
par.fxhat = fxhat ;
if Htest ==on
	fGHtest = strcat(fname,'_GHtest.mat')  ;
end
par.fxhat = fxhat ;
par.fxpar = fxpar ;

% -------------------update initial guesses --------------
if isfile(par.fname)
    load(par.fname)
end 

%---------------- inital guesses on C and O ---------------
if par.Cmodel == on 
    GC  = [data.DIC(iwet); data.POC(iwet); data.DOC(iwet); data.PIC(iwet); ...
           data.ALK(iwet); data.DOC(iwet); data.DOC(iwet);];
    GC  = real(GC) + 1e-6*randn(7*nwet,1) ;
end 
if par.Omodel == on 
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end 

%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par)  ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

%-------------------set up fminunc -------------------------
x0    = par.p0 ;
myfun = @(x) neglogpost(x, par);
objfuntolerance = 5e-11; %5e-12;
options = optimoptions(@fminunc                  , ...
                       'Algorithm','trust-region', ...
                       'GradObj','on'            , ...
                       'Hessian','on'            , ...
                       'Display','iter'          , ...
                       'MaxFunEvals',2000        , ...
                       'MaxIter',2000            , ...
                       'TolX', objfuntolerance   , ...     % 기존은 -7. decreasing
                       'TolFun',objfuntolerance  , ...     % 기존은 -7. decreasing
                       'DerivativeCheck','off'   , ...
                       'FinDiffType','central'   , ...
                       'PrecondBandWidth',Inf)   ;
%
nip = length(x0);
if(Gtest);
    fprintf('derivative test \n')
	GHtest.fx_cstep = NaN([nip,1]);
	GHtest.fxx_cstep = NaN(nip);
	GHtest.pindx = par.pindx;
    % load(fGHtest); % continuing GHtest after out of memory crash
    % display par.pindx
    fprintf('pindx = ')
    par.pindx 
    dx = sqrt(-1)*eps.^3*eye(nip);
    for ii = 1 : nip 
        x  = real(x0)+dx(:,ii);
        iter = 11; %bypasses the ResetPar in neglogpost
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
    % run model once more to save gradient and hessian in xhat without the complex step
    fprintf('Run without compex step & Save f,fx,fxx to xhat \n')
    [f,fx,fxx,data,xhat] = neglogpost(x0,par);
    xhat.pindx = par.pindx;
    xhat.f   = f   ;
    xhat.fx  = fx  ;
    xhat.fxx = fxx ;
    % save results 
    fprintf('saving optimized parameters to file: %s \n',fxhat)
    save(fxhat, 'xhat')
elseif (par.optim)
    [xsol,fval,exitflag] = fminunc(myfun,x0,options);
    fprintf('objective function tolerance = %5.1e \n',objfuntolerance);
    fprintf('----fminunc complete----\n')
    [f,fx,fxx,data,xhat] = neglogpost(xsol,par);
    fprintf('----neglogpost solved for final parameter values----\n')
    xhat.pindx = par.pindx;
    xhat.f   = f   ;
    xhat.fx  = fx  ;
    xhat.fxx = fxx ;
    % save results 
    fprintf('saving optimized parameters to file: %s \n',fxhat)
    fprintf('saving model solution to file: %s \n',par.fname)
    save(fxhat, 'xhat')
    save(par.fname, 'data')
else
    xsol = x0;
    iter = 11;
    [f,fx,fxx,data] = neglogpost(xsol,par);
    fprintf('----neglogpost complete----\n')
    %% note: skipping save for testing
    %fprintf('saving model solution to file: %s \n',par.fname)
    %save(par.fname, 'data')
end
fprintf('-------------- end! ---------------\n');
quit
