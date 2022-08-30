clc; clear all; close all 
global iter 
iter = 0 ;
on   = true  ;
off  = false ;
format short
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

Gtest = off ; 
Htest = off ;
par.optim   = on ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.LoadOpt = on  ; % if load optimial par. 
par.pscale  = 0.0 ;
par.cscale  = 1.0 ; % factor to weigh DOC in the objective function

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
fprintf('Check Point 1\n');
keyboard
% save results 
% ATTENTION: Change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('../DATA/MSK%2d/',GridVer); 
elseif isunix
    output_dir = sprintf('~/rDOC-OP/MSK%2d/',GridVer) ;
end
if ~isdir(output_dir)
    command = strcat('mkdir', " ", output_dir) ;
    system(command) ;
end

VER = strcat(output_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_nppTest_spatialGamma');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
        base_name = strcat(VER,'_PCv1');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
        base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_MODIS_CbPM_aveTeu_diffSig_O2C_uniEta');
        % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta');
        % base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta');
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
    
    GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); data.PIC(iwet); ...
           data.DIC(iwet); data.POC(iwet); data.POC(iwet);];
    GC  = real(GC) + 1e-6*randn(7*nwet,1) ;
end 
if par.Omodel == on 
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end 

%--------------------- prepare parameters ------------------
if par.optim == on 
    % load optimal parameters from a file or set them to default values 
    par = SetPar(par)  ;
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
                       'TolX',1e-6               , ...
                       'TolFun',1e-6             , ...
                       'DerivativeCheck','off'   , ...
                       'FinDiffType','central'   , ...
                       'PrecondBandWidth',Inf)   ;
%
nip = length(x0);
if(Gtest);
    dx = sqrt(-1)*eps.^3*eye(nip);
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
quit
fprintf('-------------- end! ---------------\n');