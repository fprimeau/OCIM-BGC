clc; clear all; close all 
global iter GC GC13;
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
par.optim     = on ; 
par.Pmodel    = on ; 
par.Cmodel    = on ; 
par.C13model  = on;
par.Omodel    = on ; 
par.Simodel   = off ;
par.LoadOpt   = on  ; % if load optimial par. 
par.pscale    = 0.0 ;
par.cscale    = 1.0 ; % factor to weigh DOC in the objective function

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
fprintf('Check Setup: struct data,par,etc \n');
% save results 
% ATTENTION: Change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('../DATA/MSK%2d/',GridVer); 
elseif isunix
    output_dir = sprintf('/DFS-L/DATA/primeau/oceandata/OCIM_BGC_DATA/DATA/MSK%2d/',GridVer);
    output_dir = sprintf('/DFS-L/DATA/primeau/oceandata/ModelOutput/MSK%2d/',GridVer) ;
end
if ~isdir(output_dir)
    command = strcat('mkdir', " ", output_dir) ;
    system(command) ;
end

VER = strcat(output_dir,TRdivVer);
% Create output file names based on which model(s) is(are) optimized
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
    if  ( par.C13model == on)
      fname = strcat(base_name,'_C13')
    end
end
% -------------------- Set up output files ---------------
par.fname = strcat(fname,'.mat') ; 
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ; 

% -------------------update initial guesses --------------
if isfile(par.fname)
  disp(sprintf('load optimized solution %s',par.fname));
  load(par.fname)
  disp(sprintf('xhat file exist: %s',string(logical(isfile(par.fxhat)))));
end 

%---------------- inital guesses on C, C13 and O ---------------
if par.Cmodel == on 
    DIC = data.DIC - par.dicant ;
    
    GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); data.PIC(iwet); ...
           data.DIC(iwet); data.POC(iwet); data.POC(iwet);];
    GC  = real(GC) + 1e-6*randn(7*nwet,1) ;
end 

if par.C13model == on
  GC13 = zeros(6*length(iwet),1);
end

if par.Omodel == on 
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end 

%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par)  ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

x0    = par.p0 ; % set up parameters and assign it to x0
par.skipit = on;
PCsol = sprintf('%stemp.mat',output_dir);
if ~par.skipit 

  %-------------------solve P model -------------------------
  if par.Pmodel == on 
    [par, P] = eqPcycle(x0, par) ;
    % [par, P, Px, Pxx] = eqPcycle(x0, par) ;
    % Gradient and Hessian
    % par.Px    = Px ;  par.Pxx   = Pxx ;
    
    DIP  = M3d+nan  ;  DIP(iwet)  = P(1+0*nwet:1*nwet) ;
    POP  = M3d+nan  ;  POP(iwet)  = P(1+1*nwet:2*nwet) ;
    DOP  = M3d+nan  ;  DOP(iwet)  = P(1+2*nwet:3*nwet) ;
    DOPl = M3d+nan  ;  DOPl(iwet) = P(1+2*nwet:3*nwet) ;
    
    par.DIP  = DIP(iwet) ;
    data.DIP = DIP ; data.POP  = POP  ;
    data.DOP = DOP ; data.DOPl = DOPl ;
  end
  %-------------------solve C model -------------------------
  if par.Cmodel == on 
    [par, C] = eqCcycle(x0, par) ;
    % [par, C, Cx, Cxx] = eqCcycle(x0, par) ;
    % Gradient and Hessian
    % par.Cx = Cx ;  par.Cxx = Cxx ;
    
    DIC  = M3d+nan ;  DIC(iwet)  = C(0*nwet+1:1*nwet) ;
    POC  = M3d+nan ;  POC(iwet)  = C(1*nwet+1:2*nwet) ;
    DOC  = M3d+nan ;  DOC(iwet)  = C(2*nwet+1:3*nwet) ;
    PIC  = M3d+nan ;  PIC(iwet)  = C(3*nwet+1:4*nwet) ;
    ALK  = M3d+nan ;  ALK(iwet)  = C(4*nwet+1:5*nwet) ;
    DOCl = M3d+nan ;  DOCl(iwet) = C(5*nwet+1:6*nwet) ;
    DOCr = M3d+nan ;  DOCr(iwet) = C(6*nwet+1:7*nwet) ;
    
    par.DIC = DIC(iwet) ;
    par.POC = POC(iwet) ;
    par.DOC = DOC(iwet) ;
    par.PIC = PIC(iwet);
    par.ALK = ALK(iwet) ;
    par.DOCl = DOCl(iwet) ;
    par.DOCr = DOCr(iwet) ;
    % DIC = DIC + par.dicant  ;
    
    data.DIC  = DIC  ;  data.POC  = POC ;
    data.DOC  = DOC  ;  data.PIC  = PIC ;
    data.ALK  = ALK  ;  data.DOCr = DOCr ;
    data.DOCl = DOCl ;
  end
  save(PCsol,'data','par','-v7.3');
else
  disp(sprintf('loading P and C solution from %s',PCsol));
  load(PCsol);
end

%-------------------solve C13 model -------------------------
if par.C13model == on 
  par.debug13 = off;
  [par, C13 ] = eqC13cycle(x0, par);
  % Gradient and Hessian
  %par.C13x = C13x ;  par.C13xx = C13xx ;
  
  DIC13  = M3d+nan ;  DIC13(iwet)  = C13(0*nwet+1:1*nwet) ;
  POC13  = M3d+nan ;  POC13(iwet)  = C13(1*nwet+1:2*nwet) ;
  DOC13  = M3d+nan ;  DOC13(iwet)  = C13(2*nwet+1:3*nwet) ;
  PIC13  = M3d+nan ;  PIC13(iwet)  = C13(3*nwet+1:4*nwet) ;
  DOC13l = M3d+nan ;  DOC13l(iwet) = C13(4*nwet+1:5*nwet) ;
  DOC13r = M3d+nan ;  DOC13r(iwet) = C13(5*nwet+1:6*nwet) ;

  par.DIC13 = DIC13(iwet) ;
  par.POC13 = POC13(iwet) ;
  par.DOC13 = DOC13(iwet) ;
  par.DOC13l = DOC13l(iwet) ;
  par.DOC13r = DOC13r(iwet) ;
  % DIC = DIC + par.dicant  ;

  data.DIC13  = DIC13  ;  data.POC13  = POC13 ;
  data.DOC13  = DOC13  ;  data.PIC13  = PIC13 ;
  data.DOC13r = DOC13r ;  data.DOC13l = DOC13l ;
end
par.d13 = @(dic13,dic) (dic13./dic/1.12372*1e2-1)*1e3;
save(par.fname, 'data','par','-v7.3');
