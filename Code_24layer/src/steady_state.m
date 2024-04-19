% a script for testing different carbon and isotopes 
clc; clear all;
global iter GC GC13 GC14 GO;
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
par.optim     = off ; 
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

if ismac
    output_dir = sprintf('../DATA/MSK%2d/',GridVer); 
elseif isunix
    %output_dir=sprintf('/DFS-L/DATA/primeau/salali/ModelOutput/MSK%2d/',GridVer);
   output_dir = sprintf('/DFS-L/DATA/primeau/oceandata/OCIM_BGC_DATA/DATA/MSK%2d/',GridVer);
    % output_dir = sprintf('/DFS-L/DATA/primeau/oceandata/ModelOutput/MSK%2d/',GridVer) ;
end
if ~isdir(output_dir)
    command = strcat('mkdir', " ", output_dir) ;
    system(command) ;
end

% set up model parameters and load data
% SetUp ;
% fprintf('Check Setup: struct data,par,etc \n');

TRdivVer = 'CTL_He';
VER = strcat(output_dir,TRdivVer);

fopt = '_Gamma1to3_POC2DIC_GM15_MODIS_CbPM_aveTeu_diffSig_O2C_uniEta_DOC1e+00_DOP0e+00';
par.fname = strcat(VER,'_PCO',fopt,'.mat');

%keyboard
if isfile(par.fname)
 %keyboard
  load(par.fname,'data','par');
  disp('data are updated from optimzed solutions');
  disp(sprintf('xhat file exist: %s',string(logical(isfile(par.fxhat)))));
end
par.optim = off;

par.fxhat = strcat(VER,'_PCO',fopt,'_xhat.mat');

%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par)  ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

x0    = par.p0 ; % set up parameters and assign it to x0

mf = delta1314();
par.pco2atm = 276.8; % 
par.c13.R13a = mf.d2r13(-6.61); % 
par.c14.R14a = mf.D2r14(0,-6.61); % air C14/C at 1750

% solve the P-cycle model
Ppool = {'DIP','POP','DOP','DOPl'};
% Psol = sprintf('%stempP_%i_Ppools.mat',output_dir,length(Ppool));
% [par,data] = solveM.eqP(x0,par,data,Ppool); 

%------------------------------------------------------------------------
% CAUTION: check these parameters before running an experiment
%------------------------------------------------------------------------
par.yst        = 1750;
par.yed        = 2022;
par.frpho      = 0.5;  %0.5
par.fras       = 1.00;
par.fc14       = 2.0;
par.fkw        = 0.72; %0.72
load('kw_ocim_10122023.mat');
par.kw=kw.kmean2D_af/(100*60^2); % convert from cm/hr to m/s;
par.kw_ref     = par.kw;
%------------------------------------------
% par.kw         = par.kw*0.225/0.337;
% the air-sea gas exchange is modified here
%------------------------------------------
par.kw         = par.kw*par.fkw;
% par.TRdiv      = par.TRdiv*0.5;
% modify the output dir here
%------------------------------------------
par.output_dir = sprintf('/DFS-L/DATA/primeau/salali/ModelOutput/tmp/C_%i_%i_c14aG_frac_%3.1f_kw_%4.2f_fras_%4.2f_D/',par.yst,par.yed,par.frpho,par.fkw,par.fras);
par.output_eq  = sprintf('tempPC_C13_C14_%i_frac_%3.1f_fas_%4.2f_fc14_2.0_kw_%4.2f.mat',par.yst,par.frpho,par.fras,par.fkw);
if ~isdir(par.output_dir)
  system(sprintf('mkdir %s',par.output_dir));
end
%------------------------------------------------------------------------

Cpool = {'DIC','POC','DOC','PIC','ALK','DOCl','DOCr'};
[par,data] = solveM.eqC(x0,par,data,Cpool); 


C13pool = {'DIC13','POC13','DOC13','PIC13','DOC13l','DOC13r'};
[par,data] = solveM.eqC13(x0,par,data,C13pool); 

C14pool = {'DIC14','POC14','DOC14','PIC14','DOC14l','DOC14r'};
[par,data] = solveM.eqC14(x0,par,data,C14pool); 

% GO = real(data.O2(par.iwet)) + 1e-9*randn(par.nwet,1) ;
% par.G = d0( par.alpha*par.L*par.DIP ) ;
% [par, O2] = eqOcycle(x0, par) ;
% par.O2 = O2;

% solving steady state of TT (temperature)
% par.SST_obs = nanmean(par.Temp_eSST(:,:,1:12),3);
% [f,J,par] = TTeqn([], par);
% Fp = mfactor(J);
% par.TT = mfactor(Fp,-f);

save([par.output_dir par.output_eq],'par','data','-v7.3');
