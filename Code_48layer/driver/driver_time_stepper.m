%
%
% disp(sprintf('loading steady state from %s',par.output_dir));
% load([par.output_dir par.output_eq],'par');
%
% load optimized parameters and solutions from PCO model --> get steady state solutions for C13 and C14 --> do time stepping method.
clc; clear all;
addpath('../../DATA/BGC_48layer/')
addpath('../src/MSK91/')
addpath('../src/')

%
global iter GC GC13 GC14 GO;
iter = 0 ;
on   = true  ;
off  = false ;
format short
% 
GridVer  = 91  ;
operator = 'L' ;
% GridVer: choose from 90 and 91; Ver 90 is for a Transport
% operator without diapycnal mixing but optimized using DIP ;
% Ver 91 include a bunch of operators that include diapycnal
% mixing. These operators represent sensiviity tests on He
% constraint and on mixing parameterizations (DeVries et al, 2018).
% A -> CTL_He; B -> CTL_noHe; C -> KiHIGH_He; D -> KiHIGH_noHe;
% E -> KvHIGH_KiLOW_He; F -> KvHIGH_KiLOW_noHe; G -> KiLOW_He;
% H -> KiLOW_noHe; I -> KvHIGH_He; J -> KvHIGH_noHe; K -> KvHIGH_KiHIGH_noHe
par.optim     = off ; 
par.Pmodel    = on  ; 
par.Cmodel    = on  ; 
par.Omodel    = on  ; 
par.Simodel   = off ;
par.Cisotope  = on  ;
par.LoadOpt   = on  ; % if load optimial par. 

par.dopscale = 0.0 ;
par.dipscale = 1.0 ;
par.dicscale = 1.0 ;
par.docscale = 1.0 ; % factor to weigh DOC in the objective function
par.alkscale = 1.0 ;
par.o2scale  = 1.0 ;

% par.xxxscale ------> 실제 코드에 쓰이지는 않음.. 그러나, Packpar를 위해서는 Cmodel = on이고, Cmodel = on이면 위의 scale값 있어야..
% Packpar가 꼭 쓰여야 하는건지 체크 필요.

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
par.opt_O2C_T = off ;
par.opt_rO2C  = on ;
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;

%--------------------------------------
SetUp;

output_dir = sprintf('../Results/'); 
VER = strcat(output_dir);
fopt = 'CTL_He_48layer_PCO_firstrun';
par.fname = strcat(VER, fopt,'.mat');
par.fxhat = strcat(VER, fopt, '_xhat.mat');
par.optim = off;

%--------------------- prepare parameters ------------------
% load optimal parameters and data from a file or set them to default values 
load(par.fname, 'data')
par = SetPar(par)  ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

x0    = par.p0 ; % set up parameters and assign it to x0

% check the optimzed data and solutions
mf = delta1314();
par.pco2atm = 276.8; % 
par.c13.R13a = mf.d2r13(-6.61); % 
par.c14.R14a = mf.D2r14(0,-6.61); % air C14/C at 1750


% load /DFS-L/DATA/primeau/oceandata/DATA/eSST_obs_1749_2023.mat;
% par.Temp_eSST = Temp_eSST;
% par.Temp_eSST_tt = Temp_eSST_tt;
% par.Temp_obs  = par.Temp; % keep a copy of WOA temp in case it is updated by temperature model
par.saveall = true;

% CAUTION: check these parameters before running an experiment
%------------------------------------------------------------------------
par.yst        = 1750;
par.yed        = 2022;
par.fras       = 1.00;
par.frpho      = 0.5 ;  %0.5
par.fc14       = 1.9;
% par.fkw        = 0.72; %0.72
% load('kw_ocim_10122023.mat');
% par.kw=kw.kmean2D_af/(100*60^2); % convert from cm/hr to m/s;
% par.kw_ref     = par.kw;
%------------------------------------------
% par.kw         = par.kw*0.225/0.337;
% the air-sea gas exchange is modified here
%------------------------------------------
% par.kw         = par.kw*par.fkw;
% par.TRdiv      = par.TRdiv*0.5;
% 

fprintf('solve the steady state of P...\n')
tic
Ppool = {'DIP','POP','DOP','DOPl'};
[par,data] = solveM.eqP(x0,par,data,Ppool); 
toc

fprintf('solve the steady state of C...\n')
tic
Cpool = {'DIC','POC','DOC','PIC','ALK','DOCl','DOCr'};
[par,data] = solveM.eqC(x0,par,data,Cpool); 
toc

fprintf('solve the steady state of C13...\n')
tic
C13pool = {'DIC13','POC13','DOC13','PIC13','DOC13l','DOC13r'};
[par,data] = solveM.eqC13(x0,par,data,C13pool); 
toc

fprintf('solve the steady state of C14...\n')
tic
C14pool = {'DIC14','POC14','DOC14','PIC14','DOC14l','DOC14r'};
[par,data] = solveM.eqC14(x0,par,data,C14pool); 
toc

fprintf('solve the steady state of O2...\n')
tic
O2pool = {'O2'};
[par,data] = solveM.eqO2(x0,par,data,O2pool);
toc

% save all par and data before the transient run
fileName  = 'SSdata_par_0.707xkw_fras=1_frpho=0.5_48layer_240426'
directory = '../Results'
filePath  = fullfile(directory, fileName) ;
save(filePath, 'par', 'data', '-v7.3')   ;

% check the data and par for steady-state
% drive the time stepping method for C, C13, C14, and O2
Ctype = {'C12','C13','C14','O2'}; 
for m = 1:length(Ctype)
  vn = Ctype{m};
  if contains(vn,'12')
    Xin.C12 = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr];
  elseif contains(vn, 'O2')
    Xin.O2 = [par.O2];
  else
    eval(sprintf('Xin.%s = [par.DI%s;par.PO%s;par.DO%s;par.PI%s;par.DO%sl;par.DO%sr];'...
                        ,vn,vn,vn,vn,vn,vn,vn));
  end
end
  

t0 = 1750; t1 = 2022;
fprintf('Time stepping method for C isotopes and O2')
tic
[Xout,Tout,par] = time_stepper(par,t0,t1,Xin,Ctype);
toc



