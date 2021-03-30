clc; clear all; close all
global iter
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

Gtest = off ;
Htest = off ;
par.optim   = off ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.LoadOpt = on ; % if load optimial par. 
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigma = on ; 
par.opt_kP_T  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ; 
par.opt_bP    = on ;
par.opt_alpha = on ;
par.opt_beta  = on ;
% C model parameters
par.opt_bC_T  = on ;
par.opt_bC    = on ; 
par.opt_d     = on ;
par.opt_kC_T  = on ;
par.opt_kdC   = on ; 
par.opt_R_Si  = on ; 
par.opt_rR    = on ; 
par.opt_cc    = on ;
par.opt_dd    = on ;
% O model parameters
par.opt_O2C_T = on ;
par.opt_rO2C  = on ;
par.opt_O2P_T = on ; 
par.opt_rO2P  = on ; 
% Si model parameters
par.opt_dsi   = off  ;
par.opt_at    = off ;
par.opt_bt    = off  ;
par.opt_aa    = off  ;
par.opt_bb    = off  ;
%
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
    %output_dir = sprintf(['/DFS-L/DATA/primeau/fprimeau/' ...
    %                    'TempSensi/MSK91/Zscore/'], GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
                        % 'MSK%2d/'],GridVer);
    output_dir = sprintf('/DFS-L/DATA/primeau/fprimeau/TEST_SPRING_2021/MSK%2d/', GridVer);
end
VER = strcat(output_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_Pv1');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
        base_name = strcat(VER,'_PCv2');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
        base_name = strcat(VER,'_PCOv2');
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
else
  fprintf('No optimization requested. Just run the model...\n');
  % load optimal parameters from a file or set them to default values 
  par = SetPar(par) ;
  % pack parameters into an array, assign them corresponding indices.
  par = PackPar(par) ;

  %{
  x = par.p0;
  % run the P-cycle
  [par, P, Px, Pxx] = eqPcycle(x, par) ;
  DIP = M3d+nan  ;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
  POP = M3d+nan  ;  POP(iwet) = P(1+1*nwet:2*nwet) ;
  DOP = M3d+nan  ;  DOP(iwet) = P(1+2*nwet:3*nwet) ;

  par.DIP = DIP(iwet);
  par.DOP = DOP(iwet);
  par.POP = POP(iwet);
  
  % run the C-cycle
  getfluxes = true;
  [par, C, Cx, Cxx, cfluxes] = eqCcycle(x, par, getfluxes) ;
  save('cfluxes.mat','cfluxes');
  DIC = M3d+nan ;  DIC(iwet) = C(0*nwet+1:1*nwet) ;
  POC = M3d+nan ;  POC(iwet) = C(1*nwet+1:2*nwet) ;
  DOC = M3d+nan ;  DOC(iwet) = C(2*nwet+1:3*nwet) ;
  PIC = M3d+nan ;  PIC(iwet) = C(3*nwet+1:4*nwet) ;
  ALK = M3d+nan ;  ALK(iwet) = C(4*nwet+1:5*nwet) ;
  par.DIC = DIC(iwet);
  par.ALK = ALK(iwet);
  par.POC = POC(iwet);
  par.DOC = DOC(iwet);
  
  % run the O2-cycle
  [par, O, Ox, Oxx] = eqOcycle(x, par) ;
  O2 = M3d+nan ;  O2(iwet) = O ;
  %}
  load cfluxes.mat
  
  J1 = cfluxes.J1;
  J2 = cfluxes.J2;
  J3 = cfluxes.J3;
  J4 = cfluxes.J4;

  % setup the operator for trapezoid rule transport integration
  nstep_per_year = 12;
  sec_per_year = 365.25*24*60*60;                  % (s/yr)
  dt = sec_per_year/nstep_per_year;                % (s)
  T = par.year(1):dt/(sec_per_year):par.year(end); % (yrs)
  [m,n] = size(par.TRdiv);
  I = speye(n);
  A = I+(dt/2)*par.TRdiv;
  B = I-(dt/2)*par.TRdiv;  
  FA = mfactor(A);
  %
  dic = zeros(n,length(T)); 
  load(par.fname);
  par.DIC = data.DIC(iwet)-par.dicant(iwet);
  par.ALK = data.ALK(iwet);

  % check that we are starting at an equilibrium point
  vout = Fsea2air(par,'CO2');
  JgDIC = vout.JgDIC;
  %check the the initial condition is in steady state
  %eq1 = A*par.DIC -B*par.DIC + dt*(J1 + J2 + J3 + J4 - JgDIC);
  % norm(eq1) (should be on the order of 1e-7)

  dic(:,1) = par.DIC;     % DIC in 1750
  dic_tot = zeros(length(T),1);
  dic_tot(1) = sum(dic(:,1).*par.dVt(iwet))/1000; % moles of DIC
  N = length(T);
  for k = 2:N
    t = T(k); % time (years)
    par.pco2atm = interp1(par.year,par.pco2_air,t);
    vout = Fsea2air(par,'CO2');
    JgDIC = vout.JgDIC ; % 
    rhs = B*dic(:,k-1) + ...
          dt*(-J1-J2-J3-J4+JgDIC);
    dic(:,k) = mfactor(FA,rhs);
    par.DIC = dic(:,k);

    dic_tot(k) = (sum(dic(:,k).*par.dVt(iwet)))/1000; % moles of DIC
    
    if (mod(k,12)==0)
      % make a plot
      plot(T(1:k),12*(dic_tot(1:k)-dic_tot(1))/(1e15),'-o'); 
      set(gca,'FontSize',16);
      grid on
      xlabel('time (years)');
      ylabel('anthro dic (Pg)');
      drawnow
    end
    
  end
  % save the output to a file
  fname = strcat(output_dir,'transient_co2.mat');
  save(fname,'dic','T','-v7.3');

end
fprintf('-------------- end! ---------------\n');