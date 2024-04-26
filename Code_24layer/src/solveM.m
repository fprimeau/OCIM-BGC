classdef solveM 
  methods(Static)
  % wrapper of solving different model equilibrium state

  function [par,data] = eqP(x0,par,data,Ppool)
  %-------------------solve P model -------------------------
    [par, P] = eqPcycle(x0, par) ;
    % [par, P, Px, Pxx] = eqPcycle(x0, par) ;
    % Gradient and Hessian
    % par.Px    = Px ;  par.Pxx   = Pxx ;

    % extract P model solution and pass it to par and data
    [par,data] = dispSol(P,par,data,Ppool);
  end

  function [par,data] = eqC(x0,par,data,Cpool)
    global GC;
    iwet = par.iwet; nwet = par.nwet;
    par.nC = length(Cpool);
    % make sure the ALK is the same in C and C13 model
    % par.ALK = data.ALK(iwet);
    % par.ALK(isrf) = par.co2syspar.alk;  % set GLODAP surface ALK for Fsea2air
    %---------------- inital guesses on C, C13 and O ---------------
    % DIC = data.DIC - par.dicant ;
    % DIC = data.DIC(iwet) ;

    % GC = [DIC(iwet); data.PIC(iwet);];
    % GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); data.PIC(iwet); ...
    %        data.DIC(iwet); data.POC(iwet); data.POC(iwet);];

    % the number of C pools is adjustable with the following 
    s2e = 'data.DIC(iwet);';
    for i = 2:length(Cpool)
      s2e = strcat(s2e,sprintf('data.%s(iwet);',Cpool{i}));
    end
    disp(sprintf('initial: GC = %s',s2e));
    GC = eval(sprintf('[%s]',s2e));
    % keyboard
    
    GC  = real(GC) + 1e-6*randn(par.nC*nwet,1) ;

    %-------------------solve C model -------------------------
    [par, C] = eqCcycle_v2(x0, par) ;
    % [par, C, Cx, Cxx] = eqCcycle(x0, par) ;
    % Gradient and Hessian
    % par.Cx = Cx ;  par.Cxx = Cxx ;

    [par,data] = dispSol(C,par,data,Cpool);

    % data.DIC = data.DIC + par.dicant  ;
    % save(PCsol,'data','par','-v7.3');
  end %function

  function [par,data] = eqC13(x0,par,data,Cpool)
    global GC13;
    on   = true  ;
    off  = false ;
    iwet = par.iwet; 
    nwet = par.nwet;
    par.nC = length(Cpool);
    % initialize C13 model
    GC13 = zeros(par.nC*length(iwet),1);

    par.debug13 = off;
    [par, C13 ] = eqC13cycle_v2(x0, par);
    
    % distribute C13 solution to par and ata
    [par,data] = dispSol(C13,par,data,Cpool);
  end

  function [par,data] = eqC14(x0,par,data,Cpool)
    global GC14;
    on   = true  ;
    off  = false ;
    iwet = par.iwet; 
    nwet = par.nwet;
    par.nC = length(Cpool);
    % initialize C14 model
    GC14 = zeros(par.nC*length(iwet),1);

    par.debug14 = off;
    [par, C14 ] = eqC14cycle_v2(x0, par);
    
    % distribute C14 solution to par and ata
    [par,data] = dispSol(C14,par,data,Cpool);
  end

  function [par,data] = dispSol(C,par,data,Cpool)
  % dispense the model solution to par.Variable
  % and data.Variables; Cpool is a cell containing
  % the tracers used in the individule model
    M3d = par.M3d;
    iwet = par.iwet;
    nwet = par.nwet;
    for i = 1:length(Cpool)
        eval(sprintf('%s=M3d+nan;',Cpool{i}));
        eval(sprintf('%s(iwet)=C(%i*nwet+1:%i*nwet);',Cpool{i},i-1,i)); 
        eval(sprintf('par.%s=%s(iwet);',Cpool{i},Cpool{i}));
        eval(sprintf('data.%s=%s;',Cpool{i},Cpool{i}));
    end
  end

  function par = getname(par)
    % generate file name with different choices of model combination
    % the prefix in the form of 'PC','PCO','PCOSi', etc
    Models = {'P','C','O','Si'};
    s2e = '';
    for i = 1:length(Models)
      if eval(sprintf('par.%smodel',Models{i}))
        s2e = strcat(s2e,Models{i});
      end
    end
  end

  end % end method
end % end classdef
