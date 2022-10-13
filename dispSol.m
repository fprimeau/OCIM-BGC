function [par,data] = dispSol(C,par,data,Cpool)
% function [par,data] = dispSol(C,par,data,Cpool)
% dispense the model solution to par.Variable
% and data.Variables; Cpool is a cell containing
% the tracers used in the individule model
M3d = par.M3d;
nwet = par.nwet;
iwet = par.iwet;
  for i = 1:length(Cpool)
      eval(sprintf('%s=M3d+nan;',Cpool{i}));
      eval(sprintf('%s(iwet)=C(%i*nwet+1:%i*nwet);',Cpool{i},i-1,i));
      eval(sprintf('par.%s=%s(iwet);',Cpool{i},Cpool{i}));
      eval(sprintf('data.%s=%s;',Cpool{i},Cpool{i}));
  end
end
