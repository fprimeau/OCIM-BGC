% clear all; close all;

disp(sprintf('loading steady state from %s',par.output_dir));
load([par.output_dir par.output_eq],'par');
load /DFS-L/DATA/primeau/oceandata/DATA/eSST_obs_1749_2023.mat;
par.Temp_eSST = Temp_eSST;
par.Temp_eSST_tt = Temp_eSST_tt;
par.Temp_obs  = par.Temp; % keep a copy of WOA temp in case it is updated by temperature model
par.saveall = false;

Ctype = {'C12','C13','C14'}; % Ctype = {'C12','O2'};
for m = 1:length(Ctype)
  vn = Ctype{m};
  if contains(vn,'12')
    Xin.C12 = [par.DIC;par.POC;par.DOC;par.PIC;par.ALK;par.DOCl;par.DOCr];
  elseif contains(vn,'O2')
    Xin.O2 = par.O2;
  elseif contains(vn,'TT')
    Xin.TT = par.TT; % steay state by solving TTeqn
  elseif contains(vn,'P')
    Xin.P = [par.DIP;par.POP;par.DOP;par.DOPl];
  else
    eval(sprintf('Xin.%s = [par.DI%s;par.PO%s;par.DO%s;par.PI%s;par.DO%sl;par.DO%sr];'...
                        ,vn,vn,vn,vn,vn,vn,vn));
  end
end
  
%t0 = 1750; t1 = 2022;
t0=2022; t1=2023;
[Xout,Tout] = time_stepper(par,t0,t1,Xin,Ctype);
