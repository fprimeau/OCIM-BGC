function [Cratio] = get_C14_jtr(varargin)
% function [Cratio] = get_atmC_ratio(varargin)
% varargin,i.e., tt = 1850:1/12:2015;  
% some parameters to convert from \deltaC to C13/C or C14/C ratio
  R13pdb = 1.12372*1e-2;
  R14oxa = 1.176*1e-12; 
  d2r13 = @(x) (1+x*1e-3)*R13pdb;
  D2r14 = @(x,d13) (1+x*1e-3)./((0.975./(1+d13*1e-3)).^2)*R14oxa;
  
  fprintf('read Jim C14 data...');
  % fname = '/export/nfs0home/weiweif/Research/CYCLOCIM/DATA/D14C/Graven2017/newmergedtimeseries.txt';
  fname = '/DFS-L/DATA/primeau/oceandata/DATA/D14C/Graven2017/newmergedtimeseries.txt';
  atmdata = csvread(fname,1,0);
  % disp(sprintf('[%5.1f %5.2f %5.3f %5.2f %5.3f];      ...\n',atmdata));

  fprintf('get ratios from \x03b4 value...\n');
  % R13a = d2r13(atmdata(:,5));
  R13a = d2r13(-6.6);
  % R14aN = D2r14(atmdata(:,2),R13a);
  % R14aT = D2r14(atmdata(:,2),R13a);
  % R14aS = D2r14(atmdata(:,4),R13a);
  R14a = D2r14(atmdata(:,2),R13a);

  Cratio.R13a = R13a;
  % Cratio.R14aN = R14aN;
  % Cratio.R14aT = R14aT;
  % Cratio.R14aS = R14aS;
  % Cratio.R14a = (R14aN + R14aT + R14aS)/3;
  Cratio.R14a  = R14a;
  Cratio.Ryear = atmdata(:,1);
  Cratio.d14a  = atmdata(:,2);
  % Cratio.d13a = atmdata(:,5);
  % Cratio.d14aN = atmdata(:,2);
  % Cratio.d14aT = atmdata(:,2);
  % Cratio.d14aS = atmdata(:,4);
  % Cratio.d14a = (atmdata(:,2) + atmdata(:,3) + atmdata(:,4))/3;
  Cratio.d14a = atmdata(:,2);
  if nargin > 0
    tt = varargin{1};
    disp('Interpolate data to time steps');
    % Cratio.pc13R_nstep = interp1(Cratio.Ryear,Cratio.R13a,tt);
    Cratio.pc14R_nstep = interp1(Cratio.Ryear,Cratio.R14a,tt);
    % Cratio.pc13R_nstep = interp1(Cratio.Ryear,Cratio.R13a,tt);
    Cratio.pc13R_nstep = Cratio.pc14R_nstep*0 + Cratio.R13a;
    Cratio.tt = tt;
    Cratio.pd14a_nstep = interp1(Cratio.Ryear,Cratio.d14a,tt);
  end
end
