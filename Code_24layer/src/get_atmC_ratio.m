function [Cratio] = get_atmC_ratio(varargin)
% function [Cratio] = get_atmC_ratio(varargin)
% varargin,i.e., tt = 1850:1/12:2015;  
% some parameters to convert from \deltaC to C13/C or C14/C ratio
  R13pdb = 1.12372*1e-2;
  R14oxa = 1.176*1e-12; 
  d2r13 = @(x) (1+x*1e-3)*R13pdb;
  D2r14 = @(x,d13) (1+x*1e-3)./((0.975./(1+d13*1e-3)).^2)*R14oxa;
  
  %
  fprintf('read Graven 2017 data...');
  fname = '/DFS-L/DATA/primeau/oceandata/DATA/D14C/Graven2017/TableS1_wfu.csv';
  atmdata = csvread(fname,5,0);
  
  %
  fprintf('get ratios from \x03b4 value...\n');
  R13a  = d2r13(atmdata(:,5)); % 
  R14aN = D2r14(atmdata(:,2),R13a); % North Hemisphere
  R14aT = D2r14(atmdata(:,3),R13a); % Tropic
  R14aS = D2r14(atmdata(:,4),R13a); % South Hemisphere

  Cratio.R13a  = R13a;
  Cratio.R14aN = R14aN;
  Cratio.R14aT = R14aT;
  Cratio.R14aS = R14aS;
  Cratio.R14a  = (R14aN + R14aT + R14aS)/3;
  Cratio.Ryear = atmdata(:,1);
  Cratio.d13a  = atmdata(:,5);
  Cratio.d14aN = atmdata(:,2);
  Cratio.d14aT = atmdata(:,3);
  Cratio.d14aS = atmdata(:,4);
  Cratio.d14a  = (atmdata(:,2) + atmdata(:,3) + atmdata(:,4))/3;
  if nargin > 0
    tt = varargin{1};
    disp('Interpolate data to time steps');
    pc13R_nstep = interp1(Cratio.Ryear-0.5,Cratio.R13a,tt);
    pc13R = permute(repmat(pc13R_nstep,[91 1 180]),[1 3 2]);
    Cratio.pc13R_nstep = pc13R_nstep;
    Cratio.pc13R = pc13R;
    pc14R_nstep = interp1(Cratio.Ryear-0.5,Cratio.R14a,tt);
    Cratio.pc14R_nstep = pc14R_nstep;
    pc14RN_nstep = interp1(Cratio.Ryear-0.5,Cratio.R14aN,tt);
    pc14RN = permute(repmat(pc14RN_nstep,[91 1 180]),[1 3 2]);
    pc14RT_nstep = interp1(Cratio.Ryear-0.5,Cratio.R14aT,tt);
    pc14RT = permute(repmat(pc14RT_nstep,[91 1 180]),[1 3 2]);
    pc14RS_nstep = interp1(Cratio.Ryear-0.5,Cratio.R14aS,tt);
    pc14RS = permute(repmat(pc14RS_nstep,[91 1 180]),[1 3 2]);
    pc14R = pc14RN*nan;
    pc14R(1:33,:,:) = pc14RS(1:33,:,:);
    pc14R(34:59,:,:) = pc14RT(34:59,:,:);
    pc14R(60:91,:,:) = pc14RN(60:91,:,:);
    Cratio.pc14R = pc14R;
  end
end
