function [Xout,Tout,par] = time_stepper(par,t0,t1,Xin,Ctype)
  % integrate from t0 to t1 where t0 and t1 are in units of calendar year
  
  % Unpack the useful parameters
  spa = par.spa
  dt = spa/12; %monthly time step
  tt = t0:dt/spa:t1;
  nstep = length(tt);
  Nstep = spa/dt;
  nwet = par.nwet;
  M2d  = par.M3d(:,:,1);
  isrf = find(M2d);
  disp(sprintf('run time: %i--->%i, time steps: %4.1f',t0,t1,nstep)); 

  % load the atmospheric CO2 time history
  atm = importdata('atmhist.txt');
  pco2atm = atm.data(:,5); % 그냥 atm이 아니고...?
  pco2year = atm.data(:,1);
  % interp pco2atm to match model temporal resolution
  pco2atm_nstep = interp1(pco2year,pco2atm,tt);

  % obtain C13 and C14 ratio from Graven et al. (2017)
  [Cratio] = get_atmC_ratio(tt);
  pc13R_nstep = Cratio.pc13R_nstep;
  pc14R_nstep = Cratio.pc14R_nstep; % for global
  % pc14R = Cratio.pc14R; % for North, Tropical and South
  % obtain C14 from the data provided by Jim
  % Jratio = get_C14_jtr(tt);
  % pc13R_nstep = Jratio.pc13R_nstep;
  % pc14R_nstep = Jratio.pc14R_nstep;
  % pc14R_nstep = interp1(Cratio.Ryear-0.5,Cratio.R14a,tt);
  % keyboard;
  %
  % for j=1:size(par.Temp_eSST,2)
  %  for i=1:size(par.Temp_eSST,1)
      % TT_nstep(i,j,:) = interp1((1850+1/12:1/12:2015)-0.5,squeeze(par.Temp_eSST(i,j,:)),tt);
   %   TT_nstep(i,j,:) = interp1(par.Temp_eSST_tt,squeeze(par.Temp_eSST(i,j,:)),tt);
  %  end
  % end
  % par.pco2atm = 285.91;% function of t in general
  % par.c13.R13a = 0.01116303448; % in 1850
  % par.c14.R14a = 1.220805*1e-12; % will be changed later, but needs to be initialized
  % par.SST_obs = TT_nstep(:,:,1); 

  % Ctype = {14};
  % set up A,B,FA for time-stepping C12,C13,C14 and O2
  for m = 1:length(Ctype)
    vn = Ctype{m};
    astr = sprintf('%s = Xin.%s;',vn,vn); % initial conditions 
    fprintf('...%s',astr);
    eval(astr);
    if par.saveall
      eval(sprintf('Xout.%s = zeros(length(%s),nstep);',vn,vn));
      eval(sprintf('Xout.%s(:,1) = %s;',vn,vn));
    else
      eval(sprintf('Xout.%s = zeros(length(%s),Nstep);',vn,vn));
      eval(sprintf('Xout.%s(:,1) = %s;',vn,vn));
    end
    % keyboard;
    bstr = sprintf('[A%s,B%s,FA%s,par]=setTimeStepper(%s,par,dt,"%s");',vn,vn,vn,vn,vn);
    fprintf('...%s \n',bstr);
    eval(bstr); 
  end
  % [A12,B12,FA12,Xout12,Tout12,par] = setTimeStepper(Xin12,par,dt,tt,'C12');

  Tout = zeros(1,nstep);

  t = t0;
  Tout(1) = t;
  for i = 2:nstep
    par.pco2atm = pco2atm_nstep(i);
    par.c13.R13a = pc13R_nstep(i);
    % C14 for North, Tropical and South
    % tmp = M2d*0; tmp = pc14R(:,:,i);
    % par.c14.R14a = tmp(isrf);
    % C14 uniform for global
    par.c14.R14a = pc14R_nstep(i);
    % changing SST with time
    % par.Temp(:,:,1) = TT_nstep(:,:,i);
    % par.SST_obs = TT_nstep(:,:,i);
    % Caution: the above variables are updated every time step

    fprintf('...Time: %4.2f, AtmCO2: %4.2f,R13a: %7.6g, R14a: %7.6g, %s: %7.6g\n',...
                t,max(par.pco2atm(:)), max(par.c13.R13a),max(par.c14.R14a(:)),...
                vn,eval(sprintf('mean(%s(1:nwet))',vn)));

    % time stepping C12, C13,C14, and O2
    for m = 1:length(Ctype)
      vn = Ctype{m};
      % call C12eqn(X12,par) or C13eqn(X13,par) or
      cstr = sprintf('[f%s,J%s,par]=%seqn(%s,par);',vn,vn,vn,vn);
      eval(cstr);
      % time stepping with mfactor
      dstr = sprintf('%s=mfactor(FA%s,B%s*%s-dt*f%s);',vn,vn,vn,vn,vn);
      eval(dstr);
      % assign solution X12(X13,X14) to XC12 (XC13,XC14)
      if par.saveall
        eval(sprintf('Xout.%s(:,i) = %s;',vn,vn));
      else
        ix = mod(i,Nstep);
        if ix == 0
          ix = Nstep;
        end
        eval(sprintf('Xout.%s(:,ix) = %s;',vn,vn));
        if isfield(par,'JgDIC') & isfield(par,'co2surf')
          Xout.JgDIC(:,ix) = par.JgDIC; Xout.co2surf(:,ix) = par.co2surf; Xout.pco2surf(:,ix) = par.pco2;
        end
      end
      %keyboard
      if contains(vn,'O2')
        bstr = sprintf('[A%s,B%s,FA%s,par]=setTimeStepper(%s,par,dt,"%s");',vn,vn,vn,vn,vn);
        % fprintf('...%s is updated \n',bstr);
        eval(bstr); 
      end %------------> 여기 필요 없는거 아닌가. 아니다 필요할듯. C cycle 변하면 parameter가 변하니까.
      % if contains(vn,'TT')
      %  par.Temp(par.iwet) = TT; % update 3D ocean temperature
      % end
    end

    
    if ~par.saveall
      if mod(i,Nstep) == 0
       % outname = sprintf('%smodel_%s_%i.mat',par.output_dir,strjoin(Ctype,'_'),i/Nstep);
         outname = sprintf('Cisotope_model_48layer_%s_%i.mat',strjoin(Ctype,'_'),240426);
        save(outname,'Xout', 'Tout', '-v7.3');
        sprintf('saving file')
      end
    end
    % [f12,J12,par] = C12eqn(X12,par);
    % X12 = mfactor(FB12,A12*X12+dt*f12);
    % X12 = mfactor(FA12,B12*X12-dt*f12);
    % time stepping C13
    % Caution: DIC/ALK is passed in from par.DIC/ALK, e.g., Fsea2air flux
    % time stepping C14

    % update time t
    t = t+dt/spa;
    % save the solution
    % Xout12(:,i) = X12;
    Tout(i) = t;

  end % nstep
  Xout.pc13R_nstep = pc13R_nstep;
  Xout.pc14R_nstep = pc14R_nstep;
  
end % time stepper

function [A,B,FA,par] = setTimeStepper(X,par,dt,Ctype)
  % get variable source-sink terms and Jacobian 
  % Ctype can be 12,13, or 14
  % [f,J,par] = C12eqn(Xin,par);
  fprintf('...Set up A, B, FA for %s ...',Ctype); 
  [f,J,par] = eval(sprintf('%seqn(X,par);',Ctype));
  I = speye(length(X));
  % trapezoid + euler forward
  % dX/dt + F(X,t) =0 where F(X,t) = J*X + f(X,t)
  % see dDICdt.... in C12eqn.m
  % X(n+1) - X(n) + (dt/2)*J*(X(n+1)+X(n)) + dt*fn = 0
  % ==> A*X(n+1) = B*X(n) - dt*f(n)
  % X(n+1) = A\(B*X(n)-dt*f(n))
  A = I + 0.5*dt*J;
  B = I - 0.5*dt*J;
  tic
  fprintf('Factor the FA of %s...\n',Ctype);
  % FB = mfactor(B);
  FA = mfactor(A);
  toc
  % nstep = length(tt);
  % Xout = zeros(length(X),nstep);
  % Tout = zeros(1,nstep);
  % Xout(:,1) = X;
  % Tout(1) = tt(1);
end
 
