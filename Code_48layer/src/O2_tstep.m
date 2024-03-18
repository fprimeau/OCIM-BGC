function [O2tstep, O2integral, O2mean, par] = O2_tstep(X0, par)
% compute the initial O2 field using timestepping method
% The reult would be used in initial O2 field in Newton's method (nsnew)
% output O2tstep: 4-d array: 91 x 180 x 48 x time
        
    % fixed parameters
    iwet  = par.iwet  ;
    nwet  = par.nwet  ;
    TRdiv = par.TRdiv ;
    I     = speye(nwet) ;
    PO4   = par.po4obs(iwet) ;
    dVt   = par.dVt   ;
    spa   = par.spa   ;

    % variables from C model
    POC   = par.POC   ;
    DOC   = par.DOC   ;
    DOCl  = par.DOCl  ;
    DOCr  = par.DOCr  ;
    Tz    = par.Tz    ;
    Q10C  = par.Q10C  ;
    kdC   = par.kdC   ;
    kru   = par.kru   ;
    krd   = par.krd   ;
    etau  = par.etau  ;
    etad  = par.etad  ;
    UM    = par.UM    ;
    DM    = par.DM    ;
    WM    = par.WM    ;
    cc    = par.cc    ;
    dd    = par.dd    ;
    
    % tunable parameters;
    O2C_T = par.O2C_T ; 
    rO2C  = par.rO2C  ;
    kappa_l = par.kappa_l ;
    kappa_p = par.kappa_p ;
    tf    = (par.vT - 30)/10 ;
    kC    = d0( kdC * Q10C .^ tf ) ;
    
    % O2 saturation concentration
    vout  = Fsea2air(par,'O2') ;
    KO2   = vout.KO2   ;
    o2sat = vout.o2sat ;

    % rate of o2 production
    O2C = O2C_T*Tz + rO2C ; 
    C2P = 1./(cc*PO4 + dd) ;
    G   = par.G           ;
    PO2 = d0(par.Cnpp(iwet))*O2C  ;

    % parobolic function for O2(X0) consumption
    R      = 0.5 + 0.5*tanh(X0-5)    ;
    
    % rate of o2 utilization
    kappa_r = kru*UM + krd*DM ;
    eta     = etau*WM ;
    
    % eta    = etau*UM + etad*DM ;
    LO2    = d0(eta*kC*DOC + kappa_r*DOCr + kappa_l*DOCl + kappa_p*POC)*O2C.*R  ;
    
    % Using a semi-implicit scheme for time stepping method
    % Trapezoide rule for all terms, except for Respiration (Euler Forward)
    %
    %  d
    % -- O2 +TRdiv*O2(t) = PO2 - LO2 .*R(t) + KO2*(O2sat - O2(t))  
    % dt
    %

    % time-step size schedule
    %                  1 hr  1 day  1 month 1 year 4 years   10 years
    dt_size = spa./[ 365*24   365       12     1     0.25       0.1 ];  
    nsteps  =      [  100     100      100   200     200        300 ];  

    %
    O2integral = 0;
    O2tstep    = zeros(length(X0), length(nsteps));
    j = 0; 
    t = 0;
    


    for k = 1:length(dt_size) % loop over step size schedule
        % set the step size
        dt = dt_size(k);
        % set the number of time steps
        n = nsteps(k);
        % trapezoid rule
        A = I + (dt/2)*(TRdiv+KO2);
        B = I - (dt/2)*(TRdiv+KO2);
        fprintf('factoring the big matrix...'); tic
        FA = mfactor(A);
       
        fprintf('dt = %f nsteps = %i \n', dt, n);
        for i = 1:n   % loop over time steps
            X0 = mfactor(FA, (B*X0 + dt*(-LO2.*R + PO2 +KO2*o2sat)));
            t = t+dt;
            j = j+1;
            O2tstep(:,j)  = X0;
            O2integral(j) = sum(dVt(iwet).*X0);
            O2mean(j)     = sum(dVt(iwet).*X0)/sum(dVt(iwet));
            T(j) = t;
            fprintf('.');
        end
        fprintf('\n');
    end
    
    Tyear = T/spa;
end




