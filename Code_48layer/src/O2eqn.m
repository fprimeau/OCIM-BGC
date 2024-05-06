function [f,J,par] = O2eqn(X, par)    
    on = true; off = false;
    % fixed parameters
    iwet  = par.iwet  ;
    nwet  = par.nwet  ;
    TRdiv = par.TRdiv ;
    I     = speye(nwet) ;
    PO4   = par.po4obs(iwet) ;

    %
    O2 = X(0*nwet+1:1*nwet) ;
    
    % variables from C model
    %
    if (par.Cisotope == on) ;
        POC  = par.POC  + par.POC13   ;
        DOC  = par.DOC  + par.DOC13   ;
        DOCl = par.DOCl + par.DOC13l  ;
        DOCr = par.DOCr + par.DOC13r  ;
    else
        POC  = par.POC    ;
        DOC  = par.DOC    ;
        DOCl = par.DOCl   ;
        DOCr = par.DOCr   ;
    end
    Tz   = par.Tz     ;
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
    %
    % tunable parameters;
    O2C_T = par.O2C_T ; 
    rO2C  = par.rO2C  ;
    kappa_l = par.kappa_l ;
    kappa_p = par.kappa_p ;
    tf    = (par.vT - 30)/10 ;
    kC    = d0( kdC * Q10C .^ tf ) ;
    %
    % O2 saturation concentration
    vout  = Fsea2air(par,'O2') ;
    KO2   = vout.KO2   ;
    o2sat = vout.o2sat ;
    % rate of o2 production
    O2C = O2C_T*Tz + rO2C ; 
    C2P = par.C2P         ;
    G   = par.G           ;
    PO2 = d0(par.Cnpp(iwet))*O2C  ;

    % parobolic function for o2 consumption
    R      = 0.5 + 0.5*tanh(O2-5)    ;
    dRdO   = 0.5 - 0.5*tanh(O2-5).^2 ;
    d2RdO2 = -2*d0(dRdO)*tanh(O2-5)  ;

    % rate of o2 utilization
    kappa_r = kru*UM + krd*DM ;
    eta     = etau*WM ;
    % eta    = etau*UM + etad*DM ;
    LO2    = d0(eta*kC*DOC + kappa_r*DOCr + kappa_l*DOCl + kappa_p*POC)*O2C.*R  ;
    dLdO   = d0(eta*kC*DOC + kappa_r*DOCr + kappa_l*DOCl + kappa_p*POC)*O2C.*dRdO ;
    d2LdO2 = d0(eta*kC*DOC + kappa_r*DOCr + kappa_l*DOCl + kappa_p*POC)*O2C.*d2RdO2  ;
   
    % O2 function
    F = TRdiv*O2 - PO2 + LO2 - KO2*(o2sat-O2) ;

    %F = J*X + f ;
    % Jacobian
    J = TRdiv + KO2 ;

    % source and sink term
    f = -PO2 + LO2 - KO2*(o2sat) ;
end


