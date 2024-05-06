function [f,J,par] = C12eqn(X, par)    
% unpack some useful stuff
    on = true; off = false;
    grd   = par.grd   ;
    M3d   = par.M3d   ;
    TRdiv = par.TRdiv ;
    iwet  = par.iwet  ;
    nwet  = par.nwet  ;
    dVt   = par.dVt   ;
    I     = par.I     ;
    
    Tz  = par.Tz ;
    DIC  = X(0*nwet+1:1*nwet) ; 
    POC  = X(1*nwet+1:2*nwet) ;
    DOC  = X(2*nwet+1:3*nwet) ;
    PIC  = X(3*nwet+1:4*nwet) ;
    ALK  = X(4*nwet+1:5*nwet) ;
    DOCl = X(5*nwet+1:6*nwet) ;
    DOCr = X(6*nwet+1:7*nwet) ;

    PO4 = par.po4obs(iwet) ;
    % fixed parameters
    kappa_p = par.kappa_p;
    kappa_l = par.kappa_l;
    kPIC    = par.kPIC   ;
    gamma = par.gamma    ;
    % parameters need to be optimized
    alpha = par.alpha    ;
    beta  = par.beta     ;
    sigP  = par.sigP     ;
    sigC  = par.sigC     ;
    kru   = par.kru      ;
    krd   = par.krd      ;
    etau  = par.etau     ;
    etad  = par.etad     ; 
    bC_T  = par.bC_T     ;
    bC    = par.bC       ;
    d     = par.d        ;
    Q10C  = par.Q10C     ;
    kdC   = par.kdC      ;
    R_Si  = par.R_Si     ;
    rR    = par.rR       ;
    cc    = par.cc       ;
    dd    = par.dd       ;
    pme   = par.pme      ;
    % PIC to POC rain ratio 
    vout  = mkPIC2P(par) ;
    RR    = vout.RR      ;
    RR_Si = vout.RR_Si   ;
    RR_rR = vout.RR_rR   ;
    clear vout 
    % kappa_dc ;
    tf    = (par.vT - 30)/10 ;
    kC    = d0( kdC * Q10C .^ tf ) ;
    C2P   = 1./(cc*PO4 + dd) ;
    N2C   = 18/106 ; %16/117 ; 
    par.C2P = C2P  ;

    % particle flux div_rergence [s^-1];
    PFDa = buildPFD_48layer(par, 'PIC') ;
    PFDc = buildPFD_48layer(par, 'POC') ;
    par.PFDa = PFDa ;
    par.PFDc = PFDc ;
    par.DIC  = DIC  ;
    par.ALK  = ALK  ;

    % Air-Sea gas exchange
    vout  = Fsea2air(par, 'CO2');
    G_dic = vout.G_dic ;
    G_alk = vout.G_alk ;
    JgDIC = vout.JgDIC ;
    clear vout 

    % biological DIC uptake operator
    G = uptake_C(par)  ; par.G = G ;
    
    kappa_g = par.kappa_g ;
    ALKbar  = par.ALKbar  ;
    sDICbar = par.sDICbar ;
    sALKbar = par.sALKbar ;
    
    UM = par.UM ; 
    DM = par.DM ;
    WM = par.WM ;
    kappa_r =  kru*UM +  krd*DM ;
    eta     = etau*WM ;
    % eta     = etau*UM + etad*DM ;
    
    dDICdt  = TRdiv*DIC + (1-sigC-gamma)*RR*G*C2P - eta*(kC*DOC) ...
              - kPIC*PIC - JgDIC + pme*sDICbar - kappa_r*DOCr ...
              - kappa_l*DOCl - kappa_p*POC + par.Cnpp(iwet) ;  %FDIC
    
    dPOCdt  = (PFDc+kappa_p*I)*POC - (1-sigC-gamma)*G*C2P   ; % FPOC
    
    dDOCdt  = (TRdiv+kC)*DOC - sigC*G*C2P  ; % FDOC
    
    dPICdt  = (PFDa + kPIC*I)*PIC - (1-sigC-gamma)*RR*G*C2P ; % FPIC
    
    dALKdt  = TRdiv*ALK + 2*(1-sigC-gamma)*RR*G*C2P - 2*kPIC*PIC ...
              - N2C*par.Cnpp(iwet) + N2C*(eta*(kC*DOC) + kappa_r*DOCr + ...
                                 kappa_l*DOCl + kappa_p*POC) ...
              + pme*sALKbar + kappa_g*(ALK - ALKbar) ;  % ALK 
    
    dDOCldt = (TRdiv+kappa_l*I)*DOCl - (par.Cnpp(iwet) - G*C2P) ;  % DOCl 
    
    dDOCrdt = (TRdiv+kappa_r)*DOCr - (I-eta)*(kC*DOC) ; % DOCr 
    
    zro = zeros(nwet,1);
    F   = [dDICdt; dPOCdt; dDOCdt; dPICdt; dALKdt; dDOCldt; dDOCrdt];
    f   = [(1-sigC-gamma)*RR*G*C2P - JgDIC + pme*sDICbar + par.Cnpp(iwet) ;...
           - (1-sigC-gamma)*G*C2P; ...
           - sigC*G*C2P; ...
           - (1-sigC-gamma)*RR*G*C2P; ...
           2*(1-sigC-gamma)*RR*G*C2P - N2C*par.Cnpp(iwet) ...
                 + pme*sALKbar - kappa_g*ALKbar;...
           - (par.Cnpp(iwet) - G*C2P); ...
           zro];
    

    % construct the LHS matrix for the offline model
    % disp('Preparing LHS and RHS matrix:')
    % colum 1 dFdDIC
    Jc{1,1} = TRdiv - G_dic ; 
    Jc{2,1} = 0*I ;
    Jc{3,1} = 0*I ;
    Jc{4,1} = 0*I ;
    Jc{5,1} = 0*I ;
    Jc{6,1} = 0*I ;
    Jc{7,1} = 0*I ;
    % colum 2 dFdPOC
    Jc{1,2} = -kappa_p*I ;
    Jc{2,2} = PFDc + kappa_p*I ;
    Jc{3,2} = 0*I ;
    Jc{4,2} = 0*I ;
    Jc{5,2} = N2C*kappa_p*I ;
    Jc{6,2} = 0*I ;
    Jc{7,2} = 0*I ;
    % colum 3 dFdDOC
    Jc{1,3} = -eta*kC ;
    Jc{2,3} = 0*I ;
    Jc{3,3} = TRdiv + kC ;
    Jc{4,3} = 0*I ;
    Jc{5,3} = eta*N2C*kC ;
    Jc{6,3} = 0*I ;
    Jc{7,3} = -(I-eta)*kC ;
    % colum 4 dFdPIC
    Jc{1,4} = -kPIC*I ;
    Jc{2,4} = 0*I ;
    Jc{3,4} = 0*I ;
    Jc{4,4} = PFDa + kPIC*I ;
    Jc{5,4} = -2*kPIC*I ;
    Jc{6,4} = 0*I ;
    Jc{7,4} = 0*I ;
    % column 5 dFdALK
    Jc{1,5} = -G_alk;
    Jc{2,5} = 0*I ;
    Jc{3,5} = 0*I ;
    Jc{4,5} = 0*I ;
    Jc{5,5} = TRdiv + kappa_g*I ;
    Jc{6,5} = 0*I ;
    Jc{7,5} = 0*I ;
    % column 6 dFdDOCl
    Jc{1,6} = -kappa_l*I ;
    Jc{2,6} = 0*I ;
    Jc{3,6} = 0*I ;
    Jc{4,6} = 0*I ;
    Jc{5,6} = N2C*kappa_l*I ;
    Jc{6,6} = TRdiv + kappa_l*I ;
    Jc{7,6} = 0*I ;
    % column 7 dFdDOCr
    Jc{1,7} = -kappa_r ;
    Jc{2,7} = 0*I ;
    Jc{3,7} = 0*I ;
    Jc{4,7} = 0*I ;
    Jc{5,7} = N2C*kappa_r ;
    Jc{6,7} = 0*I ;
    Jc{7,7} = TRdiv + kappa_r ;
    % Jacobian matrix
    J = cell2mat(Jc);
    % F = J*X + f;
end

