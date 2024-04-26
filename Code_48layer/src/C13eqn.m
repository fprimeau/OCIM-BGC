function [f,J,par] = C13eqn(X, par)    
% unpack some useful stuff
    on = true; off = false;
    grd   = par.grd   ;
    M3d   = par.M3d   ;
    TRdiv = par.TRdiv ;
    iwet  = par.iwet  ;
    nwet  = par.nwet  ;
    dVt   = par.dVt   ;
    I     = par.I     ;
    % get first layer index
    tmp  = M3d      ;
    tmp(:,:,2:end) = 0     ;
    isrf = find(tmp(iwet)) ;
    
    Tz  = par.Tz ;
    DIC13  = X(0*nwet+1:1*nwet) ; 
    POC13  = X(1*nwet+1:2*nwet) ;
    DOC13  = X(2*nwet+1:3*nwet) ;
    PIC13  = X(3*nwet+1:4*nwet) ;
    % ALK  = X(4*nwet+1:5*nwet) ;
    DOC13l = X(4*nwet+1:5*nwet) ;
    DOC13r = X(5*nwet+1:6*nwet) ;

    R13o = DIC13./(par.DIC); 
    dR13o = d0(1./par.DIC);

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
    PFDa = buildPFD(par, 'PIC') ;
    PFDc = buildPFD(par, 'POC') ;
    par.PFDa = PFDa ;
    par.PFDc = PFDc ;
    par.DIC13  = DIC13  ;
    %par.ALK  = ALK  ;
    
    % define C13 fractionation factors and ratios in ocean and atmosphere
    % set up fractionation factors fro C13 and R13
    %  A. Schmittner et al.: Distribution of carbon isotope ratios (δ13C) in the ocean
    % par.c13.R13a = 0.01116303448;%0.01118; % air C13/C 0.01116303448 at 1850
    par.pc13atm = par.pco2atm*par.c13.R13a;
    par.c13.R13o = R13o;
    par.c13.dR13o = dR13o;
    par.c13.alpha_k = 0.99915; % kenetic fractionation factor 
    par.c13.alpha_g2aq = 0.998764; % gas to water fractionation factor
    % temperature (in C)-dependent equilibrium fractionation factor from gaseous CO2 to DIC.
    % par.c13.alpha_g2dic = 1.01051-1.05*1e-4*par.Temp(isrf); 
    disp(sprintf('air-sea fractionation tuned by %4.2f',par.fras));
    par.c13.alpha_g2dic = (1.01051-1.05*1e-4*par.Temp-1.0)*par.fras + 1.0;
    % par.c13.alpha_g2dic = 1.01051-1.05*1e-4*par.Temp; 
    
    % Air-Sea gas exchange for total C
    vout    = Fsea2air(par, 'CO2');
    G_dic   = vout.G_dic ;
    G_alk   = vout.G_alk ;
    JgDIC   = vout.JgDIC ;
    co2surf = vout.co2surf;
    clear vout;
    
    % the equilibrium fractionation factor from aqueous CO2 to particulate organic carbon (POC) 
    co2 = M3d; nz = size(M3d,3);
    % co2(:,:,1) = co2surf; co2 = co2(:,:,ones(nz,1)); 
    co2((iwet(isrf))) = co2surf; co2 = co2(:,:,ones(nz,1)); 
    % WARNING: we should resolve the CO2 system at  all the layers  where we have
    % biological production. For now we approximate the CO2 at all layers
    % using co2surf. 
    
    alpha_aq2poc = -0.017*log(co2) + 1.0034; % check the unit of co2surf
    alpha_tmp = ( par.c13.alpha_g2aq ./ par.c13.alpha_g2dic ) .* alpha_aq2poc; 
    alpha_dic2poc = (alpha_tmp(iwet)-1.0)*par.frpho + 1.0;
    par.c13.alpha_aq2poc  = alpha_aq2poc   ;
    par.c13.alpha_tmp     = alpha_tmp      ;
    par.c13.alpha_dic2poc = alpha_dic2poc  ;

    if par.debug13
      disp('Testing fractionation factors')
      par.c13.alpha_k = 1;%0.99915; % kenetic fractionation factor 
      par.c13.alpha_g2aq = 1;%0.998764; % gas to water fractionation factor
      par.c13.alpha_g2dic = 1.01051*0 - 1.05*1e-4*par.Temp*0 + 1;
      alpha_dic2poc = alpha_dic2poc*0 + 1.0; % debug
    end
	
    % Air-Sea gas exchange for C13
    vout  = Fsea2air(par, 'C13');
    JgDIC13 = vout.JgDIC13;
    G_dic13 = vout.G_dic13;
    rhs13 = vout.rhs13;
    % biological DIC uptake operator
    G = uptake_C(par)  ; par.G = G ;
     
    kappa_g = par.kappa_g ;
    ALKbar  = par.ALKbar  ;
    % sDICbar = par.sDICbar ;
    sALKbar = par.sALKbar ;

    % recalculate sDICbar because it is now changing with time
    sDICbar = sum( par.DIC(isrf) .* dVt(iwet(isrf)) ) / sum( dVt(iwet(isrf)) );

    UM = par.UM ; 
    DM = par.DM ;
    WM = par.WM ;
    kappa_r =  kru*UM +  krd*DM ;
    eta     =  etau*WM ;
    % eta     = etau*UM + etad*DM ;

    % eq1 = TRdiv*DIC13 ...                         % advective-diffusive transport
    %       + G*d0(C2P.*alpha_dic2poc)*R13o ...     % removal of dic13 organic c13 production
    %       + (1-sigC-gamma)*RR*G*d0(C2P)*R13o  ... % removal of dic13 due to pic13 production
    %       - kPIC*PIC13 ....          % dissolution of PIC13
    %       - JgDIC13 ...              % air-sea gas exchange
    %       + sDICbar*d0(pme)*R13o ... % concentration and dillution due to precip and evaporation
    %       - eta*(kC*DOC13) ...       % respiration of DOC13
    %       - kappa_r*DOC13r ...       % respiration of DOC13r
    %       - kappa_l*DOC13l...        % respiration of DOC13l
    %       - kappa_p*POC13 ;          % respiration of PCO13
    %
    % eq2 = (PFDc+kappa_p*I)*POC13 - (1-sigC-gamma)*G*d0(C2P.*alpha_dic2poc)*R13o;   ; % FPOC
    %
    % eq3 = (TRdiv+kC)*DOC13 - sigC*G*d0(C2P.*alpha_dic2poc)*R13o  ; % FDOC
    %
    % eq4 = (PFDa + kPIC*I)*PIC13 - (1-sigC-gamma)*RR*G*d0(C2P)*R13o ; % FPIC
    %
    % eq5 = (TRdiv+kappa_l*I)*DOC13l - gamma*G*d0(C2P.*alpha_dic2poc)*R13o ;  % DOCl 
    %
    % eq6 = (TRdiv+kappa_r)*DOC13r - kC*DOC13 + eta*kC*DOC13 ; % DOCr 
    %
    % F   = [eq1; eq2; eq3; eq4; eq5; eq6];
    
    % extract varying source sink terms (depend on DIC,..)
    eq1 = d0(par.Cnpp(iwet).*alpha_dic2poc)*R13o ...     % removal of dic13 organic c13 production
          + (1-sigC-gamma)*RR*G*d0(C2P)*R13o  ... % removal of dic13 due to pic13 production
          - JgDIC13 ...              % air-sea gas exchange
          + sDICbar*d0(pme)*R13o;  % concentration and dillution due to precip and evaporation

    eq2 = - (1-sigC-gamma)*G*d0(C2P.*alpha_dic2poc)*R13o;   ; % FPOC

    eq3 = - sigC*G*d0(C2P.*alpha_dic2poc)*R13o  ; % FDOC

    eq4 = - (1-sigC-gamma)*RR*G*d0(C2P)*R13o ; % FPIC

    eq5 = - d0((par.Cnpp(iwet)-G*C2P).*alpha_dic2poc)*R13o ;  % DOCl 

    eq6 = zeros(nwet,1) ; % DOCr 

    f   = [eq1; eq2; eq3; eq4; eq5; eq6];
 
    % construct the LHS matrix for the offline model
    % disp('Preparing LHS and RHS matrix:')
    % colum 1 dFdDIC13
    Jc{1,1} = TRdiv; 
    Jc{2,1} = 0*I;
    Jc{3,1} = 0*I;
    Jc{4,1} = 0*I;
    Jc{5,1} = 0*I;
    Jc{6,1} = 0*I ;
    % colum 2 dFdPOC13
    Jc{1,2} = -kappa_p*I ;
    Jc{2,2} = PFDc + kappa_p*I ;
    Jc{3,2} = 0*I ;
    Jc{4,2} = 0*I ;
    Jc{5,2} = 0*I ;
    Jc{6,2} = 0*I ;
    % colum 3 dFdDOC13
    Jc{1,3} = -eta*kC ;
    Jc{2,3} = 0*I ;
    Jc{3,3} = TRdiv + kC ;
    Jc{4,3} = 0*I ;
    Jc{5,3} = 0*I ;
    Jc{6,3} = -kC*I + eta*kC*I;
    % colum 4 dFdPIC13
    Jc{1,4} = -kPIC*I ;
    Jc{2,4} = 0*I ;
    Jc{3,4} = 0*I ;
    Jc{4,4} = PFDa + kPIC*I ;
    Jc{5,4} = 0*I ;
    Jc{6,4} = 0*I;
    % column 6 dFdDOC13l
    Jc{1,5} = -kappa_l*I ;
    Jc{2,5} = 0*I ;
    Jc{3,5} = 0*I ;
    Jc{4,5} = 0*I ;
    Jc{5,5} = TRdiv + kappa_l*I ;
    Jc{6,5} = 0*I ;
    % column 6 dFdDOC13r
    Jc{1,6} = -kappa_r;
    Jc{2,6} = 0*I ;
    Jc{3,6} = 0*I ;
    Jc{4,6} = 0*I ;
    Jc{5,6} = 0*I ;
    Jc{6,6} = TRdiv + kappa_r ;
    J = cell2mat(Jc);
    % if nargout > 1
    %   % factorize Jacobian matrix
    %   FD = mfactor(cell2mat(Jc)) ;
    % end 
    
end

