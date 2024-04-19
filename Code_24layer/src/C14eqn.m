function [f,J,par] = C14eqn(X, par)    
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
    DIC14  = X(0*nwet+1:1*nwet) ; 
    POC14  = X(1*nwet+1:2*nwet) ;
    DOC14  = X(2*nwet+1:3*nwet) ;
    PIC14  = X(3*nwet+1:4*nwet) ;
    % ALK  = X(4*nwet+1:5*nwet) ;
    DOC14l = X(4*nwet+1:5*nwet) ;
    DOC14r = X(5*nwet+1:6*nwet) ;

    R14o = DIC14./(par.DIC); 
    dR14o = d0(1./par.DIC);

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
    par.DIC14  = DIC14  ;
    
    % define C14 fractionation factors and ratios in ocean and atmosphere
    % set up fractionation factors fro C14 and R14
    %  A. Schmittner et al.: Distribution of carbon isotope ratios (δ13C) in the ocean
    % fc14 = 2.0;  %  The fractionation  for 14C is 2.3 times the fractionation for 13C
    fc14 = par.fc14;
    lambda14 = 1/par.spa*log(2)/5730; % radiocarbon decay rate (yr^(-1) to s^(-1))
    % R14a is changing with time and passed in in par.c14
    % par.c14.R14a = 1.220805*1e-12; % air C14/C Roxa = 1.176*e-12
    par.pc14atm = par.pco2atm*par.c14.R14a;
    par.c14.R14o = R14o;
    par.c14.dR14o = dR14o;
    par.c14.fc14 = fc14;
    par.c14.alpha_k = 1 - (1 - 0.99915)*fc14 ; % kenetic fractionation factor 
    par.c14.alpha_g2aq = 1- (1 - 0.998764)*fc14 ; % gas to water fractionation factor
    % temperature (in C)-dependent equilibrium fractionation factor from gaseous CO2 to DIC.
    % par.c14.alpha_g2dic = 1.01051-1.05*1e-4*par.Temp(isrf); 
    disp(sprintf('air-sea fractionation tuned by fras %4.2f and fc14 %4.2f',par.fras,fc14));
    par.c14.alpha_g2dic = (1.01051-1.05*1e-4*par.Temp-1.0)*par.fras*fc14 + 1.0;
    % par.c14.alpha_g2dic = (1.01051-1.05*1e-4*par.Temp - 1.0)*fc14 + 1.0; 
    
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
    
    alpha_aq2poc = (-0.017*log(co2) + 1.0034 -1.0)*fc14 + 1.0; % check the unit of co2surf
    % alpha_dic2poc = ( par.c14.alpha_g2aq ./ par.c14.alpha_g2dic ) .* alpha_aq2poc; 
    alpha_tmp = ( par.c14.alpha_g2aq ./ par.c14.alpha_g2dic ) .* alpha_aq2poc; 
    alpha_dic2poc = alpha_tmp(iwet);
    % alpha_dic2poc = (alpha_tmp(iwet)-1.0)*0.60 + 1.0;
    alpha_dic2poc = (alpha_tmp(iwet)-1.0)*par.frpho*fc14 + 1.0;

    if par.debug14
      disp('Testing fractionation factors')
      par.c14.alpha_k = 1;%0.99915; % kenetic fractionation factor 
      par.c14.alpha_g2aq = 1;%0.998764; % gas to water fractionation factor
      par.c14.alpha_g2dic = 1.01051*0 - 1.05*1e-4*par.Temp*0 + 1;
      alpha_dic2poc = alpha_dic2poc*0 + 1.0; % debug
    end
	
    % Air-Sea gas exchange for C14
    vout  = Fsea2air(par, 'C14');
    JgDIC14 = vout.JgDIC14;
    G_dic14 = vout.G_dic14;
    rhs14 = vout.rhs14;
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

    % eq1 = TRdiv*DIC14 ...                         % advective-diffusive transport
    %       + G*d0(C2P.*alpha_dic2poc)*R14o ...     % removal of dic14 organic c14 production
    %       + (1-sigC-gamma)*RR*G*d0(C2P)*R14o  ... % removal of dic14 due to pic14 production
    %       - kPIC*PIC14 ....          % dissolution of PIC14
    %       - JgDIC14 ...              % air-sea gas exchange
    %       + sDICbar*d0(pme)*R14o ... % concentration and dillution due to precip and evaporation
    %       - eta*(kC*DOC14) ...       % respiration of DOC14
    %       - kappa_r*DOC14r ...       % respiration of DOC14r
    %       - kappa_l*DOC14l...        % respiration of DOC14l
    %       - kappa_p*POC14 ...        % respiration of POC14
    %       + lambda14*DIC14 ;         % decay of DIC14
    %
    % eq2 = (PFDc+kappa_p*I)*POC14 - (1-sigC-gamma)*G*d0(C2P.*alpha_dic2poc)*R14o + lambda14*POC14;   ; % FPOC
    %
    % eq3 = (TRdiv+kC)*DOC14 - sigC*G*d0(C2P.*alpha_dic2poc)*R14o + lambda14*DOC14 ; % FDOC
    %
    % eq4 = (PFDa + kPIC*I)*PIC14 - (1-sigC-gamma)*RR*G*d0(C2P)*R14o + lambda14*PIC14; % FPIC
    %
    % eq5 = (TRdiv+kappa_l*I)*DOC14l - gamma*G*d0(C2P.*alpha_dic2poc)*R14o + lambda14*DOC14l;  % DOCl 
    %
    % eq6 = (TRdiv+kappa_r)*DOC14r - kC*DOC14 + eta*kC*DOC14 + lambda14*DOC14r; % DOCr 
    %
    % F   = [eq1; eq2; eq3; eq4; eq5; eq6];
    
    % extract varying source sink terms (depend on DIC,..)
    eq1 = G*d0(C2P.*alpha_dic2poc)*R14o ...     % removal of dic14 organic c14 production
          + (1-sigC-gamma)*RR*G*d0(C2P)*R14o  ... % removal of dic14 due to pic14 production
          - JgDIC14 ...              % air-sea gas exchange
          + sDICbar*d0(pme)*R14o;  % concentration and dillution due to precip and evaporation

    eq2 = - (1-sigC-gamma)*G*d0(C2P.*alpha_dic2poc)*R14o;   ; % FPOC

    eq3 = - sigC*G*d0(C2P.*alpha_dic2poc)*R14o; % FDOC

    eq4 = - (1-sigC-gamma)*RR*G*d0(C2P)*R14o ; % FPIC

    eq5 = - gamma*G*d0(C2P.*alpha_dic2poc)*R14o ;  % DOCl 

    eq6 = zeros(nwet,1); % DOCr 

    f   = [eq1; eq2; eq3; eq4; eq5; eq6];

    % construct the LHS matrix for the offline model
    % disp('Preparing LHS and RHS matrix:')
    % colum 1 dFdDIC14
    Jc{1,1} = TRdiv + lambda14*I; 
    Jc{2,1} = 0*I;
    Jc{3,1} = 0*I;
    Jc{4,1} = 0*I;
    Jc{5,1} = 0*I;
    Jc{6,1} = 0*I ;
    % colum 2 dFdPOC14
    Jc{1,2} = -kappa_p*I ;
    Jc{2,2} = PFDc + kappa_p*I + lambda14*I;
    Jc{3,2} = 0*I ;
    Jc{4,2} = 0*I ;
    Jc{5,2} = 0*I ;
    Jc{6,2} = 0*I ;
    % colum 3 dFdDOC14
    Jc{1,3} = -eta*kC ;
    Jc{2,3} = 0*I ;
    Jc{3,3} = TRdiv + kC + lambda14*I;
    Jc{4,3} = 0*I ;
    Jc{5,3} = 0*I ;
    Jc{6,3} = -kC*I + eta*kC*I;
    % colum 4 dFdPIC14
    Jc{1,4} = -kPIC*I ;
    Jc{2,4} = 0*I ;
    Jc{3,4} = 0*I ;
    Jc{4,4} = PFDa + kPIC*I + lambda14*I;
    Jc{5,4} = 0*I ;
    Jc{6,4} = 0*I;
    % column 5 dFdDOC14l
    Jc{1,5} = -kappa_l*I ;
    Jc{2,5} = 0*I ;
    Jc{3,5} = 0*I ;
    Jc{4,5} = 0*I ;
    Jc{5,5} = TRdiv + kappa_l*I + lambda14*I;
    Jc{6,5} = 0*I ;
    % column 6 dFdDOC14r
    Jc{1,6} = -kappa_r;
    Jc{2,6} = 0*I ;
    Jc{3,6} = 0*I ;
    Jc{4,6} = 0*I ;
    Jc{5,6} = 0*I ;
    Jc{6,6} = TRdiv + kappa_r + lambda14*I ;
    J = cell2mat(Jc);
    % if nargout > 1
    %   % factorize Jacobian matrix
    %   FD = mfactor(cell2mat(Jc)) ;
    % end 
  end
