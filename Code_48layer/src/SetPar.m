function par = SetPar(par)
    on   = true    ;
    off  = false   ;
    sph  = 60^2    ;
    spd  = 24*sph  ;
    spa  = 365*spd ;

    % fixed parameters 
    par.kappa_g  = 1/(1e6*spa)  ; % geological restoring time [1/s] ;
    par.kappa_l  = 1/(12*sph )  ; % labile DOM remi time [1/s]     ;
    par.taup     = 30*spd       ; % (s) pic dissolution time-scale ;
    par.gamma    = 0 ; % 2/5 ;
    par.kappa_p  = 1/par.taup ;
    % PIC dissolution constant 0.38 day^-1 based on first-order
    % reaction kinetics according to Sarmiento
    % and Gruber book (p.271);
    par.tauPIC = 30*spd ; 
    par.kPIC   = 1/par.tauPIC ;
    % load optimal parameters if they exist
    if isfile(par.fxhat) & par.LoadOpt == on 
        load(par.fxhat)
    end
    
    if exist('xhat') & isfield(xhat,'sigP')
        par.sigP = xhat.sigP ;
    else 
        par.sigP = 3.500e-1 ;
    end
    if exist('xhat') & isfield(xhat,'Q10P')
        par.Q10P = xhat.Q10P ;
    else 
        par.Q10P = 2.28 ; 
    end 
    if exist('xhat') & isfield(xhat,'kdP')
        par.kdP = xhat.kdP ;
    else 
        par.kdP = 1.86e-08; 
    end 
    if exist('xhat') & isfield(xhat,'bP_T')
        par.bP_T = xhat.bP_T ;
    else 
        par.bP_T = 0.319124 ;
    end 
    if exist('xhat') & isfield(xhat,'bP')
        par.bP  = xhat.bP ;
    else 
        par.bP  = 0.7701130e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'alpha')
        par.alpha = xhat.alpha ;
    else 
        par.alpha = 2.45e-08   ;
    end 
    if exist('xhat') & isfield(xhat,'beta')
        par.beta = xhat.beta ;
    else
        par.beta = 4.50e-01  ;
    end 

    % C model parameters
    if exist('xhat') & isfield(xhat,'sigC')
        par.sigC = xhat.sigC ;
    else 
        par.sigC = 0.90e-1 ;
    end
    if exist('xhat') & isfield(xhat,'kru')
        par.kru = xhat.kru ;
    else
        par.kru = 5.74e-12 ; % corresponding to 2000 years.
    end
    if exist('xhat') & isfield(xhat,'krd')
        par.krd = xhat.krd ;
    else
        par.krd = 2.99e-12 ; % corresponding to 2000 years.
    end
    if exist('xhat') & isfield(xhat,'etau')
        par.etau = xhat.etau ;
    else
        par.etau = 0.9800 ; 
    end 
    if exist('xhat') & isfield(xhat,'etad')
        par.etad = xhat.etad ;
    else
        par.etad = 0.978911524882553 ; 
    end 
    if exist('xhat') & isfield(xhat,'bC_T')
        par.bC_T = xhat.bC_T ;
    else
        par.bC_T =  0.957372e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'bC')
        par.bC = xhat.bC ;
    else
        par.bC = 6.50e-01    ;
    end
    if exist('xhat') & isfield(xhat,'d')
        par.d = xhat.d   ;
    else
        par.d = 4.55e+03 ;
    end 
    if exist('xhat') & isfield(xhat,'Q10C')
        par.Q10C = xhat.Q10C ;
    else 
        par.Q10C = 1.05e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'kdC')
        par.kdC = xhat.kdC ;
    else 
        par.kdC =  5.42e-09; % from N nature paper,same as kdN;
    end
    if exist('xhat') & isfield(xhat,'R_Si')
        par.R_Si = xhat.R_Si ;
    else
        par.R_Si = 0.10 ;
    end
    if exist('xhat') & isfield(xhat,'rR')
        par.rR = xhat.rR  ;
    else
        par.rR = 2.34e-02 ;
    end
    if exist('xhat') & isfield(xhat,'cc')
        par.cc = xhat.cc  ;
    else
        par.cc = 8.38e-4 ;
    end 
    if exist('xhat') & isfield(xhat,'dd')
        par.dd = xhat.dd  ;
    else 
        par.dd = 8.83e-03 ;
    end 

    % O model parameters
    if exist('xhat') & isfield(xhat,'O2C_T')
        par.O2C_T = xhat.O2C_T ;
    else 
        par.O2C_T = 0.00 ;
    end 
    if exist('xhat') & isfield(xhat,'rO2C')
        par.rO2C = xhat.rO2C ;
    else 
        par.rO2C = 1.77e+00 ;
    end 
    %
    % Si model parameters
    if exist('xhat') & isfield(xhat,'dsi')
        par.dsi = xhat.dsi ;
    else
        par.dsi = 3300     ;
    end 
    if exist('xhat') & isfield(xhat,'at')
        par.at = xhat.at   ;
    else
        par.at = 1.32e16/spd;
    end 
    if exist('xhat') & isfield(xhat,'bt')
        par.bt = xhat.bt   ;
    else 
        par.bt = 11481     ;
    end 
    if exist('xhat') & isfield(xhat,'aa')
        par.aa = xhat.aa   ;
    else
        par.aa = 1         ;
    end 
    if exist('xhat') & isfield(xhat,'bb')
        par.bb = xhat.bb   ;
    else 
        par.bb = 0.968     ;
    end 
end

