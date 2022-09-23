function par = SetPar(par)
    on   = true    ;
    off  = false   ;
    sph  = 60^2    ;
    spd  = 24*sph  ;
    spa  = 365*spd ;

    % fixed parameters 
    par.kappa_g  = 1/(1e6*spa)  ; % geological restoring time [1/s] ;
    par.kappa_l  = 1/(12*sph )  ; % labile DOM remi time [1/s]     ;
    par.taup     = 30*spd     ; % (s) pic dissolution time-scale ;
    par.kappa_p  = 1/par.taup ;

    par.gamma = 1/3 ;
    
    % PIC dissolution constant 0.38 day^-1 based on first-order
    % reaction kinetics according to Sarmiento
    % and Gruber book (p.271);
    par.tauPIC = 30*spd ; 
    par.kPIC   = 1/par.tauPIC ;
    % load optimal parameters if they exist
    fprintf(' load optimal parameters if they exist.\n');
    if isfile(par.fxhat) & par.LoadOpt == on 
      fprintf('Loading pre-optimized parameters.\n');
      load(par.fxhat)    
    end
    
    if exist('xhat') & isfield(xhat,'sigP')
        par.sigP = xhat.sigP ;
    else 
        par.sigP = 1.00e-1 ;
    end
    if exist('xhat') & isfield(xhat,'Q10P')
        par.Q10P = xhat.Q10P ;
    else 
        par.Q10P = 1 ; 
    end 
    if exist('xhat') & isfield(xhat,'kdP')
        par.kdP = xhat.kdP ;
    else 
        par.kdP = 9.4e-08; % from N nature paper, same as kdP
    end 
    if exist('xhat') & isfield(xhat,'bP_T')
        par.bP_T = xhat.bP_T ;
    else 
        par.bP_T = 0 ;
    end 
    if exist('xhat') & isfield(xhat,'bP')
        par.bP  = xhat.bP ;
    else 
        par.bP  = 1.16e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'alpha')
        par.alpha = xhat.alpha ;
    else 
        par.alpha = 1.47e-08   ;
    end 
    if exist('xhat') & isfield(xhat,'beta')
        par.beta = xhat.beta ;
    else
        par.beta = 1.00e+00  ;
    end 

    % C model parameters
    if exist('xhat') & isfield(xhat,'sigC')
        par.sigC = xhat.sigC ;
    else 
        par.sigC = 1.00e-1 ;
    end
    if exist('xhat') & isfield(xhat,'kru')
        par.kru = xhat.kru ;
    else
        par.kru = 4.04e-12 ; % corresponding to 2000 years.
    end
    if exist('xhat') & isfield(xhat,'krd')
        par.krd = xhat.krd ;
    else
        par.krd = 4.04e-12 ; % corresponding to 2000 years.
    end
    if exist('xhat') & isfield(xhat,'etau')
        par.etau = xhat.etau ;
    else
        par.etau = 0.9973 ; 
    end 
    if exist('xhat') & isfield(xhat,'etad')
        par.etad = xhat.etad ;
    else
        par.etad = 0.99 ; 
    end 
    if exist('xhat') & isfield(xhat,'bC_T')
        par.bC_T = xhat.bC_T ;
    else
        par.bC_T =  0.00e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'bC')
        par.bC = xhat.bC ;
    else
        par.bC = 9.94e-01    ;
    end
    if exist('xhat') & isfield(xhat,'d')
        par.d = xhat.d   ;
    else
        par.d = 4.56e+03 ;
    end 
    if exist('xhat') & isfield(xhat,'Q10C')
        par.Q10C = xhat.Q10C ;
    else 
        par.Q10C = 1.00e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'kdC')
        par.kdC = xhat.kdC ;
    else 
        par.kdC =  5.7e-08; % from N nature paper,same as kdN;
    end
    if exist('xhat') & isfield(xhat,'R_Si')
        par.R_Si = xhat.R_Si ;
    else
        par.R_Si = 0.00 ;
    end
    if exist('xhat') & isfield(xhat,'rR')
        par.rR = xhat.rR  ;
    else
        par.rR = 4.23e-02 ;
    end
    if exist('xhat') & isfield(xhat,'cc')
        par.cc = xhat.cc  ;
    else
        par.cc = 6.9e-3 ;
    end 
    if exist('xhat') & isfield(xhat,'dd')
        par.dd = xhat.dd  ;
    else 
        par.dd = 6.0e-03 ;
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
        par.rO2C = 1.10e+00 ;
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

