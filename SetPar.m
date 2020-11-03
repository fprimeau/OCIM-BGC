function par = SetPar(par)
    on   = true    ;  off  = false   ;
    spd  = 24*60^2 ;  spa  = 365*spd ;
    % fixed parameters 
    par.kappa_g  = 1/(1e6*spa)  ; % geological restoring time [1/s];
    par.taup     = 720*60^2     ; % (s) pic dissolution time-scale
    % par.tau_TA   = 1./par.taup  ;
    par.kappa_p  = 1/(720*60^2) ;
    % PIC dissolution constant 0.38 day^-1 based on first-order
    % reaction kinetics according to Sarmiento
    % and Gruber book (p.271);
    par.tauPIC = 30*spd ; 
    par.kPIC   = 1/par.tauPIC ;
    % load optimal parameters if they exist
    if isfile(par.fxhat) & par.LoadOpt == on 
        load(par.fxhat)
    end
    
    if exist('xhat') & isfield(xhat,'sigma')
        par.sigma = xhat.sigma ;
    else 
        par.sigma = 1.0e-01 ;
    end
    if exist('xhat') & isfield(xhat,'kP_T')
        par.kP_T = xhat.kP_T ;
    else 
        par.kP_T = 1.18 ;
    end 
    if exist('xhat') & isfield(xhat,'kdP')
        par.kdP = xhat.kdP ;
    else 
        par.kdP = 2.42e-08 ; 
    end 
    if exist('xhat') & isfield(xhat,'bP_T')
        par.bP_T = xhat.bP_T ;
    else 
        par.bP_T = 0.00e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'bP')
        par.bP  = xhat.bP ;
    else 
        par.bP  = 9.60e-01 ;
    end 
    if exist('xhat') & isfield(xhat,'alpha')
        par.alpha = xhat.alpha ;
    else 
        par.alpha = 3.37e-08   ;
    end 
    if exist('xhat') & isfield(xhat,'beta')
        par.beta = xhat.beta ;
    else
        par.beta = 1.65e-02  ;
    end 

    % C model parameters                                      
    if exist('xhat') & isfield(xhat,'bC_T')
        par.bC_T = xhat.bC_T ;
    else
        par.bC_T =  0.00e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'bC')
        par.bC = xhat.bC ;
    else
        par.bC = 9.38e-01    ;
    end
    if exist('xhat') & isfield(xhat,'d')
        par.d = xhat.d   ;
    else
        par.d = 4.54e+03 ;
    end 
    if exist('xhat') & isfield(xhat,'kC_T')
        par.kC_T = xhat.kC_T ;
    else 
        par.kC_T = 0.00e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'kdC')
        par.kdC = xhat.kdC ;
    else 
        par.kdC = 7.35e-08 ;
    end
    if exist('xhat') & isfield(xhat,'R_Si')
        par.R_Si = xhat.R_Si ;
    else
        par.R_Si = 0.00 ;
    end
    if exist('xhat') & isfield(xhat,'rR')
        par.rR = xhat.rR  ;
    else
        par.rR = 2.64e-02 ;
    end
    if exist('xhat') & isfield(xhat,'cc')
        par.cc = xhat.cc  ;
    else
        par.cc = 7.51e-4 ;
    end 
    if exist('xhat') & isfield(xhat,'dd')
        par.dd = xhat.dd  ;
    else 
        par.dd = 5.56e-03 ;
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
    if exist('xhat') & isfield(xhat,'O2P_T')
        par.O2P_T = xhat.O2P_T ;
    else 
        par.O2P_T = 0.00 ;
    end 
    if exist('xhat') & isfield(xhat,'rO2P')
        par.rO2P = xhat.rO2P ;
    else 
        par.rO2P = 1.70e+02 ;
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

