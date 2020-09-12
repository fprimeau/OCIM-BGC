function par = SetPara(par)
    spd  = 24*60^2;
    spa  = 365*spd;

    % fixed parameters 
    par.kappa_g  = 1/(1e6*spa)  ; % geological restoring time [1/s];
    par.taup     = 720*60^2     ; % (s) pic dissolution time-scale
    par.tau_TA   = 1./par.taup  ;
    par.kappa_p  = 1/(720*60^2) ;
    % load optimal parameters if they exist
    % if isfile(par.fxhat)
    % load(par.fxhat)
    % end

    if exist('xhat') & isfield(xhat,'sigma')
        par.sigma = xhat.sigma ;
    else 
        par.sigma = 0.75 ;
    end
    if exist('xhat') & isfield(xhat,'kP_T')
        par.kP_T = xhat.kP_T ;
    else 
        par.kP_T = 1.35 ;
    end 
    if exist('xhat') & isfield(xhat,'kdP')
        par.kdP = xhat.kdP ;
    else 
        par.kdP = 1.01e-07 ; 
    end 
    if exist('xhat') & isfield(xhat,'bP_T')
        par.bP_T = xhat.bP_T ;
    else 
        par.bP_T = 6.36e-02 ;
    end 
    if exist('xhat') & isfield(xhat,'bP')
        par.bP  = xhat.bP ;
    else 
        par.bP  = 1.03 ;
    end 
    if exist('xhat') & isfield(xhat,'alpha')
        par.alpha = xhat.alpha ;
    else 
        par.alpha = 1.95e-07   ;
    end 
    if exist('xhat') & isfield(xhat,'beta')
        par.beta = xhat.beta ;
    else
        par.beta = 3.59e-02  ;
    end 

    % C model parameters                                      
    if exist('xhat') & isfield(xhat,'bC_T')
        par.bC_T = xhat.bC_T ;
    else
        par.bC_T =  8.11e-02 ; 
    end 
    if exist('xhat') & isfield(xhat,'bC')
        par.bC = xhat.bC ;
    else
        par.bC = 1.40    ;
    end 
    if exist('xhat') & isfield(xhat,'d')
        par.d = xhat.d   ;
    else
        par.d = 3.73e+03 ;
    end 
    if exist('xhat') & isfield(xhat,'kC_T')
        par.kC_T = xhat.kC_T;
    else 
        par.kC_T = 2.27 ;
    end 
    if exist('xhat') & isfield(xhat,'kdC')
        par.kdC = xhat.kdC ;
    else 
        par.kdC = 5.05e-08 ;
    end 
    if exist('xhat') & isfield(xhat,'RR')
        par.RR = xhat.RR  ;
    else
        par.RR = 8.80e-03    ;
    end
    if exist('xhat') & isfield(xhat,'cc')
        par.cc = xhat.cc  ;
    else
        par.cc = 1.0e-03 ;
    end 
    if exist('xhat') & isfield(xhat,'dd')
        par.dd = xhat.dd  ;
    else 
        par.dd = 5.35e-03 ;
    end 
    %
    % O model parameters
    if exist('xhat') & isfield(xhat,'O2C_T')
        par.O2C_T = xhat.O2C_T ;
    else 
        par.O2C_T = 0.0e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'rO2C')
        par.rO2C = xhat.rO2C ;
    else 
        par.rO2C = 1.1 ;
    end 
    if exist('xhat') & isfield(xhat,'O2P_T')
        par.O2P_T = xhat.O2P_T ;
    else 
        par.O2P_T = 0.0e+00 ;
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

