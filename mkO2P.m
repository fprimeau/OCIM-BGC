function vout = mkO2P(par)
% parameters
    if par.opt_O2P_T == on 
        O2P_T = par.O2P_T ;
    else
        O2P_T = 0 ;
    end 
    rO2P  = par.rO2P  ;
    %
    iwet = par.iwet ;
    nwet = par.nwet ;
    sal  = par.modS ;
    modT = par.modT ;
    SIL  = par.SIL  ;
    smsk = par.M3d  ;
    smsk(:,:,2:end) = 0 ;
    isrf = find(smsk(iwet))    ;
    dVs  = par.dVt(iwet(isrf)) ;
    surface_mean = @(x) sum(x(isrf).*dVs)/sum(dVs) ;

    % compute the mean of the regressor variable
    Z  = modT(iwet)      ;
    mu = surface_mean(Z) ;
    Delta = sqrt(surface_mean((Z-mu).^2)) ;

    % standardize the regressor variables
    ZR = (Z - mu) / Delta ; par.ZR = ZR ;

    vout.O2P = O2P_T*ZR + rO2P ;
    vout.dO2PdO2P_T = ZR ;
    vout.dO2PdrO2P  = sparse(nwet,1) + 1 ;

    % build CaCO3 to production ratio;
    % According to Sarmiento and Gruber, this ratio has a
    % inverse relationship with silicon acid (~0.02 in high latitude
    % and 0.09-0.1 in low latitude oceans)
    Y = 0.5-0.5*tanh((SIL(iwet) - 30)/100) ;
    if opt_R_Si == on
        R_Si = par.R_Si ;
    else 
        R_Si = 0 ;
    end 
    
    vout.RR = R_Si*Y + rR ;
    vout.dRRdSi = Y ;
    vout.dRRdrR = sparse(nwet,1) + 1 ;
end