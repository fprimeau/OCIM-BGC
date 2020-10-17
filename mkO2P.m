function vout = mkO2P(par)
% parameters
    on   = true    ; off = false    ;
    if par.opt_O2P_T == on 
        O2P_T = par.O2P_T ;
    else
        O2P_T = 0 ;
    end 
    rO2P  = par.rO2P  ;
    %
    iwet = par.iwet ;
    nwet = par.nwet ;
    modT = par.modT ;
    
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
end