function vout = mkPIC2P(par)

% build CaCO3 to production ratio;
% According to Sarmiento and Gruber, this ratio has a
% inverse relationship with silicon acid (~0.02 in high latitude
% and 0.09-0.1 in low latitude oceans)
    on   = true     ;
    off  = false    ;
    M3d  = par.M3d  ;
    iwet = par.iwet ;
    nwet = par.nwet ; 
    DSi  = par.DSi  ;

    smsk = par.M3d  ;
    smsk(:,:,3:end) = 0 ;
    isrf = find(smsk(iwet))    ;

    Y = M3d(iwet)*0 ;
    sDSi = DSi(iwet(isrf));
    Y(isrf) =  (sDSi - min(sDSi))./(max(sDSi) - min(sDSi));

    % Y = 0.5-0.5*tanh((DSi(iwet) - 30)/100) ;
    if isfield(par,'opt_R_Si') & par.opt_R_Si == on
        R_Si = par.R_Si ;
    else 
        R_Si = 0 ;
    end 
    rR         = par.rR ;
    vout.RR    = d0(R_Si*Y + rR) ;
    vout.RR_Si = d0(Y)  ;
    vout.RR_rR = d0(sparse(nwet,1) + rR) ;
end 