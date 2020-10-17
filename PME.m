function [modT,modS,pme] = PME(par)
% precipitation minus evaporation function that is used to
% creat DIC and ALk virtual flux ;
    M3d   = par.M3d  ;
    grd   = par.grd  ;
    iwet  = par.iwet ;
    nwet  = par.nwet ;
    TRdiv = par.TRdiv;
    dVt   = par.dVt  ;
    I     = speye(nwet);   % make an identity matrix;

    tmp   = M3d ;
    tmp(:,:,2:end) = 0 ;
    isrf = find(tmp(iwet)) ; 
    Isrf = d0(tmp(iwet))   ;

    B  = TRdiv+Isrf./par.taup;
    %
    fprintf('Factoring the big matrix for surface salinity 1...'); tic
    FB = mfactor(B);
    toc;
    %
    modS = M3d+nan;
    fprintf('Factoring the big matrix for surface salinity 2...'); tic
    modS(iwet) = mfactor(FB,Isrf*par.Salt(iwet)./par.taup);
    toc
    %
    modT = M3d+nan;
    fprintf('Factoring the big matrix for surface temperature 1...'); tic
    modT(iwet) = mfactor(FB,Isrf*par.Temp(iwet)./par.taup);
    toc
    fprintf('\n')
    clear memory
    clear B 
    %
    msk = M3d;
    msk(:,:,2:end) = 0;
    sg   = sum(msk(iwet).*dVt(iwet).*modS(iwet))/sum(msk(iwet).*dVt(iwet));
    test = (modS(iwet)-par.Salt(iwet))/sg;
    pme  = msk(iwet).*test.*(1/par.taup);
    for k = 1:10
        pme = pme-msk(iwet).*sum(dVt(iwet).*pme)/sum(msk(iwet).*dVt(iwet));
    end
end

