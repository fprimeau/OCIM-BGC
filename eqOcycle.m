function [par, O2, Ox, Oxx] = eqOcycle(x, par)
    on = true; off = false;
    global GO
    pindx = par.pindx ;
    % unpack the parameters to be optimized
    % O2C_T
    if (par.opt_O2C_T == on)
        par.O2C_T = x(pindx.O2C_T) ;
    end
    
    % rO2C
    if (par.opt_rO2C == on)
        lrO2C    = x(pindx.lrO2C) ;
        par.rO2C = exp(lrO2C)     ;
    end
    % O2P_T
    if (par.opt_O2P_T == on)
        par.O2P_T = x(pindx.O2P_T) ;
    end
    
    % rO2P
    if (par.opt_rO2P == on)
        lrO2P    = x(pindx.lrO2P) ;
        par.rO2P = exp(lrO2P)     ;
    end
    %
    X0  = GO;

    options.iprint = 0 ;
    options.atol   = 1e-10 ;
    options.rtol   = 1e-10 ;

    fprintf('Solving O model ...\n') ;
    [O2,ierr] = nsnew(X0,@(X) O_eqn(X, par),options) ;
    if (ierr ~= 0)
        fprintf('O2model did not converge.\n') ;

        F = O_eqn(O2, par) ;

        if (norm(F) < 1e-6)  
            % Compute the gradient of the solution wrt the parameters
            GO = real(O2) + 1e-7*randn(par.nwet,1) ;
            [F, FD, Ox, Oxx] = O_eqn(O2, par) ;
        else 
            npx  = par.npx ;
            ncx  = par.ncx ;
            nox  = par.nox ;
            nx   = npx + ncx + nox  ;
            Ox   = sparse(par.nwet, nx) ;
            Oxx  = sparse(par.nwet, nchoosek(nx,2)+nx) ;
        end 
    else
        % reset the global variable for the next call eqOcycle
        GO = real(O2) + 1e-7*randn(par.nwet,1) ;
        [F, FD, Ox, Oxx] = O_eqn(O2, par) ;
    end
end

function [F, FD, Ox, Oxx] = O_eqn(O2, par)
    on = true; off = false;
    pindx = par.pindx ;
    %
    % fixed parameters
    iwet  = par.iwet   ;
    nwet  = par.nwet   ;
    TRdiv = par.TRdiv ;
    I     = speye(nwet)   ;
    PO4   = par.po4obs(iwet) ;
    % variables from C model
    DOC  = par.DOC    ;
    Tz   = par.Tz     ;
    TZ   = par.Tz*1e8 ;
    %
    % tunable parameters;
    O2C_T = par.O2C_T   ; 
    rO2C  = par.rO2C    ;
    kC_T  = par.kC_T    ;
    kdC   = par.kdC     ;
    O2P_T = par.O2P_T  ;
    rO2P  = par.rO2P ;
    kC    = d0(kC_T * Tz + kdC) ; 
    %
    vout = mkO2P(par) ;
    O2P  = vout.O2P   ;
    dO2PdO2P_T = vout.dO2PdO2P_T ;
    dO2PdrO2P  = vout.dO2PdrO2P  ;
    %
    % O2 saturation concentration
    vout  = Fsea2air(par,'O2') ;
    KO2   = vout.KO2   ;
    o2sat = vout.o2sat ;
    % rate of o2 production
    G   = par.G        ;
    PO2 = G*O2P        ;

    % parobolic function for o2 consumption
    R      = 0.5 + 0.5*tanh(O2-10)    ;
    dRdO   = 0.5 - 0.5*tanh(O2-10).^2 ;
    d2RdO2 = -2*d0(dRdO)*tanh(O2-10)  ;
     
    O2C    = O2C_T*TZ + rO2C ; 
    % rate of o2 utilization
    LO2    = kC*DOC.*O2C.*R  ;
    dLdO   = d0(kC*DOC.*O2C.*dRdO) ;
    d2LdO2 = kC*DOC.*O2C.*d2RdO2   ;
    % O2 function
    F = TRdiv*O2 - PO2 + LO2 - KO2*(o2sat-O2) ;
    %
    if (nargout > 1)
        FD = mfactor(TRdiv + dLdO + KO2) ;
    end
    %
    %% ----------------- Gradient -------------------
    if (par.optim == off)
        Ox = [];
    elseif (par.optim & nargout > 2)
        npx  = par.npx ;
        ncx  = par.ncx ;
        nox  = par.nox ;
        nx   = npx + ncx + nox  ;
        Ox   = sparse(nwet, nx) ;
        Gx   = par.Gx ;
        DICx = par.Cx(1:nwet,:) ;
        DOCx = par.Cx(2*nwet+1:3*nwet,:) ;
        % P model only parameters
        for jj = 1:npx
            tmp = -d0(Gx(:,jj))*O2P + kC*DOCx(:,jj).*O2C.*R ;
            Ox(:,jj) = mfactor(FD, -tmp) ;
        end 

        % C model only parameters
        for jj = (npx+1) : (npx+ncx)
            if (par.opt_kC_T & jj == pindx.kC_T)
                kC_kC_T = par.kC_kC_T ; 
                tmp = kC_kC_T*DOC.*O2C.*R + ...
                      kC*DOCx(:,jj).*O2C.*R  ;
                Ox(:,jj) = mfactor(FD, -tmp) ;

            elseif (par.opt_kdC & jj == pindx.lkdC)
                kC_kdC = par.kC_kdC ;
                tmp = kC_kdC*DOC.*O2C.*R + ...
                      kC*DOCx(:,jj).*O2C.*R  ;
                Ox(:,jj) = mfactor(FD, -tmp) ;

            else
                tmp = kC*DOCx(:,jj).*O2C.*R  ;
                Ox(:,jj) = mfactor(FD, -tmp) ;
            end 
        end

        % O model only parameters
        if (par.opt_O2C_T == on)
            O2C_O2C_T = d0(TZ) ; 
            tmp = kC*DOC.*(O2C_O2C_T*R) ; 
            Ox(:,pindx.O2C_T) = mfactor(FD, -tmp) ;
        end

        if (par.opt_rO2C == on)
            O2C_rO2C = rO2C ; 
            tmp = kC*DOC.*(O2C_rO2C*R) ; 
            Ox(:,pindx.lrO2C) = mfactor(FD, -tmp) ;
        end 
        
        if (par.opt_O2P_T == on) 
            tmp = -G*dO2PdO2P_T ;  
            Ox(:,pindx.O2P_T) = mfactor(FD, -tmp) ;
        end

        if (par.opt_rO2P == on) 
            tmp = -rO2P*diag(G) ; 
            Ox(:,pindx.lrO2P) = mfactor(FD, -tmp) ;
        end
    end

    %% ------------------- Hessian -------------------------
    if (par.optim == off)
        Oxx = [];
    elseif (par.optim & nargout > 3)
        %
        Gxx   = par.Gxx ;
        DICxx = par.Cxx(1:nwet,:) ;
        DOCxx = par.Cxx(2*nwet+1:3*nwet,:) ;
        Oxx   = sparse(nwet,nchoosek(nx,2)+nx) ;
        % P model parameters
        kk = 0 ;
        for ju = 1 : npx
            for jo = ju : npx
                kk = kk + 1 ;
                tmp = -d0(Gxx(:,kk))*O2P + ...
                      kC*DOCxx(:,kk).*R.*O2C + ...
                      kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                      kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ...
                      Ox(:,ju).*d2LdO2.*Ox(:,jo) ;

                Oxx(:,kk) = mfactor(FD, -tmp) ;                 
            end 
        end
        %
        % P C model parameters
        for ju = 1 : npx
            for jo = (npx+1) : (npx+ncx)
                kk = kk + 1;
                % check if hessian associated with kdC, if
                % yes, use formula (1), otherwise use formula (2);
                if (par.opt_kC_T & jo == pindx.kC_T)
                    tmp = kC_kC_T*DOCx(:,ju).*R.*O2C + ...
                          kC_kC_T*DOC.*dRdO.*Ox(:,ju).*O2C + ...
                          kC*R.*DOCxx(:,kk).*O2C + ...
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ...
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ; 
                    
                elseif (par.opt_kdC & jo == pindx.lkdC)
                    tmp = kC_kdC*R.*DOCx(:,ju).*O2C + ...
                          kC_kdC*DOC.*dRdO.*Ox(:,ju).*O2C + ...
                          kC*R.*DOCxx(:,kk).*O2C + ...
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ...
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ; 
                    
                else
                    tmp = kC*R.*DOCxx(:,kk).*O2C + ...
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ; 

                end
                Oxx(:,kk) = mfactor(FD, -tmp)    ;
            end
        end
        
        % C model only parameters
        for ju = (npx+1) : (npx+ncx)
            for jo = ju : (npx+ncx)
                kk = kk + 1 ;
                if (par.opt_kC_T & ju == pindx.kC_T & jo == pindx.kC_T)
                    % kC_T kC_T
                    tmp = kC_kC_T*DOCx(:,ju).*R.*O2C*2 + ...
                          kC_kC_T*DOC.*dRdO.*Ox(:,ju).*O2C*2 + ... 
                          kC*DOCxx(:,kk).*R.*O2C + ... 
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ... 
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ; 
                    
                elseif (par.opt_kdC & ju == pindx.lkdC & jo == pindx.lkdC)
                    % kdC kdC
                    % kC_kdC_kdC = kC_kdC
                    tmp = kC_kdC*DOC.*R.*O2C + ...
                          kC_kdC*DOCx(:,ju).*R.*O2C*2 + ...
                          kC_kdC*DOC.*dRdO.*Ox(:,ju).*O2C*2 + ... 
                          kC*DOCxx(:,kk).*R.*O2C + ... 
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ... 
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ; 
                    
                    % kC_T kdC 
                elseif (par.opt_kC_T & par.opt_kdC & ju == pindx.kC_T ...
                        & jo == pindx.lkdC)
                    tmp = kC_kC_T*DOCx(:,jo).*R.*O2C + ...
                          kC_kC_T*DOC.*dRdO.*Ox(:,jo).*O2C + ...
                          kC_kdC*DOCx(:,ju).*R.*O2C + ...
                          kC_kdC*DOC.*dRdO.*Ox(:,ju).*O2C + ... 
                          kC*DOCxx(:,kk).*R.*O2C + ... 
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ... 
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ; 

                    %parameter pairs with only one kC_T
                elseif (par.opt_kC_T & (ju == pindx.kC_T | jo == ...
                                        pindx.kC_T) & ju ~= jo)
                    
                    if (ju < pindx.kC_T)
                        jk = ju ;
                    else
                        jk = jo ;
                    end
                    %
                    tmp = kC_kC_T*DOCx(:,jk).*R.*O2C + ...
                          kC_kC_T*DOC.*dRdO.*Ox(:,jk).*O2C + ... 
                          kC*DOCxx(:,kk).*R.*O2C + ... 
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ... 
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ; 

                    %parameter pairs with only one kdC
                elseif (par.opt_kdC & (ju == pindx.lkdC | jo == ...
                                       pindx.lkdC) & ju ~= jo)
                    
                    if (ju < pindx.lkdC)
                        jk = ju ;
                    else
                        jk = jo ;
                    end
                    %
                    tmp = kC_kdC*DOCx(:,jk).*R.*O2C + ...
                          kC_kdC*DOC.*dRdO.*Ox(:,jk).*O2C + ... 
                          kC*DOCxx(:,kk).*R.*O2C + ... 
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ... 
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                    
                    % parameter pairs without kC_T and kdC    
                else 
                    tmp = kC*DOCxx(:,kk).*R.*O2C + ...
                          kC*DOCx(:,jo).*dRdO.*Ox(:,ju).*O2C + ...
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                    
                end 
                Oxx(:,kk) = mfactor(FD, -tmp)    ;
            end 
        end 
        
        % P and O model only parameters
        dO2Pdp{1} = dO2PdO2P_T ;
        dO2Pdp{2} = rO2P*dO2PdrO2P ;
        for ju = 1 : npx
            for jo = (npx+ncx+1) : (npx+ncx+nox)
                kk = kk + 1 ;
                if (par.opt_O2C_T == on & jo == pindx.O2C_T)
                    tmp = kC*DOC.*(O2C_O2C_T*dRdO.*Ox(:,ju)) + ...
                          kC*DOCx(:,ju).*(O2C_O2C_T*R) + ...
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ;

                elseif (par.opt_rO2C == on & jo == pindx.lrO2C)
                    tmp = kC*DOC.*(O2C_rO2C*dRdO.*Ox(:,ju)) + ...
                          kC*DOCx(:,ju).*(O2C_rO2C*R) + ...
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ;

                elseif (par.opt_O2P_T == on & jo == pindx.O2P_T)
                    tmp = -d0(Gx(:,ju))*dO2Pdp{1} + ...
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ;

                elseif (par.opt_rO2P == on & jo == pindx.lrO2P)
                    tmp = -d0(Gx(:,ju))*dO2Pdp{2} + ...
                          kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                    
                end
                Oxx(:,kk) = mfactor(FD, -tmp)    ;
            end 
        end
        
        % C and O model only parameters
        for ju = (npx+1) : (npx+ncx)
            for jo = (npx+ncx+1) : (npx+ncx+nox)
                kk = kk + 1 ;
                if (par.opt_kC_T == on & ju == pindx.kC_T)
                    if (par.opt_O2C_T == on & jo == pindx.O2C_T)
                        tmp = kC_kC_T*DOC.*dRdO.*Ox(:,jo).*O2C + ...
                              kC_kC_T*DOC.*(O2C_O2C_T*R) + ...
                              kC*DOC.*(O2C_O2C_T*dRdO.*Ox(:,ju)) + ...
                              kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              kC*DOCx(:,ju).*(O2C_O2C_T*R) + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                        %
                    elseif (par.opt_rO2C == on & jo == pindx.lrO2C)
                        tmp = kC_kC_T*DOC.*dRdO.*Ox(:,jo).*O2C + ...
                              kC_kC_T*DOC.*(O2C_rO2C*R) + ...
                              kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              kC*DOCx(:,ju).*(O2C_rO2C*R) + ...
                              kC*DOC.*(O2C_rO2C*dRdO.*Ox(:,ju)) + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                        %
                    else 
                        tmp = kC_kC_T*DOC.*dRdO.*Ox(:,jo).*O2C + ...
                              kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                    end
                        %
                elseif (par.opt_kdC == on & ju == pindx.lkdC)
                    if (par.opt_O2C_T == on & jo == pindx.O2C_T)
                        tmp = kC_kdC*DOC.*dRdO.*Ox(:,jo).*O2C + ...
                              kC_kdC*DOC.*(O2C_O2C_T*R) + ...
                              kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              kC*DOCx(:,ju).*(O2C_O2C_T*R) + ...
                              kC*DOC.*(O2C_O2C_T*dRdO.*Ox(:,ju)) + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                        %
                    elseif (par.opt_rO2C == on & jo == pindx.lrO2C)
                        tmp = kC_kdC*DOC.*dRdO.*Ox(:,jo).*O2C + ...
                              kC_kdC*DOC.*(O2C_rO2C*R) + ...
                              kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              kC*DOCx(:,ju).*(O2C_rO2C*R) + ...
                              kC*DOC.*(O2C_rO2C*dRdO.*Ox(:,ju)) + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                        %
                    else 
                        tmp = kC_kdC*DOC.*dRdO.*Ox(:,jo).*O2C + ...
                              kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                    end
                    %
                else
                    if (par.opt_O2C_T & jo == pindx.O2C_T)
                        tmp = kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              kC*DOCx(:,ju).*(O2C_O2C_T*R) + ...
                              kC*DOC.*(O2C_O2C_T*dRdO.*Ox(:,ju)) + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                        
                    elseif (par.opt_rO2C & jo == pindx.lrO2C)
                        tmp = kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              kC*DOCx(:,ju).*(O2C_rO2C*R) + ...
                              kC*DOC.*(O2C_rO2C*dRdO.*Ox(:,ju)) + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;

                    else 
                        tmp = kC*DOCx(:,ju).*dRdO.*Ox(:,jo).*O2C + ...
                              Ox(:,ju).*d2LdO2.*Ox(:,jo) ;
                        
                    end 
                end
                Oxx(:,kk) = mfactor(FD, -tmp) ;
            end 
        end 
        
        % O model only parameters
        % O2C_T O2C_T
        if (par.opt_O2C_T == on)
            kk = kk + 1 ;
            tmp = kC*DOC.*(O2C_O2C_T*(dRdO.*Ox(:,pindx.O2C_T)))*2 + ... 
                  Ox(:,pindx.O2C_T).*d2LdO2.*Ox(:,pindx.O2C_T) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end             
        % O2C_T rO2C
        if (par.opt_O2C_T == on & par.opt_rO2C == on)
            kk = kk + 1 ;
            tmp = kC*DOC.*(O2C_O2C_T*dRdO.*Ox(:,pindx.lrO2C)) + ...
                  kC*DOC.*(O2C_rO2C*dRdO.*Ox(:,pindx.O2C_T)) + ... 
                  Ox(:,pindx.O2C_T).*d2LdO2.*Ox(:,pindx.lrO2C) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end             
        % O2C_T O2P_T
        if (par.opt_O2C_T == on & par.opt_O2P_T == on)
            kk = kk + 1 ;
            tmp = kC*DOC.*(O2C_O2C_T*dRdO.*Ox(:,pindx.O2P_T)) + ...
                  Ox(:,pindx.O2C_T).*d2LdO2.*Ox(:,pindx.O2P_T) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end             
        % O2C_T rO2P
        if (par.opt_O2C_T == on & par.opt_rO2P == on)
            kk = kk + 1 ;
            tmp = kC*DOC.*(O2C_O2C_T*dRdO.*Ox(:,pindx.lrO2P)) + ...
                  Ox(:,pindx.O2C_T).*d2LdO2.*Ox(:,pindx.lrO2P) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end             
        % rO2C rO2C
        if (par.opt_rO2C == on & par.opt_rO2C == on)
            kk = kk + 1 ;
            tmp = kC*DOC.*(O2C_rO2C*R) + ...
                  kC*DOC.*(O2C_rO2C*dRdO.*Ox(:,pindx.lrO2C))*2 + ...
                  Ox(:,pindx.lrO2C).*d2LdO2.*Ox(:,pindx.lrO2C) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end 
        %rO2C O2P_T
        if (par.opt_rO2C == on & par.opt_O2P_T == on)
            kk = kk + 1 ;
            tmp = kC*DOC.*(O2C_rO2C*dRdO.*Ox(:,pindx.O2P_T)) + ...
                  Ox(:,pindx.lrO2C).*d2LdO2.*Ox(:,pindx.O2P_T) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end 
        % rO2C rO2P
        if (par.opt_rO2C == on & par.opt_rO2P == on)
            kk = kk + 1 ;
            tmp = kC*DOC.*(O2C_rO2C*dRdO.*Ox(:,pindx.lrO2P)) + ...
                  Ox(:,pindx.lrO2C).*d2LdO2.*Ox(:,pindx.lrO2P) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end
        % O2P_T O2P_T
        if (par.opt_O2P_T == on)
            kk = kk + 1 ;
            tmp = Ox(:,pindx.O2P_T).*d2LdO2.*Ox(:,pindx.O2P_T) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end
        % O2P_T rO2P
        if (par.opt_O2P_T & par.opt_rO2P)
            kk = kk + 1 ;
            tmp = Ox(:,pindx.O2P_T).*d2LdO2.*Ox(:,pindx.lrO2P) ;

            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end
        % rO2P rO2P
        if (par.opt_rO2P == on) 
            kk = kk + 1 ;
            tmp = -rO2P*diag(G) + ...
                  Ox(:,pindx.lrO2P).*d2LdO2.*Ox(:,pindx.lrO2P) ;
            
            Oxx(:,kk) = mfactor(FD, -tmp) ;
        end
    end

end

