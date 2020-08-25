function [par, O2, Ox, Oxx] = eqOcycle(x, par)
    on = true; off = false;
    global GO
    pindx = par.pindx;
    % unpack the parameters to be optimized
    % slopeo
    nox = 0; % count Omodel parameters
    if (par.opt_slopeo == on)
        nox = nox + 1;
        par.slopeo = x(pindx.slopeo);
    else
        par.slopeo  = par.slopeo;
    end
    
    % interpo
    if (par.opt_interpo == on)
        nox = nox + 1;
        linterpo = x(pindx.linterpo);
        par.interpo  = exp(linterpo);
    else
        par.interpo  = par.interpo;
    end
    par.nox = nox;
    %
    X0  = GO;
    options.iprint = 1;
    options.atol = 1e-10; options.rtol = 1e-10 ;
    [O2,ierr] = nsnew(X0,@(X) O_eqn(X, x, par),options) ;

    if (ierr ~=0)
        fprintf('o2model did not converge.\n') ;
        keyboard
    else
        % reset the global variable for the next call eqOcycle
        GO = real(O2) + 1e-5*randn(par.nwet,1);
        X0 = GO;
        F = O_eqn(O2, x, par);
        if norm(F) > 1e-12
            [O2,ierr] = nsnew(X0,@(X) O_eqn(X, x, par),options);
        end
        % Compute the gradient of the solution wrt the parameters
        [F, FD, Ox, Oxx] = O_eqn(O2, x, par);
    end

function [F, FD, Ox, Oxx] = O_eqn(O2, x, par)
    on = true; off = false;
    pindx = par.pindx;
    %
    % fixed parameters
    iwet = par.iwet;
    nwet = par.nwet;
    TRdiv = par.TRdiv;
    I = speye(nwet);
    PO4 = par.po4obs(iwet);
    % variables from C model
    DOC = par.DOC;
    %
    % tunable parameters;
    kappa_dc = par.kappa_dc;
    slopeo   = par.slopeo;
    interpo  = par.interpo;
    %
    vout = mkO2P(par);
    O2P  = vout.O2P;
    dO2Pdslopeo  = vout.dO2Pdslopeo;
    dO2Pdinterpo = vout.dO2Pdinterpo;
    %
    % O2 saturation concentration
    vout = Fsea2air(par,'O2');
    KO2 = vout.KO2;
    o2sat = vout.o2sat;
    % rate of o2 production
    G = par.G;
    PO2 = G*O2P;

    % parobolic function for o2 consumption
    R = 0.5 + 0.5*tanh(O2-10);
    dRdO = 0.5 - 0.5*tanh(O2-10).^2;
    d2RdO2 = -2*d0(dRdO)*tanh(O2-10);
    O2C = 170/117;
    % rate of o2 utilization
    LO2 = kappa_dc*(DOC.*O2C.*R);
    dLdO = kappa_dc*d0(DOC.*O2C.*dRdO);
    d2LdO2 = kappa_dc*DOC.*O2C.*d2RdO2;
    % O2 function
    F = TRdiv*O2 - PO2 + LO2 - KO2*(o2sat-O2);
    %
    if (nargout > 1)
        FD = mfactor(TRdiv + dLdO + KO2);
    end
    %
    %% ----------------- Gradient -------------------
    if (par.optim == off)
        Ox = [];
    elseif (par.optim & nargout > 2)
        npx = par.npx;
        ncx = par.ncx;
        nox = par.nox;
        nx = npx + ncx + nox;
        Ox = sparse(nwet, nx);
        Gx = par.Gx;
        DICx = par.Cx(1:nwet,:);
        DOCx = par.Cx(2*nwet+1:3*nwet,:);
        % P model only parameters
        for jj = 1:npx
            tmp = -d0(Gx(:,jj))*O2P + ...
                  kappa_dc*O2C.*R.*DOCx(:,jj);
            Ox(:,jj) = mfactor(FD, -tmp);
        end 

        % C model only parameters
        for jj = (npx+1) : (npx+ncx)
            if jj  ~= pindx.lkappa_dc;
                tmp = kappa_dc*O2C.*R.*DOCx(:,jj);
                Ox(:,jj) = mfactor(FD, -tmp);
            else
                tmp = kappa_dc*O2C.*R.*DOC + ...
                      kappa_dc*O2C.*R.*DOCx(:,jj);
                Ox(:,jj) = mfactor(FD, -tmp);
            end 
        end

        % O model only parameters
        if (par.opt_slopeo);
            tmp = -G*dO2Pdslopeo; 
            Ox(:,pindx.slopeo) = mfactor(FD, -tmp);
        end

        if (par.opt_interpo);
            tmp = -interpo*diag(G); 
            Ox(:,pindx.linterpo) = mfactor(FD, -tmp);
        end
    end

    %% ------------------- Hessian -------------------------
    if (par.optim == off)
        Oxx = [];
    elseif (par.optim & nargout > 3)
        %
        Gxx = par.Gxx;
        DICxx = par.Cxx(1:nwet,:);
        DOCxx = par.Cxx(2*nwet+1:3*nwet,:);
        Oxx = sparse(nwet,nchoosek(nx,2)+nx);
        % P model parameters
        kk = 0;
        for ju = 1 : npx
            for jo = ju : npx
                kk = kk + 1;
                tmp = -d0(Gxx(:,kk))*O2P + ...
                      kappa_dc*O2C.*R.*DOCxx(:,kk) + ...
                      kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) + ...
                      kappa_dc*O2C.*DOCx(:,jo).*dRdO.*Ox(:,ju) + ...
                      Ox(:,ju).*d2LdO2.*Ox(:,jo);
                Oxx(:,kk) = mfactor(FD, -tmp);                
            end 
        end
        %
        % P C model parameters
        for ju = 1 : npx
            for jo = (npx+1) : (npx+ncx)
                kk = kk + 1;
                % check if hessian associated with kappa_dc, if
                % yes, use formula (1), otherwise use formula (2);
                if (jo == pindx.lkappa_dc)
                    tmp = kappa_dc*O2C.*R.*DOCx(:,ju) + ...
                          kappa_dc*O2C.*DOC.*dRdO.*Ox(:,ju) + ...
                          kappa_dc*O2C.*R.*DOCxx(:,kk) + ...
                          kappa_dc*O2C.*DOCx(:,jo).*dRdO.*Ox(:,ju) + ...
                          kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo); %(1)
                    Oxx(:,kk) = mfactor(FD, -tmp);
                else
                    tmp = kappa_dc*O2C.*R.*DOCxx(:,kk) + ...
                          kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) + ...
                          kappa_dc*O2C.*DOCx(:,jo).*dRdO.*Ox(:,ju) + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo); %(2)
                    Oxx(:,kk) = mfactor(FD, -tmp);
                end
            end
        end
        
        % C model only parameters
        for ju = (npx+1) : (npx+ncx)
            for jo = ju : (npx+ncx)
                kk = kk + 1;
                % parameter pairs without kappa_dc
                if (ju ~= pindx.lkappa_dc & jo ~= pindx.lkappa_dc)
                    tmp = kappa_dc*O2C.*R.*DOCxx(:,kk) + ...
                          kappa_dc*O2C.*DOCx(:,jo).*dRdO.*Ox(:,ju) + ...
                          kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo);

                    Oxx(:,kk) = mfactor(FD, -tmp);
                    % parameter pairs with one kappa_dc
                elseif (ju == pindx.lkappa_dc & jo == pindx.lkappa_dc)
                    % kappa_dc kappa_dc
                    tmp = kappa_dc*O2C.*DOC.*R + ...
                          kappa_dc*O2C.*R.*DOCx(:,ju)*2 + ...
                          kappa_dc*O2C.*DOC.*dRdO.*Ox(:,ju)*2 + ... %Fxx
                          kappa_dc*O2C.*R.*DOCxx(:,kk) + ... %dFdCxx
                          kappa_dc*O2C.*DOCx(:,jo).*dRdO.*Ox(:,ju) + ... %d2FdCdO*Cx*Ox                  
                          kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo); %d2FdO*2Oxx

                    Oxx(:,kk) = mfactor(FD, -tmp);
                else %parameter pairs with only one kappa_dc
                    if ju < pindx.lkappa_dc
                        jk = ju;
                    else
                        jk = jo;
                    end
                    %
                    tmp = kappa_dc*O2C.*R.*DOCx(:,jk) + ...
                          kappa_dc*O2C.*DOC.*dRdO.*Ox(:,jk) + ... %Fxx
                          kappa_dc*O2C.*R.*DOCxx(:,kk) + ... %dFdCxx
                          kappa_dc*O2C.*DOCx(:,jo).*dRdO.*Ox(:,ju) + ... %d2FdCdO*Cx*Ox
                          kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) + ...
                          Ox(:,ju).*d2LdO2.*Ox(:,jo); %d2FdO*2Oxx

                    Oxx(:,kk) = mfactor(FD, -tmp);
                end
            end 
        end 
        %
        % P and O model only parameters
        % P model parameters
        dO2Pdp{1} = dO2Pdslopeo;
        dO2Pdp{2} = interpo*dO2Pdinterpo;
        for ju = 1 : npx
            for jo = (npx+ncx+1) : (npx+ncx+nox)
                kk = kk + 1;
                tmp = -d0(Gx(:,ju))*dO2Pdp{jo-npx-ncx} + ...
                      kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) + ...
                      Ox(:,ju).*d2LdO2.*Ox(:,jo);

                Oxx(:,kk) = mfactor(FD, -tmp);                
            end 
        end
        
        % C and O model only parameters
        for ju = (npx+1) : (npx+ncx)
            for jo = (npx+ncx+1) : (npx+ncx+nox)
                kk = kk + 1;
                if (ju ~= pindx.lkappa_dc)
                    tmp = kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) ...
                          + Ox(:,ju).*d2LdO2.*Ox(:,jo);

                    Oxx(:,kk) = mfactor(FD, -tmp);
                else
                    tmp = kappa_dc*O2C.*DOC.*dRdO.*Ox(:,jo) + ...
                          kappa_dc*O2C.*DOCx(:,ju).*dRdO.*Ox(:,jo) ...
                          + Ox(:,ju).*d2LdO2.*Ox(:,jo);

                    Oxx(:,kk) = mfactor(FD, -tmp);
                end 
            end 
        end 
        
        % O model only parameters
        % slopeo slopeo
        if (par.opt_slopeo == on)
            kk = kk + 1;
            tmp = Ox(:,pindx.slopeo).*d2LdO2.*Ox(:,pindx.slopeo);

            Oxx(:,kk) = mfactor(FD, -tmp);
        end

        % slopeo interpo
        if (par.opt_slopeo & par.opt_interpo)
            kk = kk + 1;
            tmp = Ox(:,pindx.slopeo).*d2LdO2.*Ox(:,pindx.linterpo);

            Oxx(:,kk) = mfactor(FD, -tmp);
        end
        
        % interpo interpo
        if (par.opt_interpo == on);
            kk = kk + 1;
            tmp = -interpo*diag(G) + ...
                  Ox(:,pindx.linterpo).*d2LdO2.*Ox(:,pindx.linterpo);
            
            Oxx(:,kk) = mfactor(FD, -tmp);
        end
    end
