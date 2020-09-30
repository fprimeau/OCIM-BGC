function [par, C, Cx, Cxx] = eqCcycle(x, par)
% ip is the mapping from x to parameter names (see switch below)
% output: C is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
    on = true; off = false;
    global GC
    iwet = par.iwet;
    nwet = par.nwet;

    % unpack the parameters to be optimized

    % bC_T
    if (par.opt_bC_T == on)
        par.bC_T = x(par.pindx.bC_T);
    end

    % bC
    if (par.opt_bC == on)
        lbC = x(par.pindx.lbC) ;
        par.bC  = exp(lbC)     ;
    end

    % kPIC
    if (par.opt_kPIC == on)
        lkPIC = x(par.pindx.lkPIC) ;
        par.kPIC = exp(lkPIC)      ;
    end

    % d
    if (par.opt_d == on)
        ld = x(par.pindx.ld) ;
        par.d = exp(ld)      ;
    end

    % kC_T
    if (par.opt_kC_T == on)
        par.kC_T = x(par.pindx.kC_T);
    end

    % kdC
    if (par.opt_kdC == on)
        lkdC = x(par.pindx.lkdC);
        par.kdC  = exp(lkdC);
    end

    % R_Si
    if (par.opt_R_Si == on)
        lR_Si = x(par.pindx.lR_Si);
        par.R_Si = exp(lR_Si) ;
    end

    % rR
    if (par.opt_rR == on)
        lrR = x(par.pindx.lrR);
        par.rR = exp(lrR);
    end

    % cc
    if (par.opt_cc)
        lcc = x(par.pindx.lcc);
        par.cc = exp(lcc);
    end 

    % dd
    if (par.opt_dd)
        ldd = x(par.pindx.ldd);
        par.dd = exp(ldd);
    end
    %
    options.iprint = 0   ; 
    options.atol = 1e-10 ;
    options.rtol = 1e-10 ;
    fprintf('Solving C model ...\n') ;

    X0  = GC;
    [C,ierr] = nsnew(X0,@(X) C_eqn(X, par),options) ;
    if (ierr ~= 0)
        fprintf('eqCcycle did not converge.\n') ;
        npx  = par.npx   ;
        ncx  = par.ncx   ;
        nx   = npx + ncx ;
        Cx   = sparse(4*par.nwet, nx) ;
        Cxx  = sparse(4*par.nwet, nchoosek(nx,2)+nx) ;
        [par.G,par.Gx,par.Gxx] = uptake_C(par) ;
    else
        % reset the global variable for the next call eqCcycle
        GC = real(C) + 1e-6*randn(4*nwet,1) ;
        X0 = GC;
        F = C_eqn(C, par) ;
        % test if norm of F small enough, if now rerun nsnew;
        if norm(F) > 1e-12
            [C,ierr] = nsnew(X0,@(X) C_eqn(X, par),options);
        end 
        %
        if nargout > 2
            [F,FD,Cx,Cxx,par] = C_eqn(C, par);
        end 
    end
end

function [F,FD,Cx,Cxx,par] = C_eqn(X, par)    
% unpack some useful stuff
    on = true; off = false;
    grd   = par.grd   ;
    M3d   = par.M3d   ;
    TRdiv = par.TRdiv ;
    iwet  = par.iwet  ;
    nwet  = par.nwet  ;
    I     = par.I     ;

    Tz  = par.Tz ;
    DIC = X(0*nwet+1:1*nwet) ; 
    POC = X(1*nwet+1:2*nwet) ;
    DOC = X(2*nwet+1:3*nwet) ;
    PIC = X(3*nwet+1:4*nwet) ;
    %
    PO4 = par.po4obs(iwet);
    %
    % fixed parameters
    kappa_p = par.kappa_p ;
    % parameters need to be optimized
    alpha = par.alpha ;
    beta  = par.beta  ;
    sigma = par.sigma ;
    bC_T  = par.bC_T  ;
    bC    = par.bC    ;
    kPIC  = par.kPIC  ;
    d     = par.d     ;
    kC_T  = par.kC_T  ;
    kdC   = par.kdC   ;
    R_Si  = par.R_Si  ;
    rR    = par.rR    ;
    cc    = par.cc    ;
    dd    = par.dd    ;
    
    % PIC to production ratio 
    vout  = mkPIC2P(par) ;
    RR    = vout.RR    ;
    RR_Si = vout.RR_Si ;
    RR_rR = vout.RR_rR ;
    clear vout 
    % kappa_dc ;
    kC    = d0(kC_T * Tz + kdC) ;
    C2P   = 1./(cc*PO4 + dd) ;
    par.C2P = C2P ;
    % particle flux div_rergence [s^-1];
    PFDa = buildPFD(par,'PIC') ;
    PFDc = buildPFD(par,'POC') ;
    par.PFDa = PFDa ;
    par.PFDc = PFDc ;
    par.DIC  = DIC  ;
    % Air-Sea gas exchange
    vout  = Fsea2air(par, 'CO2');
    KG    = vout.KG;
    KGG   = vout.KGG;
    JgDIC = vout.JgDIC;
    clear vout 
    % biological DIC uptake operator
    G = uptake_C(par); par.G = G;
    
    eq1 = (I+(1-sigma)*RR)*(G*C2P) + TRdiv*DIC - kC*DOC - kPIC*PIC - JgDIC;
    eq2 = -(1-sigma)*G*C2P + (PFDc+kappa_p*I)*POC      ;
    eq3 = -sigma*G*C2P + (TRdiv+kC)*DOC - kappa_p*POC  ;
    eq4 = -(1-sigma)*RR*(G*C2P) + (PFDa+kPIC*I)*PIC ; 
    
    F   = [eq1; eq2; eq3; eq4];

    if nargout > 1
        % construct the LHS matrix for the offline model
        % disp('Preparing LHS and RHS matrix:')
        
        % colum 1 dFdDIC
        Jc{1,1} = TRdiv - KG ;
        Jc{2,1} = 0*I ;
        Jc{3,1} = 0*I ;
        Jc{4,1} = 0*I ;
        
        % colum 2 dFdPOC
        Jc{1,2} = 0*I ;
        Jc{2,2} = PFDc + kappa_p*I ;
        Jc{3,2} = -kappa_p*I ;
        Jc{4,2} = 0*I ;
        
        % colum 3 dFdDOC
        Jc{1,3} = -kC ;
        Jc{2,3} = 0*I ;
        Jc{3,3} = TRdiv + kC ;
        Jc{4,3} = 0*I ;
        
        % colum 4 dFdPIC
        Jc{1,4} = -kPIC*I ;
        Jc{2,4} = 0*I ;
        Jc{3,4} = 0*I ;
        Jc{4,4} = PFDa + kPIC*I ;
        
        % factorize Jacobian matrix
        FD = mfactor(cell2mat(Jc)) ;
    end 

    if (par.optim == off)
        Cx = [];
    elseif (par.optim & nargout > 2)
        pindx = par.pindx;
        Z = sparse(nwet,1);
        C2P_cc = -PO4./(cc*PO4 + dd).^2;
        C2P_dd = -1./(cc*PO4 + dd).^2;
        par.C2P_cc = C2P_cc;
        par.C2P_dd = C2P_dd;
        [~,Gx] = uptake_C(par);
        par.Gx = Gx;
        npx    = par.npx;
        % P model parameters
        % sigma
        if (par.opt_sigma == on)
            tmp = [sigma*RR*G*C2P; ...
                   -sigma*G*C2P; ...
                   sigma*G*C2P; ...
                   -sigma*RR*G*C2P] + ...
                  [-d0((I+(1-sigma)*RR)*Gx(:,pindx.lsigma))*C2P; ...
                   (1-sigma)*d0(Gx(:,pindx.lsigma))*C2P; ...
                   sigma*d0(Gx(:,pindx.lsigma))*C2P; ...
                   d0((1-sigma)*RR*Gx(:,pindx.lsigma))*C2P];
            
            Cx(:,pindx.lsigma) = mfactor(FD, tmp);
        end

        % kP_T
        if (par.opt_kP_T == on)
            tmp = [-d0((I+(1-sigma)*RR)*Gx(:,pindx.kP_T))*C2P; ...
                   (1-sigma)*d0(Gx(:,pindx.kP_T))*C2P; ...
                   sigma*d0(Gx(:,pindx.kP_T))*C2P; ...
                   d0((1-sigma)*RR*Gx(:,pindx.kP_T))*C2P];
            
            Cx(:,pindx.kP_T) = mfactor(FD, tmp);
        end
        
        % kdP
        if (par.opt_kdP == on)
            tmp = [-d0((I+(1-sigma)*RR)*Gx(:,pindx.lkdP))*C2P; ...
                   (1-sigma)*d0(Gx(:,pindx.lkdP))*C2P; ...
                   sigma*d0(Gx(:,pindx.lkdP))*C2P; ...
                   d0((1-sigma)*RR*Gx(:,pindx.lkdP))*C2P];
            
            Cx(:,pindx.lkdP) = mfactor(FD, tmp);
        end
        
        % bP_T
        if (par.opt_bP_T == on)
            tmp = [-d0((I+(1-sigma)*RR)*Gx(:,pindx.bP_T))*C2P; ...
                   (1-sigma)*d0(Gx(:,pindx.bP_T))*C2P; ...
                   sigma*d0(Gx(:,pindx.bP_T))*C2P; ...
                   d0((1-sigma)*RR*Gx(:,pindx.bP_T))*C2P];
            
            Cx(:,pindx.bP_T) = mfactor(FD, tmp);
        end
        
        % bP
        if (par.opt_bP == on)
            tmp = [-d0((I+(1-sigma)*RR)*Gx(:,pindx.lbP))*C2P;...
                   (1-sigma)*d0(Gx(:,pindx.lbP))*C2P;...
                   sigma*d0(Gx(:,pindx.lbP))*C2P;...
                   d0((1-sigma)*RR*Gx(:,pindx.lbP))*C2P];
            
            Cx(:,pindx.lbP) = mfactor(FD, tmp);
        end
        
        % alpha
        if (par.opt_alpha == on)
            tmp = [-d0((I+(1-sigma)*RR)*Gx(:,pindx.lalpha))*C2P; ...
                   (1-sigma)*d0(Gx(:,pindx.lalpha))*C2P; ...
                   sigma*d0(Gx(:,pindx.lalpha))*C2P; ...
                   d0((1-sigma)*RR*Gx(:,pindx.lalpha))*C2P];
            
            Cx(:,pindx.lalpha) = mfactor(FD, tmp);
        end
        
        % beta
        if (par.opt_beta == on)
            tmp = [-d0((I+(1-sigma)*RR)*Gx(:,pindx.lbeta))*C2P;...
                   (1-sigma)*d0(Gx(:,pindx.lbeta))*C2P;...
                   sigma*d0(Gx(:,pindx.lbeta))*C2P;...
                   d0((1-sigma)*RR*Gx(:,pindx.lbeta))*C2P];
            
            Cx(:,pindx.lbeta) = mfactor(FD, tmp);
        end
        % -------------------- C parameters ------------------
        % bC_T
        if (par.opt_bC_T == on)
            [~,Gout]   = buildPFD(par,'POC');
            PFD_bm     = Gout.PFD_bm;
            par.PFD_bm = PFD_bm;
            tmp = [Z; -PFD_bm*POC; Z;  Z];
            
            Cx(:,pindx.bC_T) = mfactor(FD, tmp);
        end
        
        % bC
        if (par.opt_bC == on)
            [~,Gout]   = buildPFD(par,'POC');
            PFD_bb     = Gout.PFD_bb;
            par.PFD_bb = PFD_bb;
            tmp = bC*[Z; -PFD_bb*POC; Z; Z];
            
            Cx(:,pindx.lbC) = mfactor(FD, tmp);
        end

        % kPIC
        if (par.opt_kPIC == on)
            [~,Gout]  = buildPFD(par,'PIC') ;
            PFD_k     = Gout.PFD_k ;
            par.PFD_k = PFD_k      ;
            % tmp = [Z; ...
                   % Z; ...
                   % Z; ...
                   % -kPIC*PFD_k*PIC] ;

            tmp = [kPIC*PIC; ...
                   Z; ...
                   Z; ...
                   -kPIC*PIC];
            
            Cx(:,pindx.lkPIC) = mfactor(FD, tmp) ;
        end
        
        % d
        if (par.opt_d == on)
            [~,Gout]  = buildPFD(par,'PIC');
            PFD_d     = Gout.PFD_d;
            par.PFD_d = PFD_d;
            tmp = d*[Z; Z; Z; -PFD_d*PIC];
            
            Cx(:,pindx.ld) = mfactor(FD, tmp);
        end

        % kC_T
        if (par.opt_kC_T == on)
            kC_kC_T     = d0(Tz) ; 
            par.kC_kC_T = kC_kC_T ;
            tmp = [kC_kC_T*DOC; Z; -kC_kC_T*DOC; Z];
            
            Cx(:,pindx.kC_T) = mfactor(FD, tmp);
        end
        
        % kdC
        if (par.opt_kdC == on)
            kC_kdC     = kdC ;
            par.kC_kdC = kC_kdC;
            tmp = kC_kdC*[DOC; Z; -DOC; Z];
            
            Cx(:,pindx.lkdC) = mfactor(FD, tmp);
        end

        % R_Si
        if (par.opt_R_Si == on)
            tmp = [-(1-sigma)*RR_Si*(G*C2P); ...
                   Z; ...
                   Z; ...
                   (1-sigma)*RR_Si*(G*C2P)];
            
            Cx(:,pindx.lR_Si) = mfactor(FD, tmp);
        end
        
        % rR
        if (par.opt_rR == on)
            tmp = [-(1-sigma)*RR_rR*(G*C2P); ...
                   Z; ...
                   Z; ...
                   (1-sigma)*RR_rR*(G*C2P)];
            
            Cx(:,pindx.lrR) = mfactor(FD, tmp);
        end
        
        % cc
        if (par.opt_cc == on)
            tmp = cc*[-(I+(1-sigma)*RR)*(G*C2P_cc); ...
                      (1-sigma)*G*C2P_cc; ...
                      sigma*G*C2P_cc; ...
                      (1-sigma)*RR*(G*C2P_cc)];
            
            Cx(:,pindx.lcc) = mfactor(FD, tmp);
        end
        
        % dd
        if (par.opt_dd == on)
            tmp = dd*[-(I+(1-sigma)*RR)*(G*C2P_dd); ...
                      (1-sigma)*G*C2P_dd; ...
                      sigma*G*C2P_dd; ...
                      (1-sigma)*RR*(G*C2P_dd)];
            
            Cx(:,pindx.ldd) = mfactor(FD, tmp);
        end
    end

    if (par.optim == off)
        Cxx = [];
    elseif (par.optim & nargout > 3);
        p2c       = cc*PO4 + dd;
        C2P_dd_dd = 2./p2c.^3;
        C2P_cc_cc = (2*PO4.^2)./p2c.^3;
        C2P_cc_dd = (2*PO4)./p2c.^3;
        [~,~,Gxx] = uptake_C(par);
        par.Gxx   = Gxx;
        DICx = Cx(0*nwet+1:1*nwet,:);
        POCx = Cx(1*nwet+1:2*nwet,:);
        DOCx = Cx(2*nwet+1:3*nwet,:);
        PICx = Cx(3*nwet+1:end,:);
        % ------------------------------------------------------
        % P model only parameters
        kk = 0;
        for jj = 1:npx
            for jk = jj:npx
                kk = kk + 1;
                if (par.opt_sigma == on)
                    % sigma sigma
                    if (jj == jk & jj == pindx.lsigma)
                        tmp = sigma*[RR*G*C2P + 2*d0(RR*Gx(:,jj))*C2P; ...
                                     -G*C2P - 2*d0(Gx(:,jj))*C2P;...
                                      G*C2P + 2*d0(Gx(:,jj))*C2P;...
                                     -RR*G*C2P - 2*d0(RR*Gx(:,jj))*C2P] + ...
                              [-d0((I+(1-sigma)*RR)*Gxx(:,kk))*C2P; ...
                               (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                               sigma*d0(Gxx(:,kk))*C2P; ...
                               d0((1-sigma)*RR*Gxx(:,kk))*C2P];
                        
                        Cxx(:,kk) = mfactor(FD, tmp);
                        %pairs not assciated with sigma;
                    elseif (jj ~= pindx.lsigma & jk ~= pindx.lsigma)
                        tmp = [-d0((I+(1-sigma)*RR)*Gxx(:,kk))*C2P; ...
                               (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                               sigma*d0(Gxx(:,kk))*C2P; ...
                               d0((1-sigma)*RR*Gxx(:,kk))*C2P];
                        
                        Cxx(:,kk) = mfactor(FD, tmp);
                    else 
                        tmp = [RR*sigma*d0(Gx(:,jk))*C2P; ...
                               -sigma*d0(Gx(:,jk))*C2P; ...
                               sigma*d0(Gx(:,jk))*C2P; ...
                               -RR*sigma*d0(Gx(:,jk))*C2P] + ...
                              [-d0((I+(1-sigma)*RR)*Gxx(:,kk))*C2P; ...
                               (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                               sigma*d0(Gxx(:,kk))*C2P; ...
                               d0((1-sigma)*RR*Gxx(:,kk))*C2P];
                        
                        Cxx(:,kk) = mfactor(FD, tmp);
                    end
                else 
                    tmp = [-d0((I+(1-sigma)*RR)*Gxx(:,kk))*C2P; ...
                           (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                           sigma*d0(Gxx(:,kk))*C2P; ...
                           d0((1-sigma)*RR*Gxx(:,kk))*C2P];
                    
                    Cxx(:,kk) = mfactor(FD, tmp);
                    % sigma foo
                end 
            end
        end
        % ------------------------------------------------------
        % P and  C model parameters
        % sigma bC_T
        if (par.opt_sigma & par.opt_bC_T)
            kk  = kk + 1;
            tmp =  [Z; ...
                    -PFD_bm*POCx(:,pindx.lsigma) ; ...
                    Z; ...
                    Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % sigma bC
        if (par.opt_sigma & par.opt_bC)
            kk = kk + 1;
            tmp =  [Z ; ...
                    -bC*PFD_bb*POCx(:,pindx.lsigma); ...
                    Z ; ...
                    Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % sigma kPIC
        if (par.opt_sigma & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.lsigma) ; ...
                    Z ;
                    Z ; ...
                    [-kPIC*PICx(:,pindx.lsigma)]];
            % -kPIC*PFD_k*PICx(:,pindx.lsigma) + ...
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % sigma d
        if (par.opt_sigma & par.opt_d)
            kk = kk + 1;
            tmp = [Z ; ...
                   Z ; ...
                   Z ; ...
                   -d*PFD_d*PICx(:,pindx.lsigma)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % sigma kC_T
        if (par.opt_sigma & par.opt_kC_T)
            kk  = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lsigma); ...
                   Z; ...
                   -kC_kC_T*DOCx(:,pindx.lsigma); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % sigma kdC
        if (par.opt_sigma & par.opt_kdC)
            kk  = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.lsigma); ...
                   Z; ...
                   -kC_kdC*DOCx(:,pindx.lsigma); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % sigma R_Si
        if (par.opt_sigma & par.opt_R_Si)
            kk = kk + 1;
            tmp = [sigma*RR_Si*G*C2P ; ...
                   Z ; ...
                   Z ; ...
                   -sigma*RR_Si*G*C2P] + ... 
                  [-d0((1-sigma)*RR_Si*Gx(:,pindx.lsigma))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_Si*Gx(:,pindx.lsigma))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % sigma rR
        if (par.opt_sigma & par.opt_rR)
            kk = kk + 1;
            tmp = [sigma*RR_rR*G*C2P ; ...
                   Z ; ...
                   Z ; ...
                   -sigma*RR_rR*G*C2P] + ... 
                  [-d0((1-sigma)*RR_rR*Gx(:,pindx.lsigma))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_rR*Gx(:,pindx.lsigma))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % sigma cc
        if (par.opt_sigma & par.opt_cc)
            kk  = kk + 1;
            tmp = cc*[RR*sigma*G*C2P_cc; ...
                      -sigma*G*C2P_cc; ...
                      sigma*G*C2P_cc; ...
                      -RR*sigma*G*C2P_cc] + ...
                  cc*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lsigma))*C2P_cc; ...
                      (1-sigma)*d0(Gx(:,pindx.lsigma))*C2P_cc; ... 
                      sigma*d0(Gx(:,pindx.lsigma))*C2P_cc; ... 
                      d0((1-sigma)*RR*Gx(:,pindx.lsigma))*C2P_cc];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % sigma dd
        if (par.opt_sigma & par.opt_dd)
            kk  = kk + 1;
            tmp = dd*[RR*sigma*G*C2P_dd; ...
                      -sigma*G*C2P_dd; ...
                      sigma*G*C2P_dd; ...
                      -RR*sigma*G*C2P_dd] + ...
                  dd*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lsigma))*C2P_dd; ...
                      (1-sigma)*d0(Gx(:,pindx.lsigma))*C2P_dd; ...
                      sigma*d0(Gx(:,pindx.lsigma))*C2P_dd; ... 
                      d0((1-sigma)*RR*Gx(:,pindx.lsigma))*C2P_dd];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kP_T bC_T
        if (par.opt_kP_T & par.opt_bC_T)
            kk  = kk + 1;
            tmp = [Z; ...
                   -PFD_bm*POCx(:,pindx.kP_T); ...
                   Z; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kP_T bC
        if (par.opt_kP_T & par.opt_bC)
            kk = kk + 1;
            tmp = [Z ; ...
                   -bC*PFD_bb*POCx(:,pindx.kP_T); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kP_T kPIC
        if (par.opt_kP_T & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.kP_T) ; ...
                    Z ;
                    Z ; ...
                    -kPIC*PICx(:,pindx.kP_T)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kP_T d
        if (par.opt_kP_T & par.opt_d)
            kk = kk + 1;
            tmp = [Z; ...
                   Z; ...
                   Z; ...
                   -d*PFD_d*PICx(:,pindx.kP_T)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kP_T kC_T
        if (par.opt_kP_T & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.kP_T); ...
                   Z ; ...
                   -kC_kC_T*DOCx(:,pindx.kP_T); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kP_T kdC
        if (par.opt_kP_T & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.kP_T); ...
                   Z ; ...
                   -kC_kdC*DOCx(:,pindx.kP_T); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kP_T R_Si
        if (par.opt_kP_T & par.opt_R_Si)
            kk = kk + 1;
            tmp = [-d0((1-sigma)*RR_Si*Gx(:,pindx.kP_T))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_Si*Gx(:,pindx.kP_T))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kP_T RR_rR
        if (par.opt_kP_T & par.opt_rR)
            kk = kk + 1;
            tmp = [-d0((1-sigma)*RR_rR*Gx(:,pindx.kP_T))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_rR*Gx(:,pindx.kP_T))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kP_T cc
        if (par.opt_kP_T & par.opt_cc)
            kk = kk + 1;
            tmp = cc*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.kP_T))*C2P_cc; ...
                      (1-sigma)*d0(Gx(:,pindx.kP_T))*C2P_cc; ...
                      sigma*d0(Gx(:,pindx.kP_T))*C2P_cc; ...
                      d0((1-sigma)*RR*Gx(:,pindx.kP_T))*C2P_cc];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kP_T dd
        if (par.opt_kP_T & par.opt_dd)
            kk = kk + 1;
            tmp = dd*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.kP_T))*C2P_dd; ...
                      (1-sigma)*d0(Gx(:,pindx.kP_T))*C2P_dd; ...
                      sigma*d0(Gx(:,pindx.kP_T))*C2P_dd; ...
                      d0((1-sigma)*RR*Gx(:,pindx.kP_T))*C2P_dd];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdP bC_T
        if (par.opt_kdP & par.opt_bC_T)
            kk = kk + 1;
            tmp =  [Z; ...
                    -PFD_bm*POCx(:,pindx.lkdP); ...
                    Z; ...
                    Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdP bC
        if (par.opt_kdP & par.opt_bC)
            kk = kk + 1;
            tmp = [Z ; ...
                   -bC*PFD_bb*POCx(:,pindx.lkdP); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kdP kPIC
        if (par.opt_kdP & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.lkdP) ; ...
                    Z ;
                    Z ; ...
                    -kPIC*PICx(:,pindx.lkdP)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdP d
        if (par.opt_kdP & par.opt_d)
            kk = kk + 1;
            tmp = [Z; ...
                   Z; ...
                   Z; ...
                   -d*PFD_d*PICx(:,pindx.lkdP)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kdP kC_T
        if (par.opt_kdP & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lkdP); ...
                   Z ; ...
                   -kC_kC_T*DOCx(:,pindx.lkdP); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdP kdC
        if (par.opt_kdP & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.lkdP); ...
                   Z ; ...
                   -kC_kdC*DOCx(:,pindx.lkdP); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdP R_Si
        if (par.opt_kdP & par.opt_R_Si)
            kk  = kk + 1;
            tmp = [-d0((1-sigma)*RR_Si*Gx(:,pindx.lkdP))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_Si*Gx(:,pindx.lkdP))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kdP rR
        if (par.opt_kdP & par.opt_rR)
            kk  = kk + 1;
            tmp = [-d0((1-sigma)*RR_rR*Gx(:,pindx.lkdP))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_rR*Gx(:,pindx.lkdP))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdP cc
        if (par.opt_kdP & par.opt_cc)
            kk  = kk + 1;
            tmp = cc*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lkdP))*C2P_cc; ...
                      (1-sigma)*d0(Gx(:,pindx.lkdP))*C2P_cc; ...
                      sigma*d0(Gx(:,pindx.lkdP))*C2P_cc; ...
                      d0((1-sigma)*RR*Gx(:,pindx.lkdP))*C2P_cc];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdP dd
        if (par.opt_kdP & par.opt_dd)
            kk  = kk + 1;
            tmp = dd*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lkdP))*C2P_dd; ...
                      (1-sigma)*d0(Gx(:,pindx.lkdP))*C2P_dd; ...
                      sigma*d0(Gx(:,pindx.lkdP))*C2P_dd; ...
                      d0((1-sigma)*RR*Gx(:,pindx.lkdP))*C2P_dd];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP_T bC_T
        if (par.opt_bP_T & par.opt_bC_T)
            kk  = kk + 1;
            tmp = [Z ; ...
                   -PFD_bm*POCx(:,pindx.bP_T); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP_T bC
        if (par.opt_bP_T & par.opt_bC)
            kk = kk + 1;
            tmp = [Z ; ...
                   -bC*PFD_bb*POCx(:,pindx.bP_T); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bP_T kPIC
        if (par.opt_bP_T & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.bP_T) ; ...
                    Z ;
                    Z ; ...
                    -kPIC*PICx(:,pindx.bP_T)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP_T d
        if (par.opt_bP_T & par.opt_d)
            kk = kk + 1;
            tmp = [Z ; ...
                   Z ; ...
                   Z ; ...
                   -d*PFD_d*PICx(:,pindx.bP_T)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bP_T kC_T
        if (par.opt_bP_T & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.bP_T); ...
                   Z ; ...
                   -kC_kC_T*DOCx(:,pindx.bP_T); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bP_T kdC
        if (par.opt_bP_T & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.bP_T); ...
                   Z; ...
                   -kC_kdC*DOCx(:,pindx.bP_T); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP_T R_Si
        if (par.opt_bP_T & par.opt_R_Si)
            kk = kk + 1;
            tmp = [-d0((1-sigma)*RR_Si*Gx(:,pindx.bP_T))*C2P ; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_Si*Gx(:,pindx.bP_T))*C2P] ;
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bP_T rR
        if (par.opt_bP_T & par.opt_rR)
            kk = kk + 1;
            tmp = [-d0((1-sigma)*RR_rR*Gx(:,pindx.bP_T))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_rR*Gx(:,pindx.bP_T))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP_T cc
        if (par.opt_bP_T & par.opt_cc)
            kk = kk + 1;
            tmp = cc*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.bP_T))*C2P_cc; ...
                      (1-sigma)*d0(Gx(:,pindx.bP_T))*C2P_cc; ...
                      sigma*d0(Gx(:,pindx.bP_T))*C2P_cc; ...
                      d0((1-sigma)*RR*Gx(:,pindx.bP_T))*C2P_cc];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP_T dd
        if (par.opt_bP_T & par.opt_dd)
            kk = kk + 1;
            tmp = dd*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.bP_T))*C2P_dd; ...
                      (1-sigma)*d0(Gx(:,pindx.bP_T))*C2P_dd; ...
                      sigma*d0(Gx(:,pindx.bP_T))*C2P_dd; ...
                      d0((1-sigma)*RR*Gx(:,pindx.bP_T))*C2P_dd];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP bC_T
        if (par.opt_bP & par.opt_bC_T)
            kk = kk + 1;
            tmp = [Z; ...
                   -PFD_bm*POCx(:,pindx.lbP); ...
                   Z; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP bC
        if (par.opt_bP & par.opt_bC)
            kk = kk + 1;
            tmp = [Z; ...
                   -bC*PFD_bb*POCx(:,pindx.lbP); ...
                   Z; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bP kPIC
        if (par.opt_bP & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.lbP) ; ...
                    Z ;
                    Z ; ...
                    -kPIC*PICx(:,pindx.lbP)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP d
        if (par.opt_bP & par.opt_d)
            kk = kk + 1;
            tmp = [Z; ...
                   Z; ...
                   Z; ...
                   -d*PFD_d*PICx(:,pindx.lbP)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bP kC_T
        if (par.opt_bP & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lbP); ...
                   Z ; ...
                   -kC_kC_T*DOCx(:,pindx.lbP); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bP kdC
        if (par.opt_bP & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.lbP); ...
                   Z ; ...
                   -kC_kdC*DOCx(:,pindx.lbP); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP R_Si
        if (par.opt_bP & par.opt_R_Si)
            kk = kk + 1;
            tmp = [-d0((1-sigma)*RR_Si*Gx(:,pindx.lbP))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_Si*Gx(:,pindx.lbP))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bP rR
        if (par.opt_bP & par.opt_rR)
            kk = kk + 1;
            tmp = [-d0((1-sigma)*RR_rR*Gx(:,pindx.lbP))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_rR*Gx(:,pindx.lbP))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP cc
        if (par.opt_bP & par.opt_cc)
            kk = kk + 1;
            tmp = cc*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lbP))*C2P_cc; ...
                      (1-sigma)*d0(Gx(:,pindx.lbP))*C2P_cc;...
                      sigma*d0(Gx(:,pindx.lbP))*C2P_cc;...
                      d0((1-sigma)*RR*Gx(:,pindx.lbP))*C2P_cc];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bP dd
        if (par.opt_bP & par.opt_dd)
            kk = kk + 1;
            tmp = dd*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lbP))*C2P_dd; ...
                      (1-sigma)*d0(Gx(:,pindx.lbP))*C2P_dd;...
                      sigma*d0(Gx(:,pindx.lbP))*C2P_dd;...
                      d0((1-sigma)*RR*Gx(:,pindx.lbP))*C2P_dd];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % alpha bC_T
        if (par.opt_alpha & par.opt_bC_T)
            kk = kk + 1;
            tmp = [Z ; ...
                   -PFD_bm*POCx(:,pindx.lalpha); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % alpha bC
        if (par.opt_alpha & par.opt_bC)
            kk = kk + 1;
            tmp = [Z ; ...
                   -bC*PFD_bb*POCx(:,pindx.lalpha); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % alpha kPIC
        if (par.opt_alpha & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.lalpha) ; ...
                    Z ;
                    Z ; ...
                    -kPIC*PICx(:,pindx.lalpha)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % alpha d
        if (par.opt_alpha & par.opt_d)
            kk = kk + 1;
            tmp = [Z ; ...
                   Z ; ...
                   Z ; ...
                   -d*PFD_d*PICx(:,pindx.lalpha)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % alpha kC_T
        if (par.opt_alpha & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lalpha); ...
                   Z ; ...
                   -kC_kC_T*DOCx(:,pindx.lalpha); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % alpha kdC
        if (par.opt_alpha & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.lalpha); ...
                   Z ; ...
                   -kC_kdC*DOCx(:,pindx.lalpha); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % alpha R_Si
        if (par.opt_alpha & par.opt_R_Si)
            kk  = kk + 1;
            tmp = [-d0((1-sigma)*RR_Si*Gx(:,pindx.lalpha))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_Si*Gx(:,pindx.lalpha))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % alpha rR
        if (par.opt_alpha & par.opt_rR)
            kk  = kk + 1;
            tmp = [-d0((1-sigma)*RR_rR*Gx(:,pindx.lalpha))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_rR*Gx(:,pindx.lalpha))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % alpha cc
        if (par.opt_alpha & par.opt_cc)
            kk  = kk + 1;
            tmp = cc*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lalpha))*C2P_cc; ...
                      (1-sigma)*d0(Gx(:,pindx.lalpha))*C2P_cc; ...
                      sigma*d0(Gx(:,pindx.lalpha))*C2P_cc; ...
                      d0((1-sigma)*RR*Gx(:,pindx.lalpha))*C2P_cc];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % alpha dd
        if (par.opt_alpha & par.opt_dd)
            kk = kk + 1;
            tmp = dd*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lalpha))*C2P_dd; ...
                      (1-sigma)*d0(Gx(:,pindx.lalpha))*C2P_dd; ...
                      sigma*d0(Gx(:,pindx.lalpha))*C2P_dd; ...
                      d0((1-sigma)*RR*Gx(:,pindx.lalpha))*C2P_dd];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % beta bC_T
        if (par.opt_beta & par.opt_bC_T)
            kk = kk + 1;
            tmp = [Z ; ...
                   -PFD_bm*POCx(:,pindx.lbeta); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % beta bC
        if (par.opt_beta & par.opt_bC)
            kk = kk + 1;
            tmp = [Z ; ...
                   -bC*PFD_bb*POCx(:,pindx.lbeta); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % beta kPIC
        if (par.opt_beta & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.lbeta) ; ...
                    Z ;
                    Z ; ...
                    -kPIC*PICx(:,pindx.lbeta)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % beta d
        if (par.opt_beta & par.opt_d)
            kk = kk + 1;
            tmp = [Z ; ...
                   Z ; ...
                   Z ; ...
                   -d*PFD_d*PICx(:,pindx.lbeta)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % beta kC_T
        if (par.opt_beta & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lbeta); ...
                   Z ; ...
                   -kC_kC_T*DOCx(:,pindx.lbeta); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % beta kdC
        if (par.opt_beta & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.lbeta); ...
                   Z ; ...
                   -kC_kdC*DOCx(:,pindx.lbeta); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % beta R_Si
        if (par.opt_beta & par.opt_R_Si)
            kk = kk + 1;
            tmp = [-d0((1-sigma)*RR_Si*Gx(:,pindx.lbeta))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_Si*Gx(:,pindx.lbeta))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % beta rR
        if (par.opt_beta & par.opt_rR)
            kk = kk + 1;
            tmp = [-d0((1-sigma)*RR_rR*Gx(:,pindx.lbeta))*C2P; ...
                   Z ; ...
                   Z ; ...
                   d0((1-sigma)*RR_rR*Gx(:,pindx.lbeta))*C2P];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % beta cc
        if (par.opt_beta & par.opt_cc)
            kk = kk + 1;
            tmp = cc*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lbeta))*C2P_cc; ...
                      (1-sigma)*d0(Gx(:,pindx.lbeta))*C2P_cc; ...
                      sigma*d0(Gx(:,pindx.lbeta))*C2P_cc; ...
                      d0((1-sigma)*RR*Gx(:,pindx.lbeta))*C2P_cc];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % beta dd
        if (par.opt_beta & par.opt_dd)
            kk = kk + 1;
            tmp = dd*[-d0((I+(1-sigma)*RR)*Gx(:,pindx.lbeta))*C2P_dd; ...
                      (1-sigma)*d0(Gx(:,pindx.lbeta))*C2P_dd; ...
                      sigma*d0(Gx(:,pindx.lbeta))*C2P_dd; ...
                      d0((1-sigma)*RR*Gx(:,pindx.lbeta))*C2P_dd];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % C model only parameters
        % bC_T bC_T
        if (par.opt_bC_T)
            kk = kk + 1;
            [~,~,Hout] = buildPFD(par,'POC');
            PFD_bm_bm = Hout.PFD_bm_bm;
            par.PFD_bm_bm = PFD_bm_bm;
            tmp = [Z ; ...
                   -PFD_bm_bm*POC-2*PFD_bm*POCx(:,pindx.bC_T);
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC_T bC
        if (par.opt_bC_T & par.opt_bC)
            kk = kk + 1;
            [~,~,Hout] = buildPFD(par,'POC');
            PFD_bm_bb = Hout.PFD_bm_bb;
            par.PFD_bm_bb = PFD_bm_bb;
            tmp =  [Z ; ...
                    [-bC*PFD_bm_bb*POC - ...
                     PFD_bm*POCx(:,pindx.lbC) - ...
                     bC*PFD_bb*POCx(:,pindx.bC_T)];
                    Z; Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bC_T kPIC
        if (par.opt_bC_T & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.bC_T) ; ...
                    -PFD_bm*POCx(:,pindx.lkPIC); ...
                    Z ; ...
                    -kPIC*PICx(:,pindx.bC_T)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC_T d
        if (par.opt_bC_T & par.opt_d)
            kk = kk + 1;
            tmp = [Z ; ...
                   -PFD_bm*POCx(:,pindx.ld); ...
                   Z ; ...
                   -d*PFD_d*PICx(:,pindx.bC_T)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bC_T kC_T
        if (par.opt_bC_T & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.bC_T) ; ...
                   -PFD_bm*POCx(:,pindx.kC_T); ...
                   -kC_kC_T*DOCx(:,pindx.bC_T) ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC_T kdC
        if (par.opt_bC_T & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.bC_T) ; ...
                   -PFD_bm*POCx(:,pindx.lkdC); ...
                   -kC_kdC*DOCx(:,pindx.bC_T) ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC_T R_Si
        if (par.opt_bC_T & par.opt_R_Si)
            kk = kk + 1;
            tmp = [Z ; ...
                   -PFD_bm*POCx(:,pindx.lR_Si); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bC_T rR
        if (par.opt_bC_T & par.opt_rR)
            kk = kk + 1;
            tmp = [Z ; ...
                   -PFD_bm*POCx(:,pindx.lrR); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC_T cc
        if (par.opt_bC_T & par.opt_cc)
            kk = kk + 1;
            tmp = [Z ; ...
                   -PFD_bm*POCx(:,pindx.lcc); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC_T dd
        if (par.opt_bC_T & par.opt_dd)
            kk = kk + 1;
            tmp = [Z ; ...
                   -PFD_bm*POCx(:,pindx.ldd); ...
                   Z ; ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC bC
        if (par.opt_bC)
            kk = kk + 1;
            [~,~,Hout] = buildPFD(par,'POC');
            PFD_bb_bb = Hout.PFD_bb_bb;
            par.PFD_bb_bb = PFD_bb_bb;
            tmp = bC*[Z; ...
                      -PFD_bb*POC-bC*PFD_bb_bb*POC-2* ...
                      PFD_bb*POCx(:,pindx.lbC); ...
                      Z ; ...
                      Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bC kPIC
        if (par.opt_bC & par.opt_kPIC)
            kk = kk + 1;
            tmp =  [kPIC*PICx(:,pindx.lbC) ; ...
                    -PFD_bb*POCx(:,pindx.lkPIC); ...
                    Z ; ...
                    -kPIC*PICx(:,pindx.lbC)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC d
        if (par.opt_bC & par.opt_d)
            kk = kk + 1;
            tmp = [Z ; ...
                   -bC*PFD_bb*POCx(:,pindx.ld); ...
                   Z ; ...
                   -d*PFD_d*PICx(:,pindx.lbC)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bC kC_T
        if (par.opt_bC & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lbC); ...
                   -bC*PFD_bb*POCx(:,pindx.kC_T); ...
                   -kC_kC_T*DOCx(:,pindx.lbC); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC kdC
        if (par.opt_bC & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.lbC); ...
                   -bC*PFD_bb*POCx(:,pindx.lkdC); ...
                   -kC_kdC*DOCx(:,pindx.lbC); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC R_Si
        if (par.opt_bC & par.opt_R_Si)
            kk = kk + 1;
            tmp = bC*[Z ; ...
                      -PFD_bb*POCx(:,pindx.lR_Si); ...
                      Z ; ...
                      Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % bC rR
        if (par.opt_bC & par.opt_rR)
            kk = kk + 1;
            tmp = bC*[Z ; ...
                      -PFD_bb*POCx(:,pindx.lrR); ...
                      Z ; ...
                      Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % bC cc
        if (par.opt_bC & par.opt_cc)
            kk = kk + 1;
            tmp = bC*[Z ; ...
                      -PFD_bb*POCx(:,pindx.lcc); ...
                      Z ; ...
                      Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end 
        
        % bC dd
        if (par.opt_bC & par.opt_dd)
            kk = kk + 1;
            tmp = bC*[Z ; ...
                      -PFD_bb*POCx(:,pindx.ldd); ...
                      Z ; ...
                      Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kPIC kPIC
        if (par.opt_kPIC & par.opt_kPIC)
            kk = kk + 1;
            [~,~,Hout] = buildPFD(par,'PIC');
            PFD_k_k = Hout.PFD_k_k;
            par.PFD_k_k = PFD_k_k;
            tmp =  [kPIC*PIC + 2*kPIC*PICx(:,pindx.lkPIC) ; ...
                    Z ;
                    Z ; ...
                    -kPIC*PIC - 2*kPIC*PICx(:,pindx.lkPIC)];
                     % -kPIC*PFD_k*PIC - kPIC*kPIC*PFD_k_k*PIC + ...
                     % -2*kPIC*PFD_k*PICx(:,pindx.lkPIC)]];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kPIC d
        if (par.opt_kPIC & par.opt_d)
            kk = kk + 1;
            [~,~,Hout] = buildPFD(par,'PIC');
            PFD_d_k = Hout.PFD_d_k;
            par.PFD_d_k = PFD_d_k;
            tmp =  [kPIC*PICx(:,pindx.ld) ; ...
                    Z ;
                    Z ; ...
                    [-kPIC*PICx(:,pindx.ld) + ...
                     -d*PFD_d*PICx(:,pindx.lkPIC)]];
            % -d*kPIC*PFD_d_k*PIC + ...
            % -kPIC*PFD_k*PICx(:,pindx.ld) + ...
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kPIC kC_T
        if (par.opt_kPIC & par.opt_kC_T)
            kk  = kk + 1;
            tmp = [[kC_kC_T*DOCx(:,pindx.lkPIC) + ...
                    kPIC*PICx(:,pindx.kC_T)]; ...
                   Z; ...
                   -kC_kC_T*DOCx(:,pindx.lkPIC); ...
                   -kPIC*PICx(:,pindx.kC_T)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kPIC kdC
        if (par.opt_kPIC & par.opt_kdC)
            kk = kk + 1;
            tmp = [[kC_kdC*DOCx(:,pindx.lkPIC) + ...
                    kPIC*PICx(:,pindx.lkdC)]; ...
                   Z ; ...
                   -kC_kdC*DOCx(:,pindx.lkPIC); ...
                   -kPIC*PICx(:,pindx.lkdC)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kPIC R_Si
        if (par.opt_kPIC & par.opt_R_Si)
            kk = kk + 1;
            tmp = [kPIC*PICx(:,pindx.lR_Si); ...
                   Z; ...
                   Z; ...
                   -kPIC*PICx(:,pindx.lR_Si)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kPIC rR
        if (par.opt_kPIC & par.opt_rR)
            kk = kk + 1;
            tmp = [kPIC*PICx(:,pindx.lrR); ...
                   Z; ...
                   Z; ...
                   -kPIC*PICx(:,pindx.lrR)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kPIC cc
        if (par.opt_kPIC & par.opt_cc)
            kk = kk + 1;
            tmp = [kPIC*PICx(:,pindx.lcc); ...
                   Z; ...
                   Z; ...
                   -kPIC*PICx(:,pindx.lcc)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kPIC dd
        if (par.opt_kPIC & par.opt_dd)
            kk = kk + 1;
            tmp = [kPIC*PICx(:,pindx.ldd) ; ...
                   Z; ...
                   Z; ...
                   -kPIC*PICx(:,pindx.ldd)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % d d
        if (par.opt_d & par.opt_d)
            kk = kk + 1;
            [~,~,Hout] = buildPFD(par,'PIC');
            PFD_d_d = Hout.PFD_d_d;
            par.PFD_d_d = PFD_d_d;
            tmp = d*[Z ; ...
                     Z ; ...
                     Z ; ...
                     [-PFD_d*PIC - d*PFD_d_d*PIC + ...
                      -2*PFD_d*PICx(:,pindx.ld)]];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % d kC_T
        if (par.opt_d & par.opt_kC_T)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.ld); ...
                   Z
                   -kC_kC_T*DOCx(:,pindx.ld); ...
                   -d*PFD_d*PICx(:,pindx.kC_T)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % d kdC
        if (par.opt_d & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kdC*DOCx(:,pindx.ld); ...
                   Z ; ...
                   -kC_kdC*DOCx(:,pindx.ld); ...
                   -d*PFD_d*PICx(:,pindx.lkdC)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % d R_Si
        if (par.opt_d & par.opt_R_Si)
            kk = kk + 1;
            tmp = [Z; ...
                   Z; ...
                   Z; ...
                   -d*PFD_d*PICx(:,pindx.lR_Si)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % d rR
        if (par.opt_d & par.opt_rR)
            kk = kk + 1;
            tmp = [Z; ...
                   Z; ...
                   Z; ...
                   -d*PFD_d*PICx(:,pindx.lrR)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % d cc
        if (par.opt_d & par.opt_cc)
            kk = kk + 1;
            tmp = d*[Z; ...
                     Z; ...
                     Z; ...
                     -PFD_d*PICx(:,pindx.lcc)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % d dd
        if (par.opt_d & par.opt_dd)
            kk = kk + 1;
            tmp = d*[Z; ...
                     Z; ...
                     Z; ...
                     -PFD_d*PICx(:,pindx.ldd)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kC_T kC_T
        if (par.opt_kC_T)
            kk = kk + 1;
            tmp = [2*kC_kC_T*DOCx(:,pindx.kC_T); ...
                   Z ; ...
                   -2*kC_kC_T*DOCx(:,pindx.kC_T); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kC_T kdC
        if (par.opt_kC_T & par.opt_kdC)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lkdC); ...
                   Z ; ...
                   -kC_kC_T*DOCx(:,pindx.lkdC); ...
                   Z] + ... 
                  [kC_kdC*DOCx(:,pindx.kC_T); ...
                   Z ; ...
                   -kC_kdC*DOCx(:,pindx.kC_T); ...
                   Z];

            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kC_T R_Si
        if (par.opt_kC_T & par.opt_R_Si)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lR_Si); ...
                   Z; ...
                   -kC_kC_T*DOCx(:,pindx.lR_Si); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kC_T rR
        if (par.opt_kC_T & par.opt_rR)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lrR); ...
                   Z; ...
                   -kC_kC_T*DOCx(:,pindx.lrR); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kC_T cc
        if (par.opt_kC_T & par.opt_cc)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.lcc); ...
                   Z ; ...
                   -kC_kC_T*DOCx(:,pindx.lcc); ...
                   Z] ;
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kC_T dd
        if (par.opt_kC_T & par.opt_dd)
            kk = kk + 1;
            tmp = [kC_kC_T*DOCx(:,pindx.ldd); ...
                   Z; ...
                   -kC_kC_T*DOCx(:,pindx.ldd); ...
                   Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdC kdC
        if (par.opt_kdC)
            kk = kk + 1;
            tmp = kC_kdC*[DOC + 2*DOCx(:,pindx.lkdC); ...
                          Z ; ...
                          -DOC - 2*DOCx(:,pindx.lkdC); ...
                          Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdC R_Si
        if (par.opt_kdC & par.opt_R_Si)
            kk = kk + 1;
            tmp = kC_kdC*[DOCx(:,pindx.lR_Si); ...
                          Z; ...
                          -DOCx(:,pindx.lR_Si); ...
                          Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % kdC rR
        if (par.opt_kdC & par.opt_rR)
            kk = kk + 1;
            tmp = kC_kdC*[DOCx(:,pindx.lrR); ...
                          Z; ...
                          -DOCx(:,pindx.lrR); ...
                          Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdC cc
        if (par.opt_kdC & par.opt_cc)
            kk = kk + 1;
            tmp = kC_kdC*[DOCx(:,pindx.lcc); ...
                          Z ; ...
                          -DOCx(:,pindx.lcc); ...
                          Z] ;
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % kdC dd
        if (par.opt_kdC & par.opt_dd)
            kk = kk + 1;
            tmp = kC_kdC*[DOCx(:,pindx.ldd); ...
                          Z; ...
                          -DOCx(:,pindx.ldd); ...
                          Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        %%%%%%
        % R_Si R_Si
        if (par.opt_R_Si)
            kk = kk + 1;
            tmp = [-(1-sigma)*RR_Si*(G*C2P); ...
                   Z ; ...
                   Z ; ...
                   (1-sigma)*RR_Si*(G*C2P)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end

        % R_Si rR
        if (par.opt_R_Si & par.opt_rR)
            kk = kk + 1;
            tmp = [Z; Z ; Z ; Z];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % R_Si cc
        if (par.opt_R_Si & par.opt_cc)
            kk = kk + 1;
            tmp = [-cc*(1-sigma)*RR_Si*(G*C2P_cc); ...
                   Z ; ...
                   Z ; ...
                   cc*(1-sigma)*RR_Si*(G*C2P_cc)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % R_Si dd
        if (par.opt_R_Si & par.opt_dd)
            kk = kk + 1;
            tmp = [-dd*(1-sigma)*RR_Si*(G*C2P_dd); ...
                   Z; ...
                   Z; ...
                   dd*(1-sigma)*RR_Si*(G*C2P_dd)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % rR rR
        if (par.opt_rR)
            kk = kk + 1;
            tmp = [-(1-sigma)*RR_rR*(G*C2P); ...
                   Z ; ...
                   Z ; ...
                   (1-sigma)*RR_rR*(G*C2P)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % rR cc
        if (par.opt_rR & par.opt_cc)
            kk = kk + 1;
            tmp = [-cc*(1-sigma)*RR_rR*(G*C2P_cc); ...
                   Z ; ...
                   Z ; ...
                   cc*(1-sigma)*RR_rR*(G*C2P_cc)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % rR dd
        if (par.opt_rR & par.opt_dd)
            kk = kk + 1;
            tmp = [-dd*(1-sigma)*RR_rR*(G*C2P_dd); ...
                   Z; ...
                   Z; ...
                   dd*(1-sigma)*RR_rR*(G*C2P_dd)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        %%%%%%%
        % cc cc
        if (par.opt_cc)
            kk = kk + 1;
            tmp = cc*[-(I+(1-sigma)*RR)*(G*(C2P_cc+cc*C2P_cc_cc)); ...
                      (1-sigma)*G*(C2P_cc+cc*C2P_cc_cc); ...
                      sigma*G*(C2P_cc+cc*C2P_cc_cc); ...
                      (1-sigma)*RR*(G*(C2P_cc+cc*C2P_cc_cc))];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % cc dd
        if (par.opt_cc & par.opt_dd)
            kk = kk + 1;
            tmp = cc*dd*[-(I+(1-sigma)*RR)*(G*C2P_cc_dd); ...
                         (1-sigma)*G*C2P_cc_dd; ...
                         sigma*G*C2P_cc_dd; ...
                         (1-sigma)*RR*(G*C2P_cc_dd)];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
        
        % dd dd
        if (par.opt_dd)
            kk = kk + 1;
            tmp = dd*[-(I+(1-sigma)*RR)*(G*(C2P_dd+dd*C2P_dd_dd)); ...
                      (1-sigma)*G*(C2P_dd+dd*C2P_dd_dd); ...
                      sigma*G*(C2P_dd+dd*C2P_dd_dd); ...
                      (1-sigma)*RR*(G*(C2P_dd+dd*C2P_dd_dd))];
            
            Cxx(:,kk) = mfactor(FD, tmp);
        end
    end
end


