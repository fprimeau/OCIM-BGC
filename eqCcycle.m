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
    lbC = x(par.pindx.lbC);
    par.bC  = exp(lbC);
end

% d
if (par.opt_d == on)
    ld = x(par.pindx.ld);
    par.d = exp(ld);
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

% RR
if (par.opt_RR == on)
    lRR = x(par.pindx.lRR);
    par.RR = exp(lRR);
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
options.iprint = 1 ; 
options.atol = 1e-12 ;
options.rtol = 1e-12 ;

X0  = GC;
[C,ierr] = nsnew(X0,@(X) C_eqn(X,x,par),options) ;
if (ierr ~= 0)
    fprintf('eqCcycle did not converge.\n') ;
    keyboard
else
    % reset the global variable for the next call eqCcycle
    GC = real(C) + 1e-6*randn(4*nwet,1);
    X0 = GC;
    F = C_eqn(C,x,par);
    % test if norm of F small enough, if now rerun nsnew;
    if norm(F) > 1e-12
        [C,ierr] = nsnew(X0,@(X) C_eqn(X,x,par),options);
    end 
    %
    [F,FD,Cx,Cxx,par] = C_eqn(C,x,par);
end
%% --------------------------------------------------------------
function [F,FD,Cx,Cxx,par] = C_eqn(X,x,par)    
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
sigma = par.sigma ;
d     = par.d     ;
RR    = par.RR    ;
bC_T  = par.bC_T  ;
bC    = par.bC    ;
kC_T  = par.kC_T  ;
kdC   = par.kdC   ;
alpha = par.alpha ;
beta  = par.beta  ;
cc    = par.cc    ;
dd    = par.dd    ;

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
% biological DIC uptake operator
G = uptake_C(par); par.G = G;

eq1 =     (1+RR)*G*C2P + TRdiv*DIC - kC*DOC - kappa_p*PIC - JgDIC;
eq2 = -(1-sigma)*G*C2P + (PFDc+kappa_p*I)*POC;
eq3 =     -sigma*G*C2P + (TRdiv+kC)*DOC - kappa_p*POC;
eq4 =        -RR*G*C2P + (PFDa+kappa_p*I)*PIC;

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
    Jc{2,2} = PFDc+kappa_p*I ;
    Jc{3,2} = -kappa_p*I ;
    Jc{4,2} = 0*I ;
    
    % colum 3 dFdDOC
    Jc{1,3} = -kC ;
    Jc{2,3} = 0*I ;
    Jc{3,3} = TRdiv+kC ;
    Jc{4,3} = 0*I ;
    
    % colum 4 dFdPIC
    Jc{1,4} = -kappa_p*I ;
    Jc{2,4} = 0*I ;
    Jc{3,4} = 0*I ;
    Jc{4,4} = PFDa+kappa_p*I ;
    
    % factorize Jacobian matrix
    FD = mfactor(cell2mat(Jc)) ;
end 

%% ----------------------------------------------------
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
        tmp = [Z; ...
               -sigma*G*C2P; ...
               sigma*G*C2P; ...
               Z] + ...
              [-(1+RR)*d0(Gx(:,pindx.lsigma))*C2P; ...
               (1-sigma)*d0(Gx(:,pindx.lsigma))*C2P; ...
               sigma*d0(Gx(:,pindx.lsigma))*C2P; ...
               RR*d0(Gx(:,pindx.lsigma))*C2P];
        
        Cx(:,pindx.lsigma) = mfactor(FD, tmp);
    end

    % kP_T
    if (par.opt_kP_T == on)
        tmp = [-(1+RR)*d0(Gx(:,pindx.kP_T))*C2P; ...
               (1-sigma)*d0(Gx(:,pindx.kP_T))*C2P; ...
               sigma*d0(Gx(:,pindx.kP_T))*C2P; ...
               RR*d0(Gx(:,pindx.kP_T))*C2P];
        
        Cx(:,pindx.kP_T) = mfactor(FD, tmp);
    end
    
    % kdP
    if (par.opt_kdP == on)
        tmp = [-(1+RR)*d0(Gx(:,pindx.lkdP))*C2P; ...
               (1-sigma)*d0(Gx(:,pindx.lkdP))*C2P; ...
               sigma*d0(Gx(:,pindx.lkdP))*C2P; ...
               RR*d0(Gx(:,pindx.lkdP))*C2P];
        
        Cx(:,pindx.lkdP) = mfactor(FD, tmp);
    end
    
    % bP_T
    if (par.opt_bP_T == on)
        tmp = [-(1+RR)*d0(Gx(:,pindx.bP_T))*C2P; ...
               (1-sigma)*d0(Gx(:,pindx.bP_T))*C2P; ...
               sigma*d0(Gx(:,pindx.bP_T))*C2P; ...
               RR*d0(Gx(:,pindx.bP_T))*C2P];
        
        Cx(:,pindx.bP_T) = mfactor(FD, tmp);
    end
    
    % bP
    if (par.opt_bP == on)
        tmp = [-(1+RR)*d0(Gx(:,pindx.lbP))*C2P;...
               (1-sigma)*d0(Gx(:,pindx.lbP))*C2P;...
               sigma*d0(Gx(:,pindx.lbP))*C2P;...
               RR*d0(Gx(:,pindx.lbP))*C2P];
        
        Cx(:,pindx.lbP) = mfactor(FD, tmp);
    end
    
    % alpha
    if (par.opt_alpha == on)
        tmp = [-(1+RR)*d0(Gx(:,pindx.lalpha))*C2P; ...
               (1-sigma)*d0(Gx(:,pindx.lalpha))*C2P; ...
               sigma*d0(Gx(:,pindx.lalpha))*C2P; ...
               RR*d0(Gx(:,pindx.lalpha))*C2P];
        
        Cx(:,pindx.lalpha) = mfactor(FD, tmp);
    end
    
    % beta
    if (par.opt_beta == on)
        tmp = [-(1+RR)*d0(Gx(:,pindx.lbeta))*C2P;...
               (1-sigma)*d0(Gx(:,pindx.lbeta))*C2P;...
               sigma*d0(Gx(:,pindx.lbeta))*C2P;...
               RR*d0(Gx(:,pindx.lbeta))*C2P];
        
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
    
    % RR
    if (par.opt_RR == on)
        tmp = RR*[-G*C2P; Z; Z; G*C2P];
        
        Cx(:,pindx.lRR) = mfactor(FD, tmp);
    end
    
    % cc
    if (par.opt_cc == on)
        tmp = cc*[-(1+RR)*G*C2P_cc; ...
                  (1-sigma)*G*C2P_cc; ...
                  sigma*G*C2P_cc; ...
                  RR*G*C2P_cc];
        
        Cx(:,pindx.lcc) = mfactor(FD, tmp);
    end
    
    % dd
    if (par.opt_dd == on)
        tmp = dd*[-(1+RR)*G*C2P_dd; ...
                  (1-sigma)*G*C2P_dd; ...
                  sigma*G*C2P_dd; ...
                  RR*G*C2P_dd];
        
        Cx(:,pindx.ldd) = mfactor(FD, tmp);
    end
end
%% --------------------------------------------------------------
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
                    tmp = sigma*[Z; ...
                                 -G*C2P - 2*d0(Gx(:,jj))*C2P;...
                                 G*C2P + 2*d0(Gx(:,jj))*C2P;...
                                 Z] + ...
                          [-(1+RR)*d0(Gxx(:,kk))*C2P; ...
                           (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                           sigma*d0(Gxx(:,kk))*C2P; ...
                           RR*d0(Gxx(:,kk))*C2P];
                    
                    Cxx(:,kk) = mfactor(FD, tmp);
                    %pairs not assciated with sigma;
                elseif (jj ~= pindx.lsigma & jk ~= pindx.lsigma)
                    tmp = [-(1+RR)*d0(Gxx(:,kk))*C2P; ...
                           (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                           sigma*d0(Gxx(:,kk))*C2P; ...
                           RR*d0(Gxx(:,kk))*C2P];
                    
                    Cxx(:,kk) = mfactor(FD, tmp);
                else 
                    tmp = [Z;
                           -sigma*d0(Gx(:,jk))*C2P; ...
                           sigma*d0(Gx(:,jk))*C2P; ...
                           Z] + ...
                          [-(1+RR)*d0(Gxx(:,kk))*C2P; ...
                           (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                           sigma*d0(Gxx(:,kk))*C2P; ...
                           RR*d0(Gxx(:,kk))*C2P];
                    
                    Cxx(:,kk) = mfactor(FD, tmp);
                end
            else 
                tmp = [-(1+RR)*d0(Gxx(:,kk))*C2P; ...
                       (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                       sigma*d0(Gxx(:,kk))*C2P; ...
                       RR*d0(Gxx(:,kk))*C2P];
                
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
    
    % sigma RR
    if (par.opt_sigma & par.opt_RR)
        kk = kk + 1;
        tmp = RR*[-d0(Gx(:,pindx.lsigma))*C2P; ...
                  Z ; ...
                  Z ; ...
                  d0(Gx(:,pindx.lsigma))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % sigma cc
    if (par.opt_sigma & par.opt_cc)
        kk  = kk + 1;
        tmp = cc*[-(1+RR)*d0(Gx(:,pindx.lsigma))*C2P_cc; ...
                  (1-sigma)*d0(Gx(:,pindx.lsigma))*C2P_cc-sigma*G*C2P_cc; ...
                  sigma*d0(Gx(:,pindx.lsigma))*C2P_cc+sigma*G*C2P_cc; ...
                  RR*d0(Gx(:,pindx.lsigma))*C2P_cc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % sigma dd
    if (par.opt_sigma & par.opt_dd)
        kk  = kk + 1;
        tmp = dd*[-(1+RR)*d0(Gx(:,pindx.lsigma))*C2P_dd; ...
                  (1-sigma)*d0(Gx(:,pindx.lsigma))*C2P_dd-sigma*G*C2P_dd; ...
                  sigma*d0(Gx(:,pindx.lsigma))*C2P_dd+sigma*G*C2P_dd; ...
                  RR*d0(Gx(:,pindx.lsigma))*C2P_dd];
        
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
    
    % kP_T RR
    if (par.opt_kP_T & par.opt_RR)
        kk = kk + 1;
        tmp = [-RR*d0(Gx(:,pindx.kP_T))*C2P; ...
               Z ; ...
               Z ; ...
               RR*d0(Gx(:,pindx.kP_T))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kP_T cc
    if (par.opt_kP_T & par.opt_cc)
        kk = kk + 1;
        tmp = cc*[-(1+RR)*d0(Gx(:,pindx.kP_T))*C2P_cc; ...
                  (1-sigma)*d0(Gx(:,pindx.kP_T))*C2P_cc; ...
                  sigma*d0(Gx(:,pindx.kP_T))*C2P_cc; ...
                  RR*d0(Gx(:,pindx.kP_T))*C2P_cc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kP_T dd
    if (par.opt_kP_T & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[-(1+RR)*d0(Gx(:,pindx.kP_T))*C2P_dd; ...
                  (1-sigma)*d0(Gx(:,pindx.kP_T))*C2P_dd; ...
                  sigma*d0(Gx(:,pindx.kP_T))*C2P_dd; ...
                  RR*d0(Gx(:,pindx.kP_T))*C2P_dd];
        
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
    
    % kdP RR
    if (par.opt_kdP & par.opt_RR)
        kk  = kk + 1;
        tmp = [-RR*d0(Gx(:,pindx.lkdP))*C2P; ...
               Z ; ...
               Z ; ...
               RR*d0(Gx(:,pindx.lkdP))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kdP cc
    if (par.opt_kdP & par.opt_cc)
        kk  = kk + 1;
        tmp = cc*[-(1+RR)*d0(Gx(:,pindx.lkdP))*C2P_cc; ...
                  (1-sigma)*d0(Gx(:,pindx.lkdP))*C2P_cc; ...
                  sigma*d0(Gx(:,pindx.lkdP))*C2P_cc; ...
                  RR*d0(Gx(:,pindx.lkdP))*C2P_cc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kdP dd
    if (par.opt_kdP & par.opt_dd)
        kk  = kk + 1;
        tmp = dd*[-(1+RR)*d0(Gx(:,pindx.lkdP))*C2P_dd; ...
                  (1-sigma)*d0(Gx(:,pindx.lkdP))*C2P_dd; ...
                  sigma*d0(Gx(:,pindx.lkdP))*C2P_dd; ...
                  RR*d0(Gx(:,pindx.lkdP))*C2P_dd];
        
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
    
    % bP_T RR
    if (par.opt_bP_T & par.opt_RR)
        kk = kk + 1;
        tmp = RR*[-d0(Gx(:,pindx.bP_T))*C2P; ...
                  Z ; ...
                  Z ; ...
                  d0(Gx(:,pindx.bP_T))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP_T cc
    if (par.opt_bP_T & par.opt_cc)
        kk = kk + 1;
        tmp = cc*[-(1+RR)*d0(Gx(:,pindx.bP_T))*C2P_cc; ...
                  (1-sigma)*d0(Gx(:,pindx.bP_T))*C2P_cc; ...
                  sigma*d0(Gx(:,pindx.bP_T))*C2P_cc; ...
                  RR*d0(Gx(:,pindx.bP_T))*C2P_cc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP_T dd
    if (par.opt_bP_T & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[-(1+RR)*d0(Gx(:,pindx.bP_T))*C2P_dd; ...
                  (1-sigma)*d0(Gx(:,pindx.bP_T))*C2P_dd; ...
                  sigma*d0(Gx(:,pindx.bP_T))*C2P_dd; ...
                  RR*d0(Gx(:,pindx.bP_T))*C2P_dd];
        
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
    
    % bP RR
    if (par.opt_bP & par.opt_RR)
        kk = kk + 1;
        tmp = RR*[-d0(Gx(:,pindx.lbP))*C2P; ...
                  Z ; ...
                  Z ; ...
                  d0(Gx(:,pindx.lbP))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP cc
    if (par.opt_bP & par.opt_cc)
        kk = kk + 1;
        tmp = cc*[-(1+RR)*d0(Gx(:,pindx.lbP))*C2P_cc; ...
                  (1-sigma)*d0(Gx(:,pindx.lbP))*C2P_cc;...
                  sigma*d0(Gx(:,pindx.lbP))*C2P_cc;...
                  RR*d0(Gx(:,pindx.lbP))*C2P_cc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP dd
    if (par.opt_bP & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[-(1+RR)*d0(Gx(:,pindx.lbP))*C2P_dd; ...
                  (1-sigma)*d0(Gx(:,pindx.lbP))*C2P_dd;...
                  sigma*d0(Gx(:,pindx.lbP))*C2P_dd;...
                  RR*d0(Gx(:,pindx.lbP))*C2P_dd];
        
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
    
    % alpha RR
    if (par.opt_alpha & par.opt_RR)
        kk  = kk + 1;
        tmp = RR*[-d0(Gx(:,pindx.lalpha))*C2P; ...
                  Z ; ...
                  Z ; ...
                  d0(Gx(:,pindx.lalpha))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % alpha cc
    if (par.opt_alpha & par.opt_cc)
        kk  = kk + 1;
        tmp = cc*[-(1+RR)*d0(Gx(:,pindx.lalpha))*C2P_cc; ...
                  (1-sigma)*d0(Gx(:,pindx.lalpha))*C2P_cc; ...
                  sigma*d0(Gx(:,pindx.lalpha))*C2P_cc; ...
                  RR*d0(Gx(:,pindx.lalpha))*C2P_cc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % alpha dd
    if (par.opt_alpha & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[-(1+RR)*d0(Gx(:,pindx.lalpha))*C2P_dd; ...
                  (1-sigma)*d0(Gx(:,pindx.lalpha))*C2P_dd; ...
                  sigma*d0(Gx(:,pindx.lalpha))*C2P_dd; ...
                  RR*d0(Gx(:,pindx.lalpha))*C2P_dd];
        
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
    
    % beta RR
    if (par.opt_beta & par.opt_RR)
        kk = kk + 1;
        tmp = RR*[-d0(Gx(:,pindx.lbeta))*C2P; ...
                  Z ; ...
                  Z ; ...
                  d0(Gx(:,pindx.lbeta))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % beta cc
    if (par.opt_beta & par.opt_cc)
        kk = kk + 1;
        tmp = cc*[-(1+RR)*d0(Gx(:,pindx.lbeta))*C2P_cc; ...
                  (1-sigma)*d0(Gx(:,pindx.lbeta))*C2P_cc; ...
                  sigma*d0(Gx(:,pindx.lbeta))*C2P_cc; ...
                  RR*d0(Gx(:,pindx.lbeta))*C2P_cc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % beta dd
    if (par.opt_beta & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[-(1+RR)*d0(Gx(:,pindx.lbeta))*C2P_dd; ...
                  (1-sigma)*d0(Gx(:,pindx.lbeta))*C2P_dd; ...
                  sigma*d0(Gx(:,pindx.lbeta))*C2P_dd; ...
                  RR*d0(Gx(:,pindx.lbeta))*C2P_dd];
        
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
    
    % bC_T RR
    if (par.opt_bC_T & par.opt_RR)
        kk = kk + 1;
        tmp = [Z ; ...
               -PFD_bm*POCx(:,pindx.lRR); ...
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
    
    % bC RR
    if (par.opt_bC & par.opt_RR)
        kk = kk + 1;
        tmp = bC*[Z ; ...
                  -PFD_bb*POCx(:,pindx.lRR); ...
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
    
    % d d
    if (par.opt_d & par.opt_d)
        kk = kk + 1;
        [~,~,Hout] = buildPFD(par,'PIC');
        PFD_d_d = Hout.PFD_d_d;
        par.PFD_d_d = PFD_d_d;
        tmp = d*[Z ; ...
                 Z ; ...
                 Z ; ...
                 -PFD_d*PIC-d*PFD_d_d*PIC-2*PFD_d*PICx(:,pindx.ld)];
        
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
    
    % d RR
    if (par.opt_d & par.opt_RR)
        kk = kk + 1;
        tmp = [Z; ...
               Z; ...
               Z; ...
               -d*PFD_d*PICx(:,pindx.lRR)];
        
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

    % kC_T RR
    if (par.opt_kC_T & par.opt_RR)
        kk = kk + 1;
        tmp = [kC_kC_T*DOCx(:,pindx.lRR); ...
               Z; ...
               -kC_kC_T*DOCx(:,pindx.lRR); ...
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
    
    % kdC RR
    if (par.opt_kdC & par.opt_RR)
        kk = kk + 1;
        tmp = kC_kdC*[DOCx(:,pindx.lRR); ...
                      Z; ...
                      -DOCx(:,pindx.lRR); ...
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
    
    % RR RR
    if (par.opt_RR)
        kk = kk + 1;
        tmp = -RR*[G*C2P; ...
                   Z ; ...
                   Z ; ...
                   -G*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % RR cc
    if (par.opt_RR & par.opt_cc)
        kk = kk + 1;
        tmp = -cc*RR*[G*C2P_cc; ...
                      Z ; ...
                      Z ; ...
                      -G*C2P_cc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % RR dd
    if (par.opt_RR & par.opt_dd)
        kk = kk + 1;
        tmp = -dd*RR*[G*C2P_dd; ...
                      Z; ...
                      Z; ...
                      -G*C2P_dd];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % cc cc
    if (par.opt_cc)
        kk = kk + 1;
        tmp = cc*[-(1+RR)*G*(C2P_cc+cc*C2P_cc_cc); ...
                  (1-sigma)*G*(C2P_cc+cc*C2P_cc_cc); ...
                  sigma*G*(C2P_cc+cc*C2P_cc_cc); ...
                  RR*G*(C2P_cc+cc*C2P_cc_cc)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % cc dd
    if (par.opt_cc & par.opt_dd)
        kk = kk + 1;
        tmp = cc*dd*[-(1+RR)*G*C2P_cc_dd; ...
                     (1-sigma)*G*C2P_cc_dd; ...
                     sigma*G*C2P_cc_dd; ...
                     RR*G*C2P_cc_dd];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % dd dd
    if (par.opt_dd)
        kk = kk + 1;
        tmp = dd*[-(1+RR)*G*(C2P_dd+dd*C2P_dd_dd); ...
                  (1-sigma)*G*(C2P_dd+dd*C2P_dd_dd); ...
                  sigma*G*(C2P_dd+dd*C2P_dd_dd); ...
                  RR*G*(C2P_dd+dd*C2P_dd_dd)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
end


