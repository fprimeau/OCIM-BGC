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
else
    par.bC_T  = par.bC_T;
end

% bC
if (par.opt_bC == on)
    lbC = x(par.pindx.lbC);
    par.bC  = exp(lbC);
else
    par.bC  = par.bC;
end

% d
if (par.opt_d == on)
    ld = x(par.pindx.ld);
    par.d = exp(ld);
else
    par.d = par.d;
end

% kappa_dc
if (par.opt_kappa_dc == on)
    lkappa_dc = x(par.pindx.lkappa_dc);
    par.kappa_dc  = exp(lkappa_dc);
else
    par.kappa_dc  = par.kappa_dc;
end

% RR
if (par.opt_RR == on)
    lRR = x(par.pindx.lRR);
    par.RR = exp(lRR);
else
    par.RR = par.RR;
end

% cc
if (par.opt_cc)
    lcc = x(par.pindx.lcc);
    par.cc = exp(lcc);
else
    par.cc = par.cc;
end 

% dd
if (par.opt_dd)
    ldd = x(par.pindx.ldd);
    par.dd = exp(ldd);
else
    par.dd = par.dd;
end
%
options.iprint = 1;
options.atol = 5e-9 ;
options.rtol = 5e-9 ;

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

DIC   = X(0*nwet+1:1*nwet) ; 
POC   = X(1*nwet+1:2*nwet) ;
DOC   = X(2*nwet+1:3*nwet) ;
PIC   = X(3*nwet+1:4*nwet) ;
%
PO4 = par.po4obs(iwet);
%
% fixed parameters
kappa_p = par.kappa_p ;
% parameters need to be optimized
sigma    = par.sigma   ;
d        = par.d;
RR       = par.RR;
bC_T     = par.bC_T;
bC       = par.bC;
kappa_dc = par.kappa_dc;
alpha    = par.alpha ;
beta     = par.beta  ;
cc       = par.cc;
dd       = par.dd;

C2P = 1./(cc*PO4 + dd);
par.C2P = C2P;
% particle flux div_rergence [s^-1];
PFDa = buildPFD(par,'PIC');
PFDc = buildPFD(par,'POC');
par.PFDa = PFDa;
par.PFDc = PFDc;
par.DIC  = DIC;
% Air-Sea gas exchange
vout = Fsea2air(par, 'CO2');
KG = vout.KG;
KGG = vout.KGG;
JgDIC = vout.JgDIC;
% biological DIC uptake operator
G = uptake_C(par); par.G = G;

eq1 =     (1+RR)*G*C2P + TRdiv*DIC - kappa_dc*DOC - kappa_p*PIC - JgDIC;
eq2 = -(1-sigma)*G*C2P + (PFDc+kappa_p*I)*POC;
eq3 =     -sigma*G*C2P + (TRdiv+kappa_dc*I)*DOC - kappa_p*POC;
eq4 =        -RR*G*C2P + (PFDa+kappa_p*I)*PIC;

F   = [eq1; eq2; eq3; eq4];

if nargout > 1
    % construct the LHS matrix for the offline model
    % disp('Preparing LHS and RHS matrix:')
    
    % colum 1 dFdDIC
    Jc{1,1} = TRdiv - KG;
    Jc{2,1} = 0*I;
    Jc{3,1} = 0*I;
    Jc{4,1} = 0*I;
    
    % colum 2 dFdPOC
    Jc{1,2} = 0*I;
    Jc{2,2} = PFDc+kappa_p*I;
    Jc{3,2} = -kappa_p*I;
    Jc{4,2} = 0*I;
    
    % colum 3 dFdDOC
    Jc{1,3} = -kappa_dc*I;
    Jc{2,3} = 0*I;
    Jc{3,3} = TRdiv+kappa_dc*I;
    Jc{4,3} = 0*I;
    
    % colum 4 dFdPIC
    Jc{1,4} = -kappa_p*I;
    Jc{2,4} = 0*I;
    Jc{3,4} = 0*I;
    Jc{4,4} = PFDa+kappa_p*I;
    
    % factorize Jacobian matrix
    FD = mfactor(cell2mat(Jc));
end 

%% ----------------------------------------------------
if (par.optim == off)
    Cx = [];
elseif (par.optim & nargout > 2)
    pindx = par.pindx;
    Z = sparse(nwet,1);
    dC2Pdcc = -PO4./(cc*PO4 + dd).^2;
    dC2Pddd = -1./(cc*PO4 + dd).^2;
    par.dC2Pdcc = dC2Pdcc;
    par.dC2Pddd = dC2Pddd;
    [~,Gx] = uptake_C(par);
    par.Gx = Gx;
    npx = par.npx;
    % P model parameters
    % sigma
    if (par.opt_sigma == on)
        tmp = [Z; ...
               -sigma*G*C2P; ...
               sigma*G*C2P; ...
               Z] + ...
              [  -(1+RR)*d0(Gx(:,pindx.lsigma))*C2P; ...
                 (1-sigma)*d0(Gx(:,pindx.lsigma))*C2P; ...
                 sigma*d0(Gx(:,pindx.lsigma))*C2P; ...
                 RR*d0(Gx(:,pindx.lsigma))*C2P];
        
        Cx(:,pindx.lsigma) = mfactor(FD, tmp);
    end
    
    % kappa_dp
    if (par.opt_kappa_dp == on)
        tmp = [  -(1+RR)*d0(Gx(:,pindx.lkappa_dp))*C2P; ...
                 (1-sigma)*d0(Gx(:,pindx.lkappa_dp))*C2P; ...
                 sigma*d0(Gx(:,pindx.lkappa_dp))*C2P; ...
                 RR*d0(Gx(:,pindx.lkappa_dp))*C2P];
        
        Cx(:,pindx.lkappa_dp) = mfactor(FD, tmp);
    end
    
    % bP_T
    if (par.opt_bP_T == on)
        tmp = [  -(1+RR)*d0(Gx(:,pindx.bP_T))*C2P; ...
                 (1-sigma)*d0(Gx(:,pindx.bP_T))*C2P; ...
                 sigma*d0(Gx(:,pindx.bP_T))*C2P; ...
                 RR*d0(Gx(:,pindx.bP_T))*C2P];
        
        Cx(:,pindx.bP_T) = mfactor(FD, tmp);
    end
    
    % bP
    if (par.opt_bP == on)
        tmp = [  -(1+RR)*d0(Gx(:,pindx.lbP))*C2P;...
                 (1-sigma)*d0(Gx(:,pindx.lbP))*C2P;...
                 sigma*d0(Gx(:,pindx.lbP))*C2P;...
                 RR*d0(Gx(:,pindx.lbP))*C2P];
        
        Cx(:,pindx.lbP) = mfactor(FD, tmp);
    end
    
    % alpha
    if (par.opt_alpha == on)
        tmp = [  -(1+RR)*d0(Gx(:,pindx.lalpha))*C2P; ...
                 (1-sigma)*d0(Gx(:,pindx.lalpha))*C2P; ...
                 sigma*d0(Gx(:,pindx.lalpha))*C2P; ...
                 RR*d0(Gx(:,pindx.lalpha))*C2P];
        
        Cx(:,pindx.lalpha) = mfactor(FD, tmp);
    end
    
    % beta
    if (par.opt_beta == on)
        tmp = [  -(1+RR)*d0(Gx(:,pindx.lbeta))*C2P;...
                 (1-sigma)*d0(Gx(:,pindx.lbeta))*C2P;...
                 sigma*d0(Gx(:,pindx.lbeta))*C2P;...
                 RR*d0(Gx(:,pindx.lbeta))*C2P];
        
        Cx(:,pindx.lbeta) = mfactor(FD, tmp);
    end
    % -------------------- C parameters ------------------
    % bC_T
    if (par.opt_bC_T == on)
        [~,Gout] = buildPFD(par,'POC');
        PFD_bm = Gout.PFD_bm;
        par.PFD_bm = PFD_bm;
        tmp = [Z; -PFD_bm*POC; Z;  Z];
        
        Cx(:,pindx.bC_T) = mfactor(FD, tmp);
    end
    
    % bC
    if (par.opt_bC == on)
        [~,Gout] = buildPFD(par,'POC');
        PFD_bb = Gout.PFD_bb;
        par.PFD_bb = PFD_bb;
        tmp = bC*[Z; -PFD_bb*POC; Z; Z];
        
        Cx(:,pindx.lbC) = mfactor(FD, tmp);
    end
    
    % d
    if (par.opt_d == on)
        [~,Gout] = buildPFD(par,'PIC');
        PFD_d = Gout.PFD_d;
        par.PFD_d = PFD_d;
        tmp = d*[Z; Z; Z; -PFD_d*PIC];
        
        Cx(:,pindx.ld) = mfactor(FD, tmp);
    end
    
    % kappa_dc
    if (par.opt_kappa_dc == on)
        tmp = kappa_dc*[DOC; Z; -DOC; Z];
        
        Cx(:,pindx.lkappa_dc) = mfactor(FD, tmp);
    end
    
    % RR
    if (par.opt_RR == on)
        tmp = RR*[-G*C2P; Z; Z; G*C2P];
        
        Cx(:,pindx.lRR) = mfactor(FD, tmp);
    end
    
    % cc
    if (par.opt_cc == on)
        tmp = cc*[  -(1+RR)*G*dC2Pdcc; ...
                    (1-sigma)*G*dC2Pdcc; ...
                    sigma*G*dC2Pdcc; ...
                    RR*G*dC2Pdcc];
        
        Cx(:,pindx.lcc) = mfactor(FD, tmp);
    end
    
    % dd
    if (par.opt_dd == on)
        tmp = dd*[  -(1+RR)*G*dC2Pddd; ...
                    (1-sigma)*G*dC2Pddd; ...
                    sigma*G*dC2Pddd; ...
                    RR*G*dC2Pddd];
        
        Cx(:,pindx.ldd) = mfactor(FD, tmp);
    end
end
%% --------------------------------------------------------------
if (par.optim == off)
    Cxx = [];
elseif (par.optim & nargout > 3);
    p2c = cc*PO4 + dd;
    d2C2Pddd2 = 2./p2c.^3;
    d2C2Pdcc2 = (2*PO4.^2)./p2c.^3;
    d2C2Pdccddd = (2*PO4)./p2c.^3;
    [~,~,Gxx] = uptake_C(par);
    par.Gxx = Gxx;
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
                          [  -(1+RR)*d0(Gxx(:,kk))*C2P; ...
                             (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                             sigma*d0(Gxx(:,kk))*C2P; ...
                             RR*d0(Gxx(:,kk))*C2P];
                    
                    Cxx(:,kk) = mfactor(FD, tmp);
                    %pairs not assciated with sigma;
                elseif (jj ~= pindx.lsigma & jk ~= pindx.lsigma)
                    tmp = [  -(1+RR)*d0(Gxx(:,kk))*C2P; ...
                             (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                             sigma*d0(Gxx(:,kk))*C2P; ...
                             RR*d0(Gxx(:,kk))*C2P];
                    
                    Cxx(:,kk) = mfactor(FD, tmp);
                else 
                    tmp = [Z;
                           -sigma*d0(Gx(:,jk))*C2P; ...
                           sigma*d0(Gx(:,jk))*C2P; ...
                           Z] + ...
                          [  -(1+RR)*d0(Gxx(:,kk))*C2P; ...
                             (1-sigma)*d0(Gxx(:,kk))*C2P; ...
                             sigma*d0(Gxx(:,kk))*C2P; ...
                             RR*d0(Gxx(:,kk))*C2P];
                    
                    Cxx(:,kk) = mfactor(FD, tmp);
                end
            else 
                tmp = [  -(1+RR)*d0(Gxx(:,kk))*C2P; ...
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
        kk = kk + 1;
        tmp =  [Z; ...
                -PFD_bm*POCx(:,pindx.lsigma) ; ...
                Z; ...
                Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % sigma bC
    if (par.opt_sigma & par.opt_bC)
        kk = kk + 1;
        tmp = bC * ...
              [Z ; ...
               -PFD_bb*POCx(:,pindx.lsigma); ...
               Z ; ...
               Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % sigma d
    if (par.opt_sigma & par.opt_d)
        kk = kk + 1;
        tmp = d*[Z ; ...
                 Z ; ...
                 Z ; ...
                 -PFD_d*PICx(:,pindx.lsigma)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % sigma kappa_dc
    if (par.opt_sigma & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.lsigma); ...
                        Z; ...
                        -DOCx(:,pindx.lsigma); ...
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
        kk = kk + 1;
        tmp = cc*[  -(1+RR)*d0(Gx(:,pindx.lsigma))*dC2Pdcc; ...
                    (1-sigma)*d0(Gx(:,pindx.lsigma))*dC2Pdcc-sigma*G*dC2Pdcc; ...
                    sigma*d0(Gx(:,pindx.lsigma))*dC2Pdcc+sigma*G*dC2Pdcc; ...
                    RR*d0(Gx(:,pindx.lsigma))*dC2Pdcc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % sigma dd
    if (par.opt_sigma & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[  -(1+RR)*d0(Gx(:,pindx.lsigma))*dC2Pddd; ...
                    (1-sigma)*d0(Gx(:,pindx.lsigma))*dC2Pddd-sigma*G*dC2Pddd; ...
                    sigma*d0(Gx(:,pindx.lsigma))*dC2Pddd+sigma*G*dC2Pddd; ...
                    RR*d0(Gx(:,pindx.lsigma))*dC2Pddd];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dp bC_T
    if (par.opt_kappa_dp & par.opt_bC_T)
        kk = kk + 1;
        tmp =  [Z; ...
                -PFD_bm*POCx(:,pindx.lkappa_dp); ...
                Z; ...
                Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dp bC
    if (par.opt_kappa_dp & par.opt_bC)
        kk = kk + 1;
        tmp = bC*[Z ; ...
                  -PFD_bb*POCx(:,pindx.lkappa_dp); ...
                  Z ; ...
                  Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dp d
    if (par.opt_kappa_dp & par.opt_d)
        kk = kk + 1;
        tmp = d*[Z; ...
                 Z; ...
                 Z; ...
                 -PFD_d*PICx(:,pindx.lkappa_dp)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dp kappa_dc
    if (par.opt_kappa_dp & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.lkappa_dp); ...
                        Z ; ...
                        -DOCx(:,pindx.lkappa_dp); ...
                        Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dp RR
    if (par.opt_kappa_dp & par.opt_RR)
        kk = kk + 1;
        tmp = RR*[-d0(Gx(:,pindx.lkappa_dp))*C2P; ...
                  Z ; ...
                  Z ; ...
                  d0(Gx(:,pindx.lkappa_dp))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dp cc
    if (par.opt_kappa_dp & par.opt_cc)
        kk = kk + 1;
        tmp = cc*[  -(1+RR)*d0(Gx(:,pindx.lkappa_dp))*dC2Pdcc; ...
                    (1-sigma)*d0(Gx(:,pindx.lkappa_dp))*dC2Pdcc; ...
                    sigma*d0(Gx(:,pindx.lkappa_dp))*dC2Pdcc; ...
                    RR*d0(Gx(:,pindx.lkappa_dp))*dC2Pdcc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dp dd
    if (par.opt_kappa_dp & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[  -(1+RR)*d0(Gx(:,pindx.lkappa_dp))*dC2Pddd; ...
                    (1-sigma)*d0(Gx(:,pindx.lkappa_dp))*dC2Pddd; ...
                    sigma*d0(Gx(:,pindx.lkappa_dp))*dC2Pddd; ...
                    RR*d0(Gx(:,pindx.lkappa_dp))*dC2Pddd];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP_T bC_T
    if (par.opt_bP_T & par.opt_bC_T)
        kk = kk + 1;
        tmp = [Z ; ...
               -PFD_bm*POCx(:,pindx.bP_T); ...
               Z ; ...
               Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP_T bC
    if (par.opt_bP_T & par.opt_bC)
        kk = kk + 1;
        tmp = bC*[Z ; ...
                  -PFD_bb*POCx(:,pindx.bP_T); ...
                  Z ; ...
                  Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP_T d
    if (par.opt_bP_T & par.opt_d)
        kk = kk + 1;
        tmp = d*[Z ; ...
                 Z ; ...
                 Z ; ...
                 -PFD_d*PICx(:,pindx.bP_T)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP_T kappa_dc
    if (par.opt_bP_T & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.bP_T); ...
                        Z; ...
                        -DOCx(:,pindx.bP_T); ...
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
        tmp = cc*[  -(1+RR)*d0(Gx(:,pindx.bP_T))*dC2Pdcc; ...
                    (1-sigma)*d0(Gx(:,pindx.bP_T))*dC2Pdcc; ...
                    sigma*d0(Gx(:,pindx.bP_T))*dC2Pdcc; ...
                    RR*d0(Gx(:,pindx.bP_T))*dC2Pdcc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP_T dd
    if (par.opt_bP_T & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[  -(1+RR)*d0(Gx(:,pindx.bP_T))*dC2Pddd; ...
                    (1-sigma)*d0(Gx(:,pindx.bP_T))*dC2Pddd; ...
                    sigma*d0(Gx(:,pindx.bP_T))*dC2Pddd; ...
                    RR*d0(Gx(:,pindx.bP_T))*dC2Pddd];
        
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
        tmp = bC*[Z; ...
                  -PFD_bb*POCx(:,pindx.lbP); ...
                  Z; ...
                  Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP d
    if (par.opt_bP & par.opt_d)
        kk = kk + 1;
        tmp = d*[Z; ...
                 Z; ...
                 Z; ...
                 -PFD_d*PICx(:,pindx.lbP)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP kappa_dc
    if (par.opt_bP & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.lbP); ...
                        Z ; ...
                        -DOCx(:,pindx.lbP); ...
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
        tmp = cc*[  -(1+RR)*d0(Gx(:,pindx.lbP))*dC2Pdcc; ...
                    (1-sigma)*d0(Gx(:,pindx.lbP))*dC2Pdcc;...
                    sigma*d0(Gx(:,pindx.lbP))*dC2Pdcc;...
                    RR*d0(Gx(:,pindx.lbP))*dC2Pdcc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % bP dd
    if (par.opt_bP & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[  -(1+RR)*d0(Gx(:,pindx.lbP))*dC2Pddd; ...
                    (1-sigma)*d0(Gx(:,pindx.lbP))*dC2Pddd;...
                    sigma*d0(Gx(:,pindx.lbP))*dC2Pddd;...
                    RR*d0(Gx(:,pindx.lbP))*dC2Pddd];
        
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
        tmp = bC* ...
              [Z ; ...
               -PFD_bb*POCx(:,pindx.lalpha); ...
               Z ; ...
               Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % alpha d
    if (par.opt_alpha & par.opt_d)
        kk = kk + 1;
        tmp = d*[Z ; ...
                 Z ; ...
                 Z ; ...
                 -PFD_d*PICx(:,pindx.lalpha)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % alpha kappa_dc
    if (par.opt_alpha & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.lalpha); ...
                        Z ; ...
                        -DOCx(:,pindx.lalpha); ...
                        Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % alpha RR
    if (par.opt_alpha & par.opt_RR)
        kk = kk + 1;
        tmp = RR*[-d0(Gx(:,pindx.lalpha))*C2P; ...
                  Z ; ...
                  Z ; ...
                  d0(Gx(:,pindx.lalpha))*C2P];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % alpha cc
    if (par.opt_alpha & par.opt_cc)
        kk = kk + 1;
        tmp = cc*[  -(1+RR)*d0(Gx(:,pindx.lalpha))*dC2Pdcc; ...
                    (1-sigma)*d0(Gx(:,pindx.lalpha))*dC2Pdcc; ...
                    sigma*d0(Gx(:,pindx.lalpha))*dC2Pdcc; ...
                    RR*d0(Gx(:,pindx.lalpha))*dC2Pdcc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % alpha dd
    if (par.opt_alpha & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[  -(1+RR)*d0(Gx(:,pindx.lalpha))*dC2Pddd; ...
                    (1-sigma)*d0(Gx(:,pindx.lalpha))*dC2Pddd; ...
                    sigma*d0(Gx(:,pindx.lalpha))*dC2Pddd; ...
                    RR*d0(Gx(:,pindx.lalpha))*dC2Pddd];
        
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
        tmp = bC* ...
              [Z ; ...
               -PFD_bb*POCx(:,pindx.lbeta); ...
               Z ; ...
               Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % beta d
    if (par.opt_beta & par.opt_d)
        kk = kk + 1;
        tmp = d*[Z ; ...
                 Z ; ...
                 Z ; ...
                 -PFD_d*PICx(:,pindx.lbeta)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % beta kappa_dc
    if (par.opt_beta & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.lbeta); ...
                        Z ; ...
                        -DOCx(:,pindx.lbeta); ...
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
        tmp = cc*[  -(1+RR)*d0(Gx(:,pindx.lbeta))*dC2Pdcc; ...
                    (1-sigma)*d0(Gx(:,pindx.lbeta))*dC2Pdcc; ...
                    sigma*d0(Gx(:,pindx.lbeta))*dC2Pdcc; ...
                    RR*d0(Gx(:,pindx.lbeta))*dC2Pdcc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % beta dd
    if (par.opt_beta & par.opt_dd)
        kk = kk + 1;
        tmp = dd*[  -(1+RR)*d0(Gx(:,pindx.lbeta))*dC2Pddd; ...
                    (1-sigma)*d0(Gx(:,pindx.lbeta))*dC2Pddd; ...
                    sigma*d0(Gx(:,pindx.lbeta))*dC2Pddd; ...
                    RR*d0(Gx(:,pindx.lbeta))*dC2Pddd];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
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
    
    % bC_T kappa_dc
    if (par.opt_bC_T & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = [kappa_dc*DOCx(:,pindx.bC_T) ; ...
               -PFD_bm*POCx(:,pindx.lkappa_dc); ...
               -kappa_dc*DOCx(:,pindx.bC_T) ; ...
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
        tmp = bC* ...
              [Z; ...
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
    
    % bC kappa_dc
    if (par.opt_bC & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = [kappa_dc*DOCx(:,pindx.lbC); ...
               -bC*PFD_bb*POCx(:,pindx.lkappa_dc); ...
               -kappa_dc*DOCx(:,pindx.lbC); ...
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
    
    % d kappa_dc
    if (par.opt_d & par.opt_kappa_dc)
        kk = kk + 1;
        tmp = [kappa_dc*DOCx(:,pindx.ld); ...
               Z ; ...
               -kappa_dc*DOCx(:,pindx.ld); ...
               -d*PFD_d*PICx(:,pindx.lkappa_dc)];
        
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
    
    % kappa_dc kappa_dc
    if (par.opt_kappa_dc)
        kk = kk + 1;
        tmp = kappa_dc*[DOC + 2*DOCx(:,pindx.lkappa_dc); ...
                        Z ; ...
                        -DOC - 2*DOCx(:,pindx.lkappa_dc); ...
                        Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dc RR
    if (par.opt_kappa_dc & par.opt_RR)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.lRR); ...
                        Z; ...
                        -DOCx(:,pindx.lRR); ...
                        Z];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dc cc
    if (par.opt_kappa_dc & par.opt_cc)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.lcc); ...
                        Z ; ...
                        -DOCx(:,pindx.lcc); ...
                        Z] ;
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % kappa_dc dd
    if (par.opt_kappa_dc & par.opt_dd)
        kk = kk + 1;
        tmp = kappa_dc*[DOCx(:,pindx.ldd); ...
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
        tmp = -cc*RR*[G*dC2Pdcc; ...
                      Z ; ...
                      Z ; ...
                      -G*dC2Pdcc];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % RR dd
    if (par.opt_RR & par.opt_dd)
        kk = kk + 1;
        tmp = -dd*RR*[G*dC2Pddd; ...
                      Z; ...
                      Z; ...
                      -G*dC2Pddd];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % cc cc
    if (par.opt_cc)
        kk = kk + 1;
        tmp = cc*[  -(1+RR)*G*(dC2Pdcc+cc*d2C2Pdcc2); ...
                    (1-sigma)*G*(dC2Pdcc+cc*d2C2Pdcc2); ...
                    sigma*G*(dC2Pdcc+cc*d2C2Pdcc2); ...
                    RR*G*(dC2Pdcc+cc*d2C2Pdcc2)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % cc dd
    if (par.opt_cc & par.opt_dd)
        kk = kk + 1;
        tmp = cc*dd*[  -(1+RR)*G*d2C2Pdccddd; ...
                       (1-sigma)*G*d2C2Pdccddd; ...
                       sigma*G*d2C2Pdccddd; ...
                       RR*G*d2C2Pdccddd];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
    
    % dd dd
    if (par.opt_dd)
        kk = kk + 1;
        tmp = dd*[  -(1+RR)*G*(dC2Pddd+dd*d2C2Pddd2); ...
                    (1-sigma)*G*(dC2Pddd+dd*d2C2Pddd2); ...
                    sigma*G*(dC2Pddd+dd*d2C2Pddd2); ...
                    RR*G*(dC2Pddd+dd*d2C2Pddd2)];
        
        Cxx(:,kk) = mfactor(FD, tmp);
    end
end


