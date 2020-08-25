function [par,P,Px,Pxx] = eqPcycle(x, par)
% ip is the mapping from x to parameter names (see switch below)
% output: P is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
% unpack some useful stuff
on = true; off = false;
TRdiv = par.TRdiv;
M3d   = par.M3d;
grd   = par.grd;

iwet  = par.iwet;
nwet  = par.nwet; % number of wet points;
I     = par.I   ; % make an identity matrix;

% unpack the parameters to be optimized
if (par.opt_sigma == on)
    lsigma = x(par.pindx.lsigma);
    sigma  = exp(lsigma);
else
    sigma  = par.sigma;
end
par.sigma = sigma; % pass parameter to C/Si/O models

%
if (par.opt_kappa_dp == on)
    lkappa_dp = x(par.pindx.lkappa_dp);
    kappa_dp  = exp(lkappa_dp);
else
    kappa_dp  = par.kappa_dp;
end
par.kappa_dp = kappa_dp; % pass parameter to C/Si/O models

%
if (par.opt_bP_T == on)
    bP_T = x(par.pindx.bP_T);
else
    bP_T  = par.bP_T;
end
par.bP_T = bP_T; % pass parameter to C/Si/O models

%
if (par.opt_bP == on)
    lbP = x(par.pindx.lbP);
    bP  = exp(lbP);
else
    bP  = par.bP;
end
par.bP = bP; % pass parameter to C/Si/O models

%
if (par.opt_alpha == on)
    lalpha = x(par.pindx.lalpha);
    alpha  = exp(lalpha);
else
    alpha  = par.alpha;
end
par.alpha = alpha; % pass parameter to C/Si/O models

%
if (par.opt_beta == on)
    lbeta = x(par.pindx.lbeta);
    beta  = exp(lbeta);
else
    beta  = par.beta;
end
par.beta = beta; % pass parameter to C/Si/O models

% fixed parameters
% sigma    = par.sigma;
DIPbar   = M3d(iwet)*par.DIPbar;  % gobal arerage PO4 conc.[mmol m^-3];
kappa_g  = par.kappa_g; % PO4 geological restore const.[s^-1];
kappa_p  = par.kappa_p; % POP solubilization rate constant
npp      = par.npp;     % net primary production

% build part of the biological DIP uptake operator
Lambda     = par.Lambda;
LAM        = 0*M3d;
LAM(:,:,1) = (npp.^beta).*Lambda(:,:,1);
LAM(:,:,2) = (npp.^beta).*Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
par.L     = L;

% particle flux
PFD = buildPFD(par,'POP');

% build Jacobian equations.
% column 1 dF/dDIP
Fp{1,1} = TRdiv+alpha*L+kappa_g*I;
Fp{2,1} = -(1-sigma)*alpha*L;
Fp{3,1} = -sigma*alpha*L;

% column 2 dF/dPOP
Fp{1,2} = 0*I;
Fp{2,2} = PFD+kappa_p*I;
Fp{3,2} = -kappa_p*I;

% column 3 dF/dDOP
Fp{1,3} = -kappa_dp*I;
Fp{2,3} = 0*I;
Fp{3,3} = TRdiv+kappa_dp*I;

% right hand side of phosphate equations
RHS = [kappa_g*DIPbar;...
       sparse(nwet,1);...
       sparse(nwet,1)];

% dP/dt + Fp*P = RHS
% factoring Jacobian matrix
FFp = mfactor(cell2mat(Fp));
% solve for P-cycle model state
P = mfactor(FFp,RHS);

%% ---------------------------------------------------------
if (par.optim == off)
    Px = [];
elseif (par.optim & nargout > 2)
    %
    % Compute the gradient of the solution wrt the parameters
    %
    Z   = sparse(nwet,1);
    DIP = P(1:nwet);
    POP = P(nwet+1:2*nwet);
    DOP = P(2*nwet+1:end);

    % sigma
    if (par.opt_sigma == on)
        tmp = sigma*[Z; alpha*L*DIP; -alpha*L*DIP];
        Px(:,par.pindx.lsigma) = mfactor(FFp,-tmp);
    end
    
    % kappa_dp
    if (par.opt_kappa_dp == on)
        tmp = kappa_dp*[-DOP; Z; DOP];
        Px(:,par.pindx.lkappa_dp) = mfactor(FFp,-tmp);
    end

    % bP_T
    if (par.opt_bP_T == on)
        [~,Gout] = buildPFD(par,'POP');
        PFD_bm = Gout.PFD_bm;
        tmp =  [Z; PFD_bm*POP; Z];
        Px(:,par.pindx.bP_T) = mfactor(FFp, -tmp);
    end

    % bP
    if (par.opt_bP == on)
        [~,Gout] = buildPFD(par,'POP');
        PFD_bb = Gout.PFD_bb;
        tmp = bP*[Z; PFD_bb*POP;  Z];
        Px(:,par.pindx.lbP) = mfactor(FFp,-tmp);
    end
    
    % alpha
    if (par.opt_alpha == on)
        tmp = alpha*[L*DIP;...
                     -(1-sigma)*L*DIP;...
                     -sigma*L*DIP];
        Px(:,par.pindx.lalpha) = mfactor(FFp,-tmp);
    end
    
    %beta
    if (par.opt_beta == on)
        dLambdadbeta = 0*Lambda;
        dLambdadbeta(:,:,1) = log(npp).*LAM(:,:,1);
        dLambdadbeta(:,:,2) = log(npp).*LAM(:,:,2);
        iz = find(isinf(dLambdadbeta(:)));
        dLambdadbeta(iz) = 0;
        inan = find(isnan(dLambdadbeta(:)));
        dLambdadbeta(inan) = 0;
        dLdbeta = d0(dLambdadbeta(iwet));
        tmp = beta*[ alpha*dLdbeta*DIP;...
                     -(1-sigma)*alpha*dLdbeta*DIP;...
                     -sigma*alpha*dLdbeta*DIP];
        % will need L and dLdbeta for gradients of other
        % biogeochemical cycles
        par.L = L;
        par.dLdbeta = dLdbeta;
        Px(:,par.pindx.lbeta) = mfactor(FFp,-tmp);
    end
end
%% ---------------------------------------------------------
if (par.optim == off)
    Pxx = [];
elseif (par.optim & nargout > 3) 
    % Compute the hessian of the solution wrt the parameters
    DIPx = Px(0*nwet+1:1*nwet,:);
    POPx = Px(1*nwet+1:2*nwet,:);
    DOPx = Px(2*nwet+1:end,:);
    % sigma sigma
    kk = 1;
    if (par.opt_sigma == on)
        tmp = [Z; ... % d2Jdsigma2 * DIP 
               sigma*alpha*L*DIP; ...
               -sigma*alpha*L*DIP] + ...
              2*[Z; ... % 2 * dJdigmadPdsigma
                 sigma*alpha*L*DIPx(:,par.pindx.lsigma); ...
                 -sigma*alpha*L*DIPx(:,par.pindx.lsigma)];

        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % sigma kappa_dp
    if (par.opt_sigma == on & par.opt_kappa_dp == on)
        tmp = [Z; Z; Z] + ... % d2Jdsigmadkappa * DIP
              [Z; ... % dJdsigmadPdkappa
               sigma*alpha*L*DIPx(:,par.pindx.lkappa_dp);...
               -sigma*alpha*L*DIPx(:,par.pindx.lkappa_dp)] + ...
              [-kappa_dp*I*DOPx(:,par.pindx.lsigma); ... % dJdkappadPdsigma
               Z; ...
               kappa_dp*I*DOPx(:,par.pindx.lsigma)];

        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % sigma bP_T
    if (par.opt_sigma == on & par.opt_bP_T == on)
        tmp = [Z; Z; Z] + ... % d2Jdsigmadslope * P
              [Z; ... % dJdsigmadPdslope
               sigma*alpha*L*DIPx(:,par.pindx.bP_T); ...
               -sigma*alpha*L*DIPx(:,par.pindx.bP_T)] + ...
              [Z;...  % dJdbP_TdPdsigma
               PFD_bm*POPx(:,par.pindx.lsigma); ...
               Z];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % sigma bP
    if (par.opt_sigma == on & par.opt_bP == on)
        tmp = [Z; Z; Z] + ... % d2Jdsigmadinterp * P
              [Z; ... % dJdsigmadPdslope
               sigma*alpha*L*DIPx(:,par.pindx.lbP); ...
               -sigma*alpha*L*DIPx(:,par.pindx.lbP)] + ...
              [Z; ... % dJdbPdPdsigma
               bP*PFD_bb*POPx(:,par.pindx.lsigma); ...
               Z]; 
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % sigma alpha
    if (par.opt_sigma == on & par.opt_alpha == on)
        tmp = [Z; ... % d2Jdsigmadalpha
               sigma*alpha*L*DIP;...
               -sigma*alpha*L*DIP] + ...
              [Z; ... % dJdsigmadPdalpha
               sigma*alpha*L*DIPx(:,par.pindx.lalpha); ...
               -sigma*alpha*L*DIPx(:,par.pindx.lalpha)] + ...
              [alpha*L*DIPx(:,par.pindx.lsigma);... % dJdalphadsigma
               -(1-sigma)*alpha*L*DIPx(:,par.pindx.lsigma)
               -sigma*alpha*L*DIPx(:,par.pindx.lsigma)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % sigma beta
    if (par.opt_sigma == on & par.opt_beta == on)
        tmp = [Z; ... % d2Jdsigmadbeta
               sigma*alpha*beta*dLdbeta*DIP;...
               -sigma*alpha*beta*dLdbeta*DIP] + ...
              [Z; ... %dJdsigmadPdbeta
               sigma*alpha*L*DIPx(:,par.pindx.lbeta); ...
               -sigma*alpha*L*DIPx(:,par.pindx.lbeta)] + ...
              [beta*alpha*dLdbeta*DIPx(:,par.pindx.lsigma);...  % dJdbetadPdsigma
               -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,par.pindx.lsigma);...
               -sigma*alpha*beta*dLdbeta*DIPx(:,par.pindx.lsigma)];;
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp kappa_dp
    if (par.opt_kappa_dp == on)
        tmp = [-kappa_dp*DOP; ... % d2Jdkappa2 * P
               Z;...
               kappa_dp*DOP] + ...
              2*[-kappa_dp*DOPx(:,par.pindx.lkappa_dp); ... % dJdkappadPdkappa
                 Z;...
                 kappa_dp*DOPx(:,par.pindx.lkappa_dp)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp bP_T
    if (par.opt_kappa_dp == on & par.opt_bP_T == on)
        tmp = [Z; Z; Z] + ... % d2Jdkappadslope
              [-kappa_dp*DOPx(:,par.pindx.bP_T); ... % dJdkappadPdslope
               Z;...
               kappa_dp*DOPx(:,par.pindx.bP_T)] + ...
              [Z; ...  % dJdslopedPdkappa
               PFD_bm*POPx(:,par.pindx.lkappa_dp); ...
               Z]; 
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp bP
    if (par.opt_kappa_dp == on & par.opt_bP == on)
        tmp = [Z; Z; Z] + ... % d2Jdkappadinterp
              [-kappa_dp*DOPx(:,par.pindx.lbP); ... % dJdkappadPdinterp
               Z;...
               kappa_dp*DOPx(:,par.pindx.lbP)] + ...
              [Z; ...  % dJdinterpdPdkappa
               bP*PFD_bb*POPx(:,par.pindx.lkappa_dp); ...
               Z]; 
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp alpha
    if (par.opt_kappa_dp == on & par.opt_alpha == on)
        tmp = [Z; Z; Z] + ... % d2Jdkappadalpha
              kappa_dp*[-DOPx(:,par.pindx.lalpha); ... % dJdkappadPdalpha
                        Z;...
                        DOPx(:,par.pindx.lalpha)] + ...
              alpha*[L*DIPx(:,par.pindx.lkappa_dp);...  % dJdalphadPdkappa
                     -(1-sigma)*L*DIPx(:,par.pindx.lkappa_dp);...
                     -sigma*L*DIPx(:,par.pindx.lkappa_dp)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp beta
    if (par.opt_kappa_dp == on & par.opt_beta == on)
        tmp = [Z; Z; Z] + ... % d2Jdkappadbeta
              kappa_dp*[-DOPx(:,par.pindx.lbeta); ... % dJdkappadPdbeta
                        Z;...
                        DOPx(:,par.pindx.lbeta)] + ...
              beta*[alpha*dLdbeta*DIPx(:,par.pindx.lkappa_dp);... % dJdbetadPdkappa
                    -(1-sigma)*alpha*dLdbeta*DIPx(:,par.pindx.lkappa_dp);...
                    -sigma*alpha*dLdbeta*DIPx(:,par.pindx.lkappa_dp)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % bP_T bP_T
    if (par.opt_bP_T == on)
        [~,~,Hout] = buildPFD(par,'POP');
        PFD_bm_bm = Hout.PFD_bm_bm;
        tmp = [Z; ... % d2Jdslope2
               PFD_bm_bm*POP; ...
               Z] + ...
              2*[Z; ... % dJdslopedPslope
               PFD_bm*POPx(:,par.pindx.bP_T);...
               Z];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
    
    % bP_T bP
    if (par.opt_bP_T == on & par.opt_bP == on)
        [~,~,Hout] = buildPFD(par,'POP');
        PFD_bm_bb = Hout.PFD_bm_bb;
        tmp = bP*[Z; ...  % d2Jdslopedinterp
                  PFD_bm_bb*POP; ...
                  Z] + ...
              [Z; ... % dJdslopedPdinterp
               PFD_bm*POPx(:,par.pindx.lbP);...
               Z] + ...
              bP*[Z; ...  % dJdinterpdPdslope
                       PFD_bb*POPx(:,par.pindx.bP_T); ...
                       Z]; 
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % bP_T alpha
    if (par.opt_bP_T == on & par.opt_alpha == on)
        tmp = [Z; Z; Z] + ... % d2Jdslopedalpha
              [Z; ... % dJdslopedPdalpha
               PFD_bm*POPx(:,par.pindx.lalpha);...
               Z] + ...
              [alpha*L*DIPx(:,par.pindx.bP_T); ...  % dJdalphadPdslope
               -(1-sigma)*alpha*L*DIPx(:,par.pindx.bP_T); ...
               -sigma*alpha*L*DIPx(:,par.pindx.bP_T)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % bP_T beta
    if (par.opt_bP_T == on & par.opt_beta == on)
        tmp = [Z; Z; Z] + ... % d2Jdslopedbeta
              [Z; ... % dJdslopedPdbeta
               PFD_bm*POPx(:,par.pindx.lbeta);...
               Z] + ...
              [alpha*beta*dLdbeta*DIPx(:,par.pindx.bP_T);...  % dJdbetadPdkappa
               -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,par.pindx.bP_T);...
               -sigma*alpha*beta*dLdbeta*DIPx(:,par.pindx.bP_T)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
    
    % bP bP
    if (par.opt_bP == on)
        [~,~,Hout] = buildPFD(par,'POP');
        PFD_bb_bb = Hout.PFD_bb_bb;
        tmp = [Z; ... % d2Jdinterp2
               bP*bP*PFD_bb_bb*POP; ...
               Z] + ...
              [Z; ... % d2Jdinterp2
               bP*PFD_bb*POP; ...
               Z] + ...
              2*[Z; ... % dJdinterpdPinterp
                 bP*PFD_bb*POPx(:,par.pindx.lbP);...
                 Z];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end    

    % bP alpha
    if (par.opt_bP == on & par.opt_alpha == on)
        tmp = [Z; Z; Z] + ... % d2Jdinterpdalpha
              [Z; ... % dJdinterpdPdalpha
               bP*PFD_bb*POPx(:,par.pindx.lalpha);...
               Z] + ...
              [alpha*L*DIPx(:,par.pindx.lbP); ... % dJdalphadPdinterp
               -(1-sigma)*alpha*L*DIPx(:,par.pindx.lbP); ...
               -sigma*alpha*L*DIPx(:,par.pindx.lbP)];  
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % bP beta
    if (par.opt_bP == on & par.opt_beta == on)
        tmp = [Z; ...
               Z; ...
               Z] + ... % d2Jdinterpdbeta
              [Z; ... % dJdinterpdPdbeta
               bP*PFD_bb*POPx(:,par.pindx.lbeta);...
               Z] + ...
              beta*[alpha*dLdbeta*DIPx(:,par.pindx.lbP); ... % dJdbetadPdinterp
                    -(1-sigma)*alpha*dLdbeta*DIPx(:,par.pindx.lbP); ...
                    -sigma*alpha*dLdbeta*DIPx(:,par.pindx.lbP)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
    
    % alpha alpha
    if (par.opt_alpha == on)
        tmp = [alpha*L*DIP;... % d2Jdalpha2*P
               -(1-sigma)*alpha*L*DIP;...
               -alpha*sigma*L*DIP] + ...
              2*[alpha*L*DIPx(:,par.pindx.lalpha); ...  % dJdalphadPdalpha
                 -(1-sigma)*alpha*L*DIPx(:,par.pindx.lalpha); ...
                 -sigma*alpha*L*DIPx(:,par.pindx.lalpha)];

        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % alpha beta
    if (par.opt_alpha == on & par.opt_beta == on)
        tmp = [alpha*beta*dLdbeta*DIP;... % d2Jdbeta2*P
               -(1-sigma)*alpha*beta*dLdbeta*DIP;...
               -sigma*alpha*beta*dLdbeta*DIP] + ...
              alpha*[L*DIPx(:,par.pindx.lbeta);... % dJdalphadPdbeta
                     -(1-sigma)*L*DIPx(:,par.pindx.lbeta);...
                     -sigma*L*DIPx(:,par.pindx.lbeta)] + ...
              beta*[alpha*dLdbeta*DIPx(:,par.pindx.lalpha);... %dJdbetadPdalpha
                    -(1-sigma)*alpha*dLdbeta*DIPx(:,par.pindx.lalpha);...
                    -sigma*alpha*dLdbeta*DIPx(:,par.pindx.lalpha)];

        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
    
    % beta beta
    if (par.opt_beta == on)
        d2Lambdadbetadbeta = 0*Lambda;
        d2Lambdadbetadbeta(:,:,1) = log(npp).*log(npp).*LAM(:,:,1);
        d2Lambdadbetadbeta(:,:,2) = log(npp).*log(npp).*LAM(:,:,2);
        iz = find(isinf(d2Lambdadbetadbeta(:)));
        d2Lambdadbetadbeta(iz) = 0;
        inan = find(isnan(d2Lambdadbetadbeta(:)));
        d2Lambdadbetadbeta(inan) = 0;
        d2Ldbetadbeta = d0(d2Lambdadbetadbeta(iwet));
        par.d2Ldbetadbeta = d2Ldbetadbeta;
        tmp = [alpha*beta*dLdbeta*DIP;... % d2Jdbeta2 * P
               -(1-sigma)*alpha*beta*dLdbeta*DIP;...
               -sigma*alpha*beta*dLdbeta*DIP] + ...
              [alpha*beta*beta*d2Ldbetadbeta*DIP;... % d2Jdbeta2 * P
               -(1-sigma)*alpha*beta*beta*d2Ldbetadbeta*DIP; ...
               -sigma*alpha*beta*beta*d2Ldbetadbeta*DIP] + ...
              2*[ alpha*beta*dLdbeta*DIPx(:,par.pindx.lbeta);...
                  -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,par.pindx.lbeta);...
                  -sigma*alpha*beta*dLdbeta*DIPx(:,par.pindx.lbeta)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
end

