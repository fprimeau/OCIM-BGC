function [P,Px,Pxx,parm] = eqPcycle(par,parm,x)
% ip is the mapping from x to parameter names (see switch below)
% output: P is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
% unpack some useful stuff
on = true; off = false;
TRdiv = parm.TRdiv;
M3d   = parm.M3d;
grd   = parm.grd;

iwet  = parm.iwet;
nwet  = parm.nwet; % number of wet points;
I     = parm.I   ; % make an identity matrix;

% unpack the parameters to be optimized
if (par.biogeochem.opt_sigma == on)
    lsigma = x(par.pindx.lsigma);
    sigma  = exp(lsigma);
else
    sigma  = par.biogeochem.sigma;
end

if (par.biogeochem.opt_kappa_dp == on)
    lkappa_dp = x(par.pindx.lkappa_dp);
    kappa_dp  = exp(lkappa_dp);
else
    kappa_dp  = par.biogeochem.kappa_dp;
end

if (par.biogeochem.opt_slopep == on)
    lslopep = x(par.pindx.lslopep);
    slopep  = exp(lslopep);
else
    slopep  = par.biogeochem.slopep;
end

if (par.biogeochem.opt_interpp == on)
    linterpp = x(par.pindx.linterpp);
    interpp  = exp(linterpp);
else
    interpp  = par.biogeochem.interpp;
end

if (par.biogeochem.opt_alpha == on)
    lalpha = x(par.pindx.lalpha);
    alpha  = exp(lalpha);
else
    alpha  = par.biogeochem.alpha;
end

if (par.biogeochem.opt_beta == on)
    lbeta = x(par.pindx.lbeta);
    beta  = exp(lbeta);
else
    beta  = par.biogeochem.beta;
end

% fixed parameters
% sigma    = parm.sigma;
DIPbar   = M3d(iwet)*parm.DIPbar;  % gobal arerage PO4 conc.[mmol m^-3];
kappa_g  = parm.kappa_g; % PO4 geological restore const.[s^-1];
kappa_p  = parm.kappa_p; % POP solubilization rate constant
npp      = parm.npp;     % net primary production

% build part of the biological DIP uptake operator
Lambda     = parm.Lambda;
LAM        = 0*M3d;
LAM(:,:,1) = (npp.^beta).*Lambda(:,:,1);
LAM(:,:,2) = (npp.^beta).*Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];
parm.L     = L;
% particle flux
PFD = buildPFD(M3d,grd,parm,slopep,interpp);

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

if (nargout>1)
    %
    % Compute the gradient of the solution wrt the parameters
    %
    Z   = sparse(nwet,1);
    DIP = P(1:nwet);
    POP = P(nwet+1:2*nwet);
    DOP = P(2*nwet+1:end);

    % sigma
    if (par.biogeochem.opt_sigma == on)
        tmp = sigma*[Z; alpha*L*DIP; -alpha*L*DIP];
        Px(:,par.pindx.lsigma) = mfactor(FFp,-tmp);
    end

    % kappa_dp
    if (par.biogeochem.opt_kappa_dp == on)
        tmp = kappa_dp*[-DOP; Z; DOP];
        Px(:,par.pindx.lkappa_dp) = mfactor(FFp,-tmp);
    end

    % slopep
    if (par.biogeochem.opt_slopep == on)
        [~,vout] = buildPFD(M3d,grd,parm,slopep,interpp);
        dPFDdslope = vout.dPFDdslope;
        tmp =  [Z; dPFDdslope*POP; Z];
        Px(:,par.pindx.lslopep) = mfactor(FFp, -tmp);
    end

    % interpp
    if (par.biogeochem.opt_interpp == on)
        [~,vout] = buildPFD(M3d,grd,parm,slopep,interpp);
        dPFDdinterp = vout.dPFDdinterp;
        tmp = interpp*[Z; dPFDdinterp*POP;  Z];
        Px(:,par.pindx.linterpp) = mfactor(FFp,-tmp);
    end
    
    % alpha
    if (par.biogeochem.opt_alpha == on)
        tmp = alpha*[L*DIP;...
                     -(1-sigma)*L*DIP;...
                     -sigma*L*DIP];
        Px(:,par.pindx.lalpha) = mfactor(FFp,-tmp);
    end
    
    %beta
    if (par.biogeochem.opt_beta == on)
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
        parm.L = L;
        parm.dLdbeta = dLdbeta;
        Px(:,par.pindx.lbeta) = mfactor(FFp,-tmp);
    end
end

DIPx = Px(0*nwet+1:1*nwet,:);
POPx = Px(1*nwet+1:2*nwet,:);
DOPx = Px(2*nwet+1:end,:);
if (nargout>2)
    %
    % Compute the hessian of the solution wrt the parameters
    %
    % sigma sigma
    kk = 1;
    if (par.biogeochem.opt_sigma == on)
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
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_kappa_dp == on)
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

    % sigma slopep
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_slopep == on)
        tmp = [Z; Z; Z] + ... % d2Jdsigmadslope * P
              [Z; ... % dJdsigmadPdslope
               sigma*alpha*L*DIPx(:,par.pindx.lslopep); ...
               -sigma*alpha*L*DIPx(:,par.pindx.lslopep)] + ...
              [Z;...  % dJdslopepdPdsigma
               dPFDdslope*POPx(:,par.pindx.lsigma); ...
               Z];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % sigma interpp
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_interpp == on)
        tmp = [Z; Z; Z] + ... % d2Jdsigmadinterp * P
              [Z; ... % dJdsigmadPdslope
               sigma*alpha*L*DIPx(:,par.pindx.linterpp); ...
               -sigma*alpha*L*DIPx(:,par.pindx.linterpp)] + ...
              [Z; ... % dJdinterppdPdsigma
               interpp*dPFDdinterp*POPx(:,par.pindx.lsigma); ...
               Z]; 
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % sigma alpha
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_alpha == on)
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
    if (par.biogeochem.opt_sigma == on & par.biogeochem.opt_beta == on)
        tmp = [Z; ... % d2Jdsigmadbeta
               sigma*alpha*beta*dLdbeta*DIP;...
               -sigma*lalpha*beta*dLdbeta*DIP] + ...
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
    if (par.biogeochem.opt_kappa_dp == on)
        tmp = [-kappa_dp*DOP; ... % d2Jdkappa2 * P
               Z;...
               kappa_dp*DOP] + ...
              2*[-kappa_dp*DOPx(:,par.pindx.lkappa_dp); ... % dJdkappadPdkappa
                 Z;...
                 kappa_dp*DOPx(:,par.pindx.lkappa_dp)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp slopep
    if (par.biogeochem.opt_kappa_dp == on & par.biogeochem.opt_slopep == on)
        tmp = [Z; Z; Z] + ... % d2Jdkappadslope
              [-kappa_dp*DOPx(:,par.pindx.slopep); ... % dJdkappadPdslope
               Z;...
               kappa_dp*DOPx(:,par.pindx.slopep)] + ...
              [Z; ...  % dJdslopedPdkappa
               dPFDdslope*POPx(:,par.pindx.lkappa_dp); ...
               Z]; 
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp interpp
    if (par.biogeochem.opt_kappa_dp == on & par.biogeochem.opt_interpp == on)
        tmp = [Z; Z; Z] + ... % d2Jdkappadinterp
              [-kappa_dp*DOPx(:,par.pindx.linterpp); ... % dJdkappadPdinterp
               Z;...
               kappa_dp*DOPx(:,par.pindx.linterpp)] + ...
              [Z; ...  % dJdinterpdPdkappa
               interpp*dPFDdinterp*POPx(:,par.pindx.lkappa_dp); ...
               Z]; 
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp alpha
    if (par.biogeochem.opt_kappa_dp == on & par.biogeochem.opt_alpha == on)
        tmp = [Z; Z; Z] + ... % d2Jdkappadalpha
              [-kappa_dp*DOPx(:,par.pindx.lalpha); ... % dJdkappadPdalpha
               Z;...
               kappa_dp*DOPx(:,par.pindx.lalpha)] + ...
              [alpha*L*DIPx(:,par.pindx.lkappa_dp);...  % dJdalphadPdkappa
               -(1-sigma)*alpha*L*DIPx(:,par.pindx.lkappa_dp);...
               -sigma*alpha*L*DIPx(:,par.pindx.lkappa_dp)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % kappa_dp beta
    if (par.biogeochem.opt_kappa_dp == on & par.biogeochem.opt_beta == on)
        tmp = [Z; Z; Z] + ... % d2Jdkappadbeta
              [-kappa_dp*DOPx(:,par.pindx.lbeta); ... % dJdkappadPdbeta
               Z;...
               kappa_dp*DOPx(:,par.pindx.lbeta)] + ...
              [alpha*beta*dLdbeta*DIPx(:,par.pindx.lkappa_dp);... % dJdbetadPdkappa
               -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,par.pindx.lkappa_dp);...
               -sigma*alpha*beta*dLdbeta*DIPx(:,par.pindx.lkappa_dp)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % slopep slopep
    if (par.biogeochem.opt_slopep == on)
        [~,vout] = buildPFD(M3d,grd,parm,slopep,interpp);
        d2PFDdslope2 = vout.d2PFDdslope2;
        tmp = [Z; ... % d2Jdslope2
               d2PDFdslope2*DOP; ...
               Z] + ...
              2*[Z; ... % dJdslopedPslope
                 dPFDdslope*POPx(:,par.pindx.slopep);...
                 Z];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
    
    % slopep interpp
    if (par.biogeochem.opt_slopep == on & par.biogeochem.opt_interpp ...
        == on)
        [~,vout] = buildPFD(M3d,grd,parm,slopep,interpp);
        d2PFDdslopedinterp = vout.d2PFDdslopedinterp;
        tmp = [Z; ...  % d2Jdslopedinterp
               d2PFDdslopedinterp*POP; ...
               Z] + ...
              [Z; ... % dJdslopedPdinterp
               dPFDdslope*POPx(:,par.pindx.linterpp);...
               Z] + ...
              [Z; ...  % dJdinterpdPdslope
               dPFDdinterp*POPx(:,par.pindx.slopep); ...
               Z]; 
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % slopep alpha
    if (par.biogeochem.opt_slopep == on & par.biogeochem.opt_alpha == on)
        tmp = [Z; Z; Z] + ... % d2Jdslopedalpha
              [Z; ... % dJdslopedPdalpha
               dPFDdslope*POPx(:,par.pindx.lalpha);...
               Z] + ...
              [alpha*L*DIPx(:,par.pindx.slopep); ...  % dJdalphadPdslope
               -(1-sigma)*alpha*L*DIPx(:,par.pindx.slopep); ...
               -sigma*alpha*L*DIPx(:,par.pindx.slopep)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % slopep beta
    if (par.biogeochem.opt_slopep == on & par.biogeochem.opt_beta == on)
        tmp = [Z; Z; Z] + ... % d2Jdslopedbeta
              [Z; ... % dJdslopedPdbeta
               dPFDdslope*POPx(:,par.pindx.lbeta);...
               Z] + ...
              [alpha*beta*dLdbeta*DIPx(:,par.pindx.slopep);...  % dJdbetadPdkappa
               -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,par.pindx.slopep);...
               -sigma*alpha*beta*dLdbeta*DIPx(:,par.pindx.slopep)];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
    
    % interpp interpp
    if (par.biogeochem.opt_interpp == on)
        [~,vout] = buildPFD(M3d,grd,parm,slopep,interpp);
        d2PFDdinterp2 = vout.d2PFDdinterp2;
        tmp = [Z; ... % d2Jdinterp2
               interpp*interpp*d2PFDdinterp2*POP; ...
               Z] + ...
              [Z; ... % d2Jdinterp2
               interpp*dPFDdinterp*POP; ...
               Z] + ...
              2*[Z; ... % dJdinterpdPinterp
                 interpp*dPFDdinterp*POPx(:,par.pindx.linterpp);...
                 Z];
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end    

    % interpp alpha
    if (par.biogeochem.opt_interpp == on & par.biogeochem.opt_alpha == on)
        tmp = [Z; Z; Z] + ... % d2Jdinterpdalpha
              [Z; ... % dJdinterpdPdalpha
               interpp*dPFDdinterp*POPx(:,par.pindx.lalpha);...
               Z] + ...
              [alpha*L*DIPx(:,par.pindx.linterpp); ... % dJdalphadPdinterp
               -(1-sigma)*alpha*L*DIPx(:,par.pindx.linterpp); ...
               -sigma*alpha*L*DIPx(:,par.pindx.linterpp)];  
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end

    % interpp beta
    if (par.biogeochem.opt_interpp == on & par.biogeochem.opt_beta == on)
        tmp = [Z; Z; Z] + ... % d2Jdinterpdbeta
              [Z; ... % dJdinterpdPdbeta
               interpp*dPFDdinterp*POPx(:,par.pindx.lbeta);...
               Z] + ...
              [alpha*L*DIPx(:,par.pindx.linterpp); ... % dJdbetadPdinterp
               -(1-sigma)*alpha*L*DIPx(:,par.pindx.linterpp); ...
               -sigma*alpha*L*DIPx(:,par.pindx.linterpp)];  
        
        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
    
    % alpha alpha
    if (par.biogeochem.opt_alpha == on)
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
    if (par.biogeochem.opt_alpha == on & par.biogeochem.opt_beta == on)
        tmp = [alpha*beta*dLdbeta*DIP;... % d2Jdbeta2*P
               -(1-sigma)*alpha*beta*dLdbeta*DIP;...
               -sigma*alpha*beta*dLdbeta*DIP] + ...
              [alpha*L*DIPx(:,par.pindx.lbeta);... % dJdalphadPdbeta
               -(1-sigma)*alpha*L*DIPx(:,par.pindx.lbeta);...
               -sigma*alpha*L*DIPx(:,par.pindx.lbeta)] + ...
              [alpha*beta*dLdbeta*DIPx(:,par.pindx.lalpha);... %dJdbetadPdalpha
               -(1-sigma)*alpha*beta*dLdbeta*DIPx(:,par.pindx.lalpha);...
               -sigma*alpha*beta*dLdbeta*DIPx(:,par.pindx.lalpha)];

        Pxx(:,kk) = mfactor(FFp,-tmp);
        kk = kk + 1;
    end
    
    % beta beta
    if (par.biogeochem.opt_beta == on)
        d2Lambdadbetadbeta = 0*Lambda;
        d2Lambdadbetadbeta(:,:,1) = log(npp).*log(npp).*LAM(:,:,1);
        d2Lambdadbetadbeta(:,:,2) = log(npp).*log(npp).*LAM(:,:,2);
        iz = find(isinf(d2Lambdadbetadbeta(:)));
        d2Lambdadbetadbeta(iz) = 0;
        inan = find(isnan(d2Lambdadbetadbeta(:)));
        d2Lambdadbetadbeta(inan) = 0;
        d2Ldbetadbeta = d0(d2Lambdadbetadbeta(iwet));

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