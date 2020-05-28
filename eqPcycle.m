function [P,Px,parm] = eqPcycle(par,parm,x)
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
        [~,~,dPFDdslope] = buildPFD(M3d,grd,parm,slopep,interpp);
        tmp =  [Z; dPFDdslope*POP; Z];
        Px(:,par.pindx.lslopep) = mfactor(FFp, -tmp);
    end

    % interpp
    if (par.biogeochem.opt_interpp == on)
        [~,~,~,dPFDdinterp] = buildPFD(M3d,grd,parm,slopep,interpp);
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


