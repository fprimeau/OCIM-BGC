function [C, Cx, parm] = eqCcycle(par, parm, x)
% ip is the mapping from x to parameter names (see switch below)
% output: C is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
on = true; off = false;
global GC
nx = length(x) ;
parm.nx = nx   ;
iwet = parm.iwet;
nwet = parm.nwet;

% unpack the parameters to be optimized
%sigma
if (par.biogeochem.opt_sigma == on)
    lsigma = x(par.pindx.lsigma);
    parm.sigma  = exp(lsigma);
else
    parm.sigma  = par.biogeochem.sigma;
end

% kappa_dp
if (par.biogeochem.opt_kappa_dp == on)
    lkappa_dp = x(par.pindx.lkappa_dp);
    parm.kappa_dp  = exp(lkappa_dp);
else
    parm.kappa_dp  = par.biogeochem.kappa_dp;
end

% slopep
if (par.biogeochem.opt_slopep == on)
    slopep = x(par.pindx.slopep);
else
    parm.slopep  = par.biogeochem.slopep;
end

% interpp
if (par.biogeochem.opt_interpp == on)
    linterpp = x(par.pindx.linterpp);
    parm.interpp  = exp(linterpp);
else
    parm.interpp  = par.biogeochem.interpp;
end

% alpha
if (par.biogeochem.opt_alpha == on)
    lalpha = x(par.pindx.lalpha);
    parm.alpha  = exp(lalpha);
else
    parm.alpha  = par.biogeochem.alpha;
end

% beta
if (par.biogeochem.opt_beta == on)
    lbeta = x(par.pindx.lbeta);
    parm.beta  = exp(lbeta);
else
    parm.beta  = par.biogeochem.beta;
end

% kappa_dc
if (par.biogeochem.opt_kappa_dc == on)
    lkappa_dc = x(par.pindx.lkappa_dc);
    parm.kappa_dc  = exp(lkappa_dc);
else
    parm.kappa_dc  = par.biogeochem.kappa_dc;
end

% slopec
if (par.biogeochem.opt_slopec == on)
    slopec = x(par.pindx.slopec);
else
    parm.slopec  = par.biogeochem.slopec;
end

% interpc
if (par.biogeochem.opt_interpc == on)
    linterpc = x(par.pindx.linterpc);
    parm.interpc  = exp(linterpc);
else
    parm.interpc  = par.biogeochem.interpc;
end

% d
if (par.biogeochem.opt_d == on)
    ld = x(par.pindx.ld);
    parm.d  = exp(ld);
else
    parm.d = par.biogeochem.d;
end

% RR
if (par.biogeochem.opt_RR == on)
    lRR = x(par.pindx.lRR);
    parm.RR = exp(lRR);
else
    parm.RR = par.biogeochem.RR;
end

options.iprint = 1;
options.atol = 5e-10 ;
options.rtol = 5e-10 ;

DIC = GC(0*nwet+1:1*nwet) ; 
POC = GC(1*nwet+1:2*nwet) ;
DOC = GC(2*nwet+1:3*nwet) ;
CaC = GC(3*nwet+1:4*nwet) ;
X0  = [DIC;POC;DOC;CaC]   ;

[C,ierr] = nsnew(X0,@(X) C_eqn(X,parm,x,par),options) ;

if (ierr ~=0)
    fprintf('eqCcycle did not converge.\n') ;
    keyboard
    else
        % reset the global variable for the next call eqCcycle
        GC = real(C) + 1e-7*randn(4*nwet,1);
    if nargout>1     
        %
        % Compute the gradient of the solution wrt the parameters
        [F,FD,Cx] = C_eqn(C,parm,x,par);
        %
    end
end

function [F,FD,Cx] = C_eqn(X,parm,x,par)    
% unpack some useful stuff
on = true; off = false;
nx    = parm.nx    ;
grd   = parm.grd   ;
dVt   = parm.dVt   ;
M3d   = parm.M3d   ;
TRdiv = parm.TRdiv ;
iwet  = parm.iwet  ;
nwet  = parm.nwet  ;
I     = parm.I     ;

DIC   = X(0*nwet+1:1*nwet) ; 
POC   = X(1*nwet+1:2*nwet) ;
DOC   = X(2*nwet+1:3*nwet) ;
CaC   = X(3*nwet+1:4*nwet) ;

save tmpC DIC POC DOC CaC
% fixed parameters
kappa_p = parm.kappa_p ;  
sigma   = parm.sigma   ;

% parameters need to be optimized
d        = parm.d;
RR       = parm.RR;
slopec   = parm.slopec;
interpc  = parm.interpc;
kappa_dc = parm.kappa_dc;
alpha    = parm.alpha ;
beta     = parm.beta  ;

% particle flux div_rergence [s^-1];
PFDa = buildPFD_CaCO3(parm,grd,M3d);
PFDc = buildPFD(M3d,grd,parm,slopec,interpc);
% Air-Sea gas exchange %%%%%%%%%%%%%%%%%%%%%
[JgDIC,KG] = Fsea2air(parm,DIC);

% biological DIC uptake operator
[G,Gx] = uptake(parm, par);

eq1 =   (1 + RR)*G + TRdiv * DIC - kappa_dc * DOC - kappa_p*CaC - JgDIC;
eq2 = -(1-sigma)*G + (PFDc + kappa_p*I)   *   POC;
eq3 =     -sigma*G + (TRdiv + kappa_dc*I) *   DOC - kappa_p*POC;
eq4 =        -RR*G + (PFDa  + kappa_p*I)  *   CaC;

F   = [eq1;eq2;eq3;eq4];

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

% colum 4 dFdCaC
Jc{1,4} = -kappa_p*I;
Jc{2,4} = 0*I;
Jc{3,4} = 0*I;
Jc{4,4} = PFDa+kappa_p*I;

% factorize Jacobian matrix
FD = mfactor(cell2mat(Jc));

if nargout > 2
    
    Z = zeros(nwet,1);
    
    % interpp
    if (par.biogeochem.opt_interpp == on)
        tmp = [(1+RR)*Gx(:,par.pindx.linterpp) ;...
               -(1-sigma)*Gx(:,par.pindx.linterpp) ;...
               -sigma*Gx(:,par.pindx.linterpp)     ;...
               -RR*Gx(:,par.pindx.linterpp)]   ;
        Cx(:,par.pindx.linterpp) = mfactor(FD, -tmp);
    end
    
    % slopep
    if (par.biogeochem.opt_slopep == on)
        tmp = [(1+RR)*Gx(:,par.pindx.slopp) ;...
               -(1-sigma)*Gx(:,par.pindx.slopp) ;...
               -sigma*Gx(:,par.pindx.slopp)     ;...
               -RR*Gx(:,par.pindx.slopp)]   ;
        Cx(:,par.pindx.slopp) = mfactor(FD, -tmp);
    end
    
    % interpc
    if (par.biogeochem.opt_interpc == on)
        [~,vout] = buildPFD(M3d,grd,parm,slopec,interpc);
        dPFDdinterp = vout.dPFDdinterp;
        tmp = interpc * ...
              [Z ; ...
               dPFDdinterp*POC ; ...
               Z ; ...
               Z];
        Cx(:,par.pindx.linterpc) = mfactor(FD, -tmp);
    end
    
    % slopec
     if (par.biogeochem.opt_slopec == on)
        tmp = [Z ; ...
               dPFDdslope*POC ; ...
               Z ; ...
               Z];
        Cx(:,par.pindx.slopec) = mfactor(FD, -tmp);
    end
    
    % sigma
    if (par.biogeochem.opt_sigma == on)
        tmp = sigma *  ...
              [(1+RR)*Gx(:,par.pindx.lsigma) ;...
               -(1-sigma)*Gx(:,par.pindx.lsigma) + G;...
               -sigma*Gx(:,par.pindx.lsigma) - G  ;...
               -RR*Gx(:,par.pindx.lsigma)]   ;
        Cx(:,par.pindx.lsigma) = mfactor(FD, -tmp);
    end
    
    % kappa_dp
    if (par.biogeochem.opt_kappa_dp == on)
        tmp = [(1+RR)*Gx(:,par.pindx.lkappa_dp) ;...
               -(1-sigma)*Gx(:,par.pindx.lkappa_dp) ;...
               -sigma*Gx(:,par.pindx.lkappa_dp)     ;...
               -RR*Gx(:,par.pindx.lkappa_dp)]   ;
        Cx(:,par.pindx.lkappa_dp) = mfactor(FD, -tmp);
    end
    
    % alpha
    if (par.biogeochem.opt_alpha == on)
        tmp = [(1+RR)*Gx(:,par.pindx.lalpha) ;...
               -(1-sigma)*Gx(:,par.pindx.lalpha) ;...
               -sigma*Gx(:,par.pindx.lalpha)     ;...
               -RR*Gx(:,par.pindx.lalpha)]   ;
        Cx(:,par.pindx.lalpha) = mfactor(FD, -tmp);
    end
    
    % beta
    if (par.biogeochem.opt_beta == on)
        tmp = [(1+RR)*Gx(:,par.pindx.lbeta) ;...
               -(1-sigma)*Gx(:,par.pindx.lbeta) ;...
               -sigma*Gx(:,par.pindx.lbeta)     ;...
               -RR*Gx(:,par.pindx.lbeta)]   ;
        Cx(:,par.pindx.lbeta) = mfactor(FD, -tmp);
    end
    
    % d
    if (par.biogeochem.opt_d == on)
        [~,dPFDdd] = buildPFD_CaCO3(parm,grd,M3d);
        tmp = d * [Z ; ...
                   Z ; ...
                   Z ; ...
                   dPFDdd*CaC];
        Cx(:,par.pindx.ld) = mfactor(FD, -tmp);
    end
    
    % kappa_dc
    if (par.biogeochem.opt_kappa_dc == on)
        tmp = kappa_dc * [-DOC ; ...
                          Z ; ...
                          DOC ; ...
                          Z] ;
        Cx(:,par.pindx.lkappa_dc) = mfactor(FD, -tmp);
    end
    
    % RR
    if (par.biogeochem.opt_RR == on)
        tmp = RR* [G ; ...
                   Z ; ...
                   Z ; ...
                   -G] ;
        Cx(:,par.pindx.lRR) = mfactor(FD, -tmp);
    end            
end

