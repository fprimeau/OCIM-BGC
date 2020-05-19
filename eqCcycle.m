function [C,Cx] = eqCcycle(x,parm,ic)
% ip is the mapping from x to parameter names (see switch below)
% output: C is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
global GC
nx = length(x) ;
parm.nx = nx   ;
% unpack the parameters to be optimized
for ik1 = 1:nx
    switch ic(ik1)
      case 1
        parm.interpp  = exp(x(ik1)) ;
      case 2
        parm.slopep   = x(ik1)      ;
      case 3
        parm.interpc  = exp(x(ik1)) ;
      case 4
        parm.slopec   = x(ik1)      ;
      case 5
        parm.sigma    = exp(x(ik1)) ;
      case 6
        parm.kappa_dp = exp(x(ik1)) ;
      case 7
        parm.alpha    = exp(x(ik1)) ;
      case 8
        parm.beta     = exp(x(ik1)) ;
      case 9
        parm.d        = exp(x(ik1)) ;
      case 10
        parm.kappa_dc = exp(x(ik1)) ;
      case 11
        parm.RR       = exp(x(ik1)) ;
    end
end
parm.slopec = parm.x(4);
parm.sigma  = parm.x(5);

iwet = parm.iwet;
nwet = parm.nwet;

% PME part;
[sst,ss] = PME(parm) ;
parm.ss  = ss        ;
parm.sst = sst       ;

options.atol   = 5e-10 ;
options.rtol   = 5e-10 ;
options.iprint = 1     ;

DIC = GC(0*nwet+1:1*nwet) ; 
POC = GC(1*nwet+1:2*nwet) ;
DOC = GC(2*nwet+1:3*nwet) ;
CaC = GC(3*nwet+1:4*nwet) ;
X0  = [DIC;POC;DOC;CaC]   ;

[C,ierr] = nsnew(X0,@(X) C_eqn(X,parm,x,ic),options) ;

if (ierr ~=0)
    fprintf('eqCcycle did not converge.\n') ;
    keyboard
    else
        % reset the global variable for the next call eqCcycle
    GC = 0.99999*real(C); 
    if nargout>1     
        %
        % Compute the gradient of the solution wrt the parameters
        [F,FD,Cx] = C_eqn(C,parm,x,ic);
        %
    end
end

function [F,FD,Cx] = C_eqn(X,parm,x,ic)    
% unpack some useful stuff
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

save tmpC DIC
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
[G,Gx] = uptake(parm);

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
    % particle flux div_rergence [s^-1];
    [~,dPFDdd] = buildPFD_CaCO3(parm,grd,M3d);
    [~,~,dPFDdslope,dPFDdinterp] = buildPFD(M3d,grd,parm,slopec,interpc);
    Z = zeros(nwet,1);
    % bp
    for ik1 = 1:length(x)
        switch ic(ik1)
          case 1 % interpp
            Fx(:,ik1) = [(1+RR)*Gx(:,ik1) ;...
                         -(1-sigma)*Gx(:,ik1) ;...
                         -sigma*Gx(:,ik1)     ;...
                         -RR*Gx(:,ik1)]   ;
            
          case 2 % slopep
            Fx(:,ik1) = [(1+RR)*Gx(:,ik1) ;...
                         -(1-sigma)*Gx(:,ik1) ;...
                         -sigma*Gx(:,ik1)     ;...
                         -RR*Gx(:,ik1)]   ;
            
          case 3 % interpc
            Fx(:,ik1) = exp(x(ik1)) * ...
                [Z ; ...
                 dPFDdinterp*POC ; ...
                 Z ; ...
                 Z];

          case 4 % slopec
            Fx(:,ik1) = [Z ; ...
                         dPFDdslope*POC ; ...
                         Z ; ...
                         Z];
            
          case 5 % sigma
            Fx(:,ik1) = exp(x(ik1)) *  ...
                [(1+RR)*Gx(:,ik1) ;...
                 -(1-sigma)*Gx(:,ik1) ;...
                 -sigma*Gx(:,ik1)     ;...
                 -RR*Gx(:,ik1)]   ;
            
          case 6 % kappa_dp
            Fx(:,ik1) = [(1+RR)*Gx(:,ik1) ;...
                         -(1-sigma)*Gx(:,ik1) ;...
                         -sigma*Gx(:,ik1)     ;...
                         -RR*Gx(:,ik1)]   ;
            
          case 7 % alpha
            Fx(:,ik1) = [(1+RR)*Gx(:,ik1) ;...
                         -(1-sigma)*Gx(:,ik1) ;...
                         -sigma*Gx(:,ik1)     ;...
                         -RR*Gx(:,ik1)]   ;
            
          case 8 % beta
            Fx(:,ik1) = [(1+RR)*Gx(:,ik1) ;...
                         -(1-sigma)*Gx(:,ik1) ;...
                         -sigma*Gx(:,ik1)     ;...
                         -RR*Gx(:,ik1)]   ;
            
          case 9 % d
            Fx(:,ik1) = exp(x(ik1)) * [Z ; ...
                                Z ; ...
                                Z ; ...
                                dPFDdd*CaC];
            
          case 10 % kappa_dc
            Fx(:,ik1) = exp(x(ik1)) * [-DOC ; ...
                                Z ; ...
                                DOC ; ...
                                Z] ;

          case 11 % RR
            Fx(:,ik1) = exp(x(ik1)) * [G ; ...
                                Z ; ...
                                Z ; ...
                                -G] ;
            
        end
    end
    Cx = mfactor(FD,-Fx);
end
