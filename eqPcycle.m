function [P,Px,parm] = eqPcycle(parm,ip,x)
% ip is the mapping from x to parameter names (see switch below)
% output: P is model prediction of DIP,POP,and DOP
% output: F partial derivative of P model w.r.t. model parameters x
% output: Fxx hessian matrix of P model w.r.t.  model parameters x
% unpack some useful stuff
TRdiv = parm.TRdiv;
M3d   = parm.M3d;
grd   = parm.grd;

iwet  = parm.iwet;
nwet  = parm.nwet; % number of wet points;
I     = parm.I   ; % make an identity matrix;

% unpack the parameters to be optimized
if (nargin>1)
    for ik1 = 1:length(ip)
        switch ip(ik1)
          case 1
            parm.interpp = exp(x(ik1)); % interp of temperature dependence
                                        %
          case 2
            parm.slopep = x(ik1); % slope of temperature dependence

          case 3
            parm.interpc = x(ik1); % interp of temperature dependence
                                   %
          case 4
            parm.slopec  = x(ik1); % slope of temperature dependence
                                  %
          case 5
            parm.sigma = exp(x(ik1)); % fraction of organic P
                                      % allocated to dissolved pool
          case 6
            parm.kappa_dp = exp(x(ik1)); % DOP remineralization
                                         % const.[s^-1];
          case 7
            parm.alpha = exp(x(ik1)); % npp scaling factor for DIP
                                      % uptake rate
          case 8
            parm.beta = exp(x(ik1)); % npp scaling exponent for DIP
                                     % uptake rate
        end
    end
end
slopep = parm.x(2);
sigma  = parm.x(5);

% slopep   = parm.slopep;
interpp  = parm.interpp;
kappa_dp = parm.kappa_dp;
alpha    = parm.alpha;
beta     = parm.beta;

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
    Fx = zeros(3*nwet,length(ip));
    [~,~,dPFDdslope,dPFDdinterp] = buildPFD(M3d,grd,parm,slopep,interpp);
    for ik1 = 1:length(ip)
        switch ip(ik1)
          case 1 % interpp
            Fx(:,ik1) = exp(x(ik1))*[Z;...
                                dPFDdinterp*POP;...
                                Z];
          case 2 % slopep
            Fx(:,ik1) =  [Z;...
                          dPFDdslope*POP;...
                          Z];
          case 5 % sigma
            Fx(:,ik1) =  exp(x(ik1))*[Z;...
                                alpha*L*DIP;...
                                -alpha*L*DIP];
            
          case 6 % kappa_dp
            Fx(:,ik1) = exp(x(ik1))*[-DOP;...
                                Z;...
                                DOP];
          case 7 % alpha
            Fx(:,ik1) = exp(x(ik1))*[L*DIP;...
                                -(1-sigma)*L*DIP;...
                                -sigma*L*DIP];
          case 8 %beta
            dLambdadbeta = 0*Lambda;
            dLambdadbeta(:,:,1) = log(npp).*LAM(:,:,1);
            dLambdadbeta(:,:,2) = log(npp).*LAM(:,:,2);
            iz = find(isinf(dLambdadbeta(:)));
            dLambdadbeta(iz) = 0;
            inan = find(isnan(dLambdadbeta(:)));
            dLambdadbeta(inan) = 0;
            dLdbeta = d0(dLambdadbeta(iwet));
            Fx(:,ik1) = exp(x(ik1))*[ alpha*dLdbeta*DIP;...
                                -(1-sigma)*alpha*dLdbeta*DIP;...
                                -sigma*alpha*dLdbeta*DIP];
            % will need L and dLdbeta for gradients of other
            % biogeochemical cycles
            parm.L = L;
            parm.dLdbeta = dLdbeta;
        end
    end
    % Fx is the derivative of the solution wrt to the parameters
    Px = mfactor(FFp,-Fx);
end
