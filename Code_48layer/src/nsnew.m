function [sol, ierr, jac, it_hist, x_hist] = nsnew(x,f,options)
% NSNEW  Newton-Armijo nonlinear solver
%
% Compute Newton direction with backslash
%
% Hybrid of Newton, Shamanskii, Chord
%
% C. T. Kelley, April 1, 2003.
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sol, it_hist, ierr, x_hist] = nsnew(x,f,tol,parms)
%
% inputs:
%        initial iterate = x
%        function = f (2nd output argument is Jacobian)
%        tol = [atol, rtol] relative/absolute
%                           error tolerances
%        parms = [maxit, isham, rsham]
%        maxit = maxmium number of iterations
%                default = 100
%        isham, rsham: The Jacobian matrix is
%                computed and factored after isham
%                updates of x or whenever the ratio
%                of successive l2 norms of the
%                 nonlinear residual exceeds rsham.
%
%            isham = -1, rsham = .5 is the default
%            isham =  1, rsham = 0 is Newton's method,
%            isham = -1, rsham = 1 is the chord method,
%            isham =  m, rsham = 1 is the Shamanskii method with
%                        m steps per Jacobian evaluation
%
%                       The Jacobian is computed and factored
%                       whenever the stepsize
%                       is reduced in the line search.
%
%
% output:
%    sol = solution
%    it_hist = array of iteration history, useful for tables and plots
%                The two columns are the residual norm and
%                number of step size reductions done in the line search.
%
%        ierr = 0 upon successful termination
%        ierr = 1 if after maxit iterations
%             the termination criterion is not satsified
%        ierr = 2 failure in the line search. The iteration
%             is terminated if too many steplength reductions
%             are taken.
%        ierr = 3 error in computing Jacobian
%
%    x_hist = matrix of the entire interation history.
%             The columns are the nonlinear iterates. This
%             is useful for making movies, for example, but
%             can consume way too much storage. This is an
%             OPTIONAL argument. Storage is only allocated
%             if x_hist is in the output argument list.

% set the iteration parameters.
iprint = 0;
iplot = 0;
fid = 1;
maxarm = 100;
maxit = 100;
isham = -1;
rsham = .5;
atol = 1e-8;
rtol = 1e-8;
if( ~isempty(options) )
    if( isfield(options,'iprint') ); iprint = options.iprint; end
    if( isfield(options,'iplot') );  iplot  = options.iplot; end
    if( isfield(options,'fid') );    fid    = options.fid; end
    if( isfield(options,'maxarm') ); maxarm = options.maxarm; end
    if( isfield(options,'maxit') );  maxit  = options.maxit; end
    if( isfield(options,'atol') );   atol   = options.atol; end
    if( isfield(options,'rtol') );   rtol   = options.rtol; end
    if( isfield(options,'isham') );  isham  = options.isham; end
    if( isfield(options,'rsham') );  rsham  = options.rsham; end
end

% initialize flags and counters
ierr = 0;
n = length(x);
if nargout == 5
    x_hist = x;
end
itc = 0;

% evaluate f at the initial iterate and compute the stop tolerance
f0 = feval(f,x);
fnrm = norm(f0);
fnrmo = 1;
itsham = isham;
stop_tol = atol+rtol*fnrm;

% check display options
if iprint == 1
    fprintf(fid,' Newton-Armijo solver \n')
    fprintf(fid,['     Iteration     ||F(x_k)||    Steps   \n'])
    fprintf(fid,'      %4d        %12.4e     %4d  \n',...
            [itc, fnrm, 0] )
end

% check to see if initial iterate also final iterate
if fnrm<=stop_tol
    [f0,jac] = feval(f,x);
end
   
%
% MAIN ITERATION LOOP
%
while(fnrm > stop_tol & itc < maxit)
    
    % keep track of the ratio (rat = fnrm/frnmo) of successive residual
    % norms and the iteration counter (itc)
    rat = fnrm/fnrmo;
    fnrmo = fnrm;
    itc = itc+1;
    
    % evaluate and factor the Jacobian on the first iteration, every
    % isham iterates, or if the ratio of successive residual norm is too
    % large
    if(itc == 1 | rat > rsham | itsham == 0 | armflag == 1)
      clear jac
      itsham = isham;
      jac_age = -1;
      [fv,jac] = feval(f,x);
    end
    itsham = itsham-1;
    
    % compute the Newton direcion
    direction = mfactor(jac,-f0);
    
    % Add one to the age of the Jacobian after the factors have been used
    % in a solve. A fresh Jacobian has an age of -1 at birth.
    jac_age = jac_age+1;
    xold = x; fold = f0;
    
    % line search along Newton direction
    [step,iarm,x,f0,armflag] = armijo(direction,x,f0,f,maxarm);
    
    % If the line search fails and the Jacobian is old, update it.
    % If the Jacobian is fresh; you're dead.
    if armflag == 1  
        if jac_age > 0
            sol = xold;
            x = xold; f0 = fold;
            if iprint==1
                disp('Armijo failure; recompute Jacobian.');
            end
        else
            if iprint==1
                disp('Complete Armijo failure.');
            end
            sol = xold;
            ierr = 2;
            return
        end
    end
    fnrm = norm(f0);
    if nargout == 5, x_hist = [x_hist,x]; end
    rat = fnrm/fnrmo;
    if iprint == 1
        fprintf(fid,'      %4d        %12.4e    %4d  \n', [itc, fnrm, iarm] )
    end
    
    % output iteraton history
    it_hist(itc,:) = [itc fnrm iarm];
end
sol = x;

% on failure, set the error flag
if (fnrm>stop_tol|isnan(fnrm))
    ierr = 1;
end

function [step,iarm,xp,fp,armflag] = armijo(direction,x,f0,f,maxarm)
iarm = 0;
sigma1 = .5;
alpha = 1.d-4;
armflag = 0;
xp = x;
fp = f0; 
%
xold = x;
lambda = 1; lamm = 1; lamc = lambda; iarm = 0;
step = lambda*direction;
xt = x + step;
ft = feval(f,xt);
nft = norm(ft);
nf0 = norm(f0);
ff0 = nf0*nf0;
ffc = nft*nft;
ffm = nft*nft;
while nft >= (1 - alpha*lambda) * nf0
    
    %   Apply the three point parabolic model.
    if iarm == 0
        lambda = sigma1*lambda;
    else
        lambda = parab3p(lamc, lamm, ff0, ffc, ffm);
    end
    
    % Update x; keep the books on lambda.
    step = lambda*direction;
    xt = x + step;
    lamm = lamc;
    lamc = lambda;
    
    % Keep the books on the function norms.
    ft = feval(f,xt);
    nft = norm(ft);
    ffm = ffc;
    ffc = nft*nft;
    iarm = iarm+1;
    if iarm > maxarm
        disp(' Armijo failure, too many reductions ');
        armflag = 1;
        sol = xold;
        return;
    end
end

xp = xt;
fp = ft;


function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
% Apply three-point safeguarded parabolic model for a line search.
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
%
% input:
%       lambdac = current steplength
%       lambdam = previous steplength
%       ff0 = value of \| F(x_c) \|^2
%       ffc = value of \| F(x_c + \lambdac d) \|^2
%       ffm = value of \| F(x_c + \lambdam d) \|^2
%
% output:
%       lambdap = new value of lambda given parabolic model
%
% internal parameters:
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%

% Set internal parameters.
sigma0 = .1; sigma1 = .5;

% Compute coefficients of interpolation polynomial.
%
% p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
%
% d1 = (lambdac - lambdam)*lambdac*lambdam < 0
%      so, if c2 > 0 we have negative curvature and default to
%      lambdap = sigam1 * lambda.
c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
if c2 >= 0
    lambdap = sigma1*lambdac; return
end
c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
lambdap = -c1*.5/c2;
if lambdap < sigma0*lambdac
  lambdap = sigma0*lambdac;
end
if lambdap > sigma1*lambdac
  lambdap = sigma1*lambdac;
end
