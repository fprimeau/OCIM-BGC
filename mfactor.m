function [result] = mfactor (arg1, arg2, arg3)
% MFACTOR factorize a matrix, or use the factors to solve Ax=b.
% Uses sparse LU to factorize A, or uses a previously computed factorization to
% solve a linear system.
%
% Example:
%   F = mfactor (A) ;     % factorizes A into the object F
%   x = mfactor (F,b) ;   % uses F to solve Ax=b
%   norm (A*x-b)
%
% This function is based on LINFACTOR by Timothy A. Davis
%
% Requires MATLAB 7.3 (R2006b) or later.
%
% Copyright 2007, Timothy A. Davis, University of Florida
% VERSION 1.1.0, Nov 1, 2007

if (nargin < 1 | nargin > 3 | nargout > 1)          %#ok
    error ('Usage: F=linfactor(A) or x=linfactor(F,b)') ;
end

if (nargin == 1)
  
  %---------------------------------------------------------------------------
  % F = mfactor (A) ;
  %---------------------------------------------------------------------------
  
  A = arg1 ;
  [m n] = size (A) ;
  if (m ~= n)
    error ('linfactor: A must be square') ;
  end
  
  % check to see if partitioned matrix and factor accordingly
  if isstruct(arg1) 
    % A = [A.P A.Q; A.R A.S];
    [p,s] = size(A.Q);
    result.p = p ;
    result.s = s ;
    result.R = A.R;
    % factor the P block
    [L, U, P, Q, R] = lu (A.P) ;
    % compute (P^-1)*Q
    result.pinvq = Q * (U \ (L \ (P * (R \ A.Q)))) ;
    % store the factorization
    result.P.L = L ;
    result.P.U = U ;
    result.P.P = P ;
    result.P.Q = Q ;
    result.P.R = R ;
    % factor S-R*(P^-1*Q)
    [L, U, P, Q, R] = lu (A.S-A.R*result.pinvq) ;
    result.S.L = L ;
    result.S.U = U ;
    result.S.P = P ;
    result.S.Q = Q ;
    result.S.R = R ;
    result.code = 1 ;
  
  else
    
    % use sparse LU (UMFPACK, with row scaling): L*U = P*(R\A)*Q
    [L, U, P, Q, R] = lu (A) ;
    result.L = L ;
    result.U = U ;
    result.P = P ;
    result.Q = Q ;
    result.R = R ;
    result.kind = 'sparse LU: L*U = P*(R\A)*Q where R is diagonal' ;
    result.code = 2 ;
  end
  
else

  %---------------------------------------------------------------------------
  % x = linfactor (F,b)
  %---------------------------------------------------------------------------
  
  F = arg1 ;
  b = arg2 ;
  
  if nargin<3
    method = 'notranspo' ;
  else
    method = arg3 ;
  end
  
  if F.code == 1
    bp = b(1:F.p);
    bs = b(F.p+1:F.p+F.s);
    %
    pinvb = F.P.Q * (F.P.U \ (F.P.L \ (F.P.P * (F.P.R \ bp)))) ;
    rbp = -F.S.Q * (F.S.U \ (F.S.L \ (F.S.P * (F.S.R \ (F.R*pinvb))))) ;
    pbp = pinvb - F.pinvq*rbp ;
    %
    sbs = F.S.Q * (F.S.U \ (F.S.L \ (F.S.P * (F.S.R \ bs)))) ;
    qbs = -F.pinvq*sbs;
    %
    result = [pbp+qbs;rbp+sbs];
  end
  
  if F.code == 2
    if method == 'notranspo'
      result = F.Q * (F.U \ (F.L \ (F.P * (F.R \ b)))) ;
    else
      result = F.R' \ (F.P' * (F.L' \ (F.U' \ (F.Q' * b)))) ;
    end
  end
  
end

