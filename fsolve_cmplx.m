function z = fsolve_cmplx(f,z0)
% z = fsolve_cmplx(f,z0);
% wrapper function for matlab's fsolve function; finds the real
% and imaginary part of the root of a function whose argument has
% an imaginary part of O(eps^3)
% inputs
%      f:  function handle
%      z0: initial guess at the root
options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);
%options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-16,'OptimalityTolerance',1e-16);
    i = sqrt(-1);
    x = real(z0);
    y = imag(z0);
    f2 = @(x,y) [real(f(x+i*y));
          imag(f(x+i*y))/eps^3]; % we need to rescale by 1/eps^3 otherwise the imaginary part is too small
    r = fsolve(@(r) f2(r(1),r(2)),[x;y],options);
    z = r(1)+i*r(2);
end
