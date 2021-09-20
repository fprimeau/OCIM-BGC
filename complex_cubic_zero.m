function [root] = complex_cubic_zero(c3,c2,c1,c0,x0)
%complex_cubic_zero: only work for complex step functions. dealing with
%tiny imaginary component (eps^3)
%   inputs:
%   x0 = initial guess
%   c3,c2,c,1,c0 = coefficients of cubic function
%   f = c3*x^3 + c2*x^2 + c1*x * c0
%   outputs:
%   root = (real part of root) + (imag part of root)*i

xr = real(x0);
xi = imag(x0);

c3r = real(c3);
c3i = imag(c3);

c2r = real(c2);
c2i = imag(c2);

c1r = real(c1);
c1i = imag(c1);

c0r = real(c0);
c0i = imag(c0);

options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);

fr = @(x) (c3r*x^3 + c2r*x^2 + c1r*x +c0r)/1e-7 ;

[xr, ~, exitflag] = fsolve(fr, xr, options);

xi = -(c3i*xr^3 + c2i*xr^2 + c1i*xr + c0i)/(3*c3r*xr^2 + 2*c2r*xr + c1r);

root = xr + i*xi;


if exitflag<=0
		fprintf('exitflag = % i \n', exitflag)
	end

end
