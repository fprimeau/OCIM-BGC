function z = fsolve_cmplx(f,z0)
% z = fsolve_cmplx(f,z0);
% wrapper function for matlab's fsolve function; finds the real
% and imaginary part of the root of a function whose argument has
% an imaginary part of O(eps^3)
% inputs
%      f:  function handle
%      z0: initial guess at the root
%options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);
options = optimoptions('fsolve','Display','iter-detailed','FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'MaxIter',200,'MaxFunEvals',200,'SpecifyObjectiveGradient',true,'PlotFcn','optimplotx');
    i = sqrt(-1);
    x = real(z0);
    y = imag(z0)/eps^3/1000;

    % f2 = @(x,y) [real(f(x+i*y));
    %       imag(f(x+i*y))/eps^3]; % we need to rescale by 1/eps^3 otherwise the imaginary part is too small

    [r,fval,exitflag] = fsolve(@(r) f2(r(1),r(2)),[x;y],options);
    z = r(1)+i*r(2)*eps^3*1000;

	if exitflag<=0
		fprintf('exitflag = % i \n', exitflag)
	end

	function [out,Jac] = f2(x,y)
		[fval, J] = f(x+y*i*eps^3*1000);
		out = [real(fval); imag(fval)/eps^3]; % we need to rescale by 1/eps^3 otherwise the imaginary part is too small

		if nargout > 1
			% note: Cauchy–Riemann equations
			J1 = real(J);
	        J2 = -imag(J)*eps^3;
	        J3 = imag(J)/eps^3*1000;
	        J4 = real(J)*1000;
			Jac = [J1, J2 ;
				   J3, J4]  ;
		end
	end
end
