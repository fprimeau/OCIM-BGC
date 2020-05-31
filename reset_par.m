function  x = reset_par(x, parm, par);
%
on = true; off = false;
%++++++++++ check if the optimization routine suggests strange
% parameter values
A = exist('x0');
if A == 0
    x0 = parm.p0;
end

if (par.biogeochem.opt_sigma == on)
    isigma = par.pindx.lsigma;
    xnew = exp(x(isigma));
    xold = exp(x0(isigma));
    if (xnew > 0.75 | xnew < 0.1)
        x(isigma) = x0(isigma);
    end
    fprintf('current sigma is %3.2e \n', exp(x(isigma)));
end

if (par.biogeochem.opt_slopep == on)
    islopep = par.pindx.lslopep;
    xnew = x(islopep);
    xold = x0(islopep);
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(islopep) = x0(islopep);
    end
    fprintf('current slopep is %3.2e \n', x(islopep));
end

if (par.biogeochem.opt_interpp == on)
    iinterpp = par.pindx.linterpp;
    xnew = exp(x(iinterpp));
    xold = exp(x0(iinterpp));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(iinterpp) = x0(iinterpp);
    end
    fprintf('current interpp is %3.2e \n', exp(x(iinterpp)));
end

if (par.biogeochem.opt_alpha == on)
    ialpha = par.pindx.lalpha;
    xnew = exp(x(ialpha));
    xold = exp(x0(ialpha));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(ialpha) = x0(ialpha);
    end
    fprintf('current alpha is %3.2e \n', exp(x(ialpha)));
end

if (par.biogeochem.opt_beta == on)
    ibeta = par.pindx.lbeta;
    xnew = exp(x(ibeta));
    xold = exp(x0(ibeta));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(ibeta) = x0(ibeta);
    end
    fprintf('current beta is %3.2e \n', exp(x(ibeta)));
end

if (par.biogeochem.opt_kappa_dp == on)
    ikappa_dp = par.pindx.lkappa_dp;
    xnew = exp(x(ikappa_dp));
    xold = exp(x0(ikappa_dp));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(ikappa_dp) = x0(ikappa_dp);
    end
    fprintf('current kappa_dp is %3.2e \n', exp(x(ikappa_dp)));
end

if (par.biogeochem.opt_slopec == on)
    islopec = par.pindx.lslopec;
    xnew = x(islopec);
    xold = x0(islopec);
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(islopec) = x0(islopec);
    end
    fprintf('current slopec is %3.2e \n', x(islopec));
end

if (par.biogeochem.opt_interpc == on)
    iinterpc = par.pindx.linterpc;
    xnew = exp(x(iinterpc));
    xold = exp(x0(iinterpc));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(iinterpc) = x0(iinterpc);
    end
    fprintf('current interpc is %3.2e \n', exp(x(iinterpc)));
end

if (par.biogeochem.opt_d == on)
    id = par.pindx.ld;
    xnew = exp(x(id));
    xold = exp(x0(id));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(id) = x0(id);
    end
    fprintf('current d is %3.2e \n', exp(x(id)));
end

if (par.biogeochem.opt_kappa_dc == on)
    ikappa_dc = par.pindx.lkappa_dc;
    xnew = exp(x(ikappa_dc));
    xold = exp(x0(ikappa_dc));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(ikappa_dc) = x0(ikappa_dc);
    end
    fprintf('current kappa_dc is %3.2e \n', exp(x(ikappa_dc)));
end

if (par.biogeochem.opt_RR == on)
    iRR = par.pindx.lRR;
    xnew = exp(x(iRR));
    xold = exp(x0(iRR));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(iRR) = x0(iRR);
    end
    fprintf('current RR is %3.2e \n', exp(x(iRR)));
end

if (par.biogeochem.opt_slopeo == on)
    islopeo = par.pindx.slopeo;
    xnew = x(islopeo);
    xold = x0(islopeo);
    if (xnew > 2 | xnew < -2);
        x(islopeo) = x0(islopeo);
    end
    fprintf('current slopeo is %3.2e \n', x(islopeo));
end

if (par.biogeochem.opt_interpo == on)
    iinterpo = par.pindx.linterpo;
    xnew = exp(x(iinterpo));
    xold = exp(x0(iinterpo));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(iinterpo) = x0(iinterpo);
    end
    fprintf('current interpo is %3.2e \n', exp(x(iinterpo)));
end

%+++++++++restore the parameter values back to their original ones.
