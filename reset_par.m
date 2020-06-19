function  [parm, x] = reset_par(x, parm, par);
%
on = true; off = false;
%++++++++++ check if the optimization routine suggests strange
% parameter values
A = exist('x0');
if A == 0
    x0 = parm.p0;
end

if (par.opt_sigma == on)
    isigma = par.pindx.lsigma;
    xnew = exp(x(isigma));
    xold = exp(x0(isigma));
    if (xnew > 0.75 | xnew < 0.1)
        x(isigma) = x0(isigma);
    end
    fprintf('current sigma is %3.2e \n', exp(x(isigma)));
end

if (par.opt_kappa_dp == on)
    ikappa_dp = par.pindx.lkappa_dp;
    xnew = exp(x(ikappa_dp));
    xold = exp(x0(ikappa_dp));
    if (xnew > 5*xold | xnew < 0.2*xold);
        x(ikappa_dp) = x0(ikappa_dp);
    end
    fprintf('current kappa_dp is %3.2e \n', exp(x(ikappa_dp)));
end

if (par.opt_slopep == on)
    islopep = par.pindx.slopep;
    xnew = x(islopep);
    xold = x0(islopep);
    if (xnew > 1e-2 | xnew < -1e-2);
        x(islopep) = x0(islopep);
    end
    fprintf('current slopep is %3.2e \n', x(islopep));
end

if (par.opt_interpp == on)
    iinterpp = par.pindx.linterpp;
    xnew = exp(x(iinterpp));
    xold = exp(x0(iinterpp));
    if (xnew > 2 | xnew < 0.5);
        x(iinterpp) = x0(iinterpp);
    end
    fprintf('current interpp is %3.2e \n', exp(x(iinterpp)));
end

if (par.opt_alpha == on)
    ialpha = par.pindx.lalpha;
    xnew = exp(x(ialpha));
    xold = exp(x0(ialpha));
    if (xnew > 2*xold | xnew < 0.5*xold);
        x(ialpha) = x0(ialpha);
    end
    fprintf('current alpha is %3.2e \n', exp(x(ialpha)));
end

if (par.opt_beta == on)
    ibeta = par.pindx.lbeta;
    xnew = exp(x(ibeta));
    xold = exp(x0(ibeta));
    if (xnew > 2*xold | xnew < 0.5*xold);
        x(ibeta) = x0(ibeta);
    end
    fprintf('current beta is %3.2e \n', exp(x(ibeta)));
end

if par.Cmodel == on 
    if (par.opt_slopec == on)
        islopec = par.pindx.slopec;
        xnew = x(islopec);
        xold = x0(islopec);
        if (xnew > 1e-2 | xnew < -1e-2);
            x(islopec) = x0(islopec);
        end
        fprintf('current slopec is %3.2e \n', x(islopec));
    end
    
    if (par.opt_interpc == on)
        iinterpc = par.pindx.linterpc;
        xnew = exp(x(iinterpc));
        xold = exp(x0(iinterpc));
        if (xnew > 2 | xnew < 0.5);
            x(iinterpc) = x0(iinterpc);
        end
        fprintf('current interpc is %3.2e \n', exp(x(iinterpc)));
    end
    
    if (par.opt_d == on)
        id = par.pindx.ld;
        xnew = exp(x(id));
        xold = exp(x0(id));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(id) = x0(id);
        end
        fprintf('current d is %3.2e \n', exp(x(id)));
    end
    
    if (par.opt_kappa_dc == on)
        ikappa_dc = par.pindx.lkappa_dc;
        xnew = exp(x(ikappa_dc));
        xold = exp(x0(ikappa_dc));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(ikappa_dc) = x0(ikappa_dc);
        end
        fprintf('current kappa_dc is %3.2e \n', exp(x(ikappa_dc)));
    end
    
    if (par.opt_RR == on)
        iRR = par.pindx.lRR;
        xnew = exp(x(iRR));
        xold = exp(x0(iRR));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(iRR) = x0(iRR);
        end
        fprintf('current RR is %3.2e \n', exp(x(iRR)));
    end
    
    if (par.opt_cc == on)
        icc = par.pindx.lcc;
        xnew = exp(x(icc));
        xold = exp(x0(icc));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(icc) = x0(icc);
        end
        fprintf('current cc is %3.2e \n', exp(x(icc)));
    end
    
    if (par.opt_dd == on)
        idd = par.pindx.ldd;
        xnew = exp(x(idd));
        xold = exp(x0(idd));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(idd) = x0(idd);
        end
        fprintf('current dd is %3.2e \n', exp(x(idd)));
    end
end 
% ------------------------------------------------------------
if (par.Omodel == on)
    if (par.opt_slopeo == on)
        islopeo = par.pindx.slopeo;
        xnew = abs(x(islopeo));
        xold = abs(x0(islopeo));
        if (xnew > 5 | xnew < -5);
            x(islopeo) = x0(islopeo);
        end
        fprintf('current slopeo is %3.2e \n', x(islopeo));
    end
    
    if (par.opt_interpo == on)
        iinterpo = par.pindx.linterpo;
        xnew = exp(x(iinterpo));
        xold = exp(x0(iinterpo));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(iinterpo) = x0(iinterpo);
        end
        fprintf('current interpo is %3.2e \n', exp(x(iinterpo)));
    end
end
% ------------------------------------------------------------
% bsi
if (par.Simodel==on)
    if (par.opt_bsi == on)
        ibsi = par.pindx.lbsi;
        xnew = exp(x(ibsi));
        xold = exp(x0(ibsi));
        if (xnew > 1.5 | xnew < 0.2);
            x(ibsi) = x0(ibsi);
        end
        fprintf('current interpo is %3.2e \n', exp(x(ibsi)));
    end
    
    % at
    if (par.opt_at == on)
        iat = par.pindx.lat;
        xnew = exp(x(iat));
        xold = exp(x0(iat));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(iat) = x0(iat);
        end
        fprintf('current interpo is %3.2e \n', exp(x(iat)));
    end
    
    % bt
    if (par.opt_bt == on)
        ibt = par.pindx.lbt;
        xnew = exp(x(ibt));
        xold = exp(x0(ibt));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(ibt) = x0(ibt);
        end
        fprintf('current bt is %3.2e \n', exp(x(ibt)));
    end
    
    % aa
    if (par.opt_aa == on)
        iaa = par.pindx.aa;
        xnew = x(iaa);
        xold = x0(iaa);
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(iaa) = x0(iaa);
        end
        fprintf('current aa is %3.2e \n', x(iaa));
    end
    
    % bb
    if (par.opt_bb == on)
        ibb = par.pindx.lbb;
        xnew = exp(x(ibb));
        xold = exp(x0(ibb));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(ibb) = x0(ibb);
        end
        fprintf('current bb is %3.2e \n', exp(x(ibb)));
    end
    
    % kappa_gs
    if (par.opt_kappa_gs == on)
        ikappa_gs = par.pindx.lkappa_gs;
        xnew = exp(x(ikappa_gs));
        xold = exp(x0(ikappa_gs));
        if (xnew > 10*xold | xnew < 0.1*xold);
            x(ikappa_gs) = x0(ikappa_gs);
        end
        fprintf('current kappa_gs is %3.2e \n', exp(x(ikappa_gs)));
    end
end

% reset parm.p0;
parm.p0 = x;
%+++++++++restore the parameter values back to their original ones.
