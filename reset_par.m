function  [par, x] = reset_par(x, par);
    on = true; off = false;
    %++++++++++ check if the optimization routine suggests strange
    % parameter values
    A = isfield(par,'x0');
    if A == 0
        x0 = par.p0;
    end
    
    if (par.opt_sigma == on)
        isigma = par.pindx.lsigma;
        xnew = exp(x(isigma));
        xold = exp(x0(isigma));
        if (xnew > 1 | xnew < 0)
            x(isigma) = x0(isigma);
        end
    end

    if (par.opt_kappa_dp == on)
        ikappa_dp = par.pindx.lkappa_dp;
        xnew = exp(x(ikappa_dp));
        xold = exp(x0(ikappa_dp));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(ikappa_dp) = x0(ikappa_dp);
        end
    end

    if (par.opt_bP_T == on)
        ibP_T = par.pindx.bP_T;
        xnew = x(ibP_T);
        xold = x0(ibP_T);
        if (xnew > 5e-2 | xnew < -5e-2);
            x(ibP_T) = x0(ibP_T);
        end
    end

    if (par.opt_bP == on)
        ibP = par.pindx.lbP;
        xnew = exp(x(ibP));
        xold = exp(x0(ibP));
        if (xnew > 5*xold | xnew < 0.2*xold);
            x(ibP) = x0(ibP);
        end
    end

    if (par.opt_alpha == on)
        ialpha = par.pindx.lalpha;
        xnew = exp(x(ialpha));
        xold = exp(x0(ialpha));
        if (xnew > 2*xold | xnew < 0.5*xold);
            x(ialpha) = x0(ialpha);
        end
    end

    if (par.opt_beta == on)
        ibeta = par.pindx.lbeta;
        xnew = exp(x(ibeta));
        xold = exp(x0(ibeta));
        if (xnew > 2*xold | xnew < 0.5*xold);
            x(ibeta) = x0(ibeta);
        end
    end

    if par.Cmodel == on 
        if (par.opt_bC_T == on)
            ibC_T = par.pindx.bC_T;
            xnew = x(ibC_T);
            xold = x0(ibC_T);
            if (xnew > 5e-2 | xnew < -5e-2);
                x(ibC_T) = x0(ibC_T);
            end
        end
        
        if (par.opt_bC == on)
            ibC = par.pindx.lbC;
            xnew = exp(x(ibC));
            xold = exp(x0(ibC));
            if (xnew > 10*xold | xnew < 0.1*xold);
                x(ibC) = x0(ibC);
            end
        end
        
        if (par.opt_d == on)
            id = par.pindx.ld;
            xnew = exp(x(id));
            xold = exp(x0(id));
            if (xnew > 5*xold | xnew < 0.2*xold);
                x(id) = x0(id);
            end
        end
        
        if (par.opt_kappa_dc == on)
            ikappa_dc = par.pindx.lkappa_dc;
            xnew = exp(x(ikappa_dc));
            xold = exp(x0(ikappa_dc));
            if (xnew > 5*xold | xnew < 0.2*xold);
                x(ikappa_dc) = x0(ikappa_dc);
            end
        end
        
        if (par.opt_RR == on)
            iRR = par.pindx.lRR;
            xnew = exp(x(iRR));
            xold = exp(x0(iRR));
            if (xnew > 5*xold | xnew < 0.2*xold);
                x(iRR) = x0(iRR);
            end
        end
        
        if (par.opt_cc == on)
            icc = par.pindx.lcc;
            xnew = exp(x(icc));
            xold = exp(x0(icc));
            if (xnew > 5*xold | xnew < 0.2*xold);
                x(icc) = x0(icc);
            end
        end
        
        if (par.opt_dd == on)
            idd = par.pindx.ldd;
            xnew = exp(x(idd));
            xold = exp(x0(idd));
            if (xnew > 5*xold | xnew < 0.2*xold);
                x(idd) = x0(idd);
            end
        end
    end 
    % ------------------------------------------------------------
    if (par.Omodel == on)
        if (par.opt_slopeo == on)
            islopeo = par.pindx.slopeo;
            xnew = abs(x(islopeo));
            xold = abs(x0(islopeo));
            if (xnew > 50 | xnew < -50);
                x(islopeo) = x0(islopeo);
            end
        end
        
        if (par.opt_interpo == on)
            iinterpo = par.pindx.linterpo;
            xnew = exp(x(iinterpo));
            xold = exp(x0(iinterpo));
            if (xnew > 5*xold | xnew < 0.2*xold);
                x(iinterpo) = x0(iinterpo);
            end
        end
    end
    % ------------------------------------------------------------
    % dsi
    if (par.Simodel==on)
        if (par.opt_dsi == on)
            idsi = par.pindx.ldsi;
            xnew = exp(x(idsi));
            xold = exp(x0(idsi));
            if (xnew > 2*xold | xnew < 0.5*xold);
                x(idsi) = x0(idsi);
            end
        end
        % at
        if (par.opt_at == on)
            iat = par.pindx.lat;
            xnew = exp(x(iat));
            xold = exp(x0(iat));
            if (xnew > 5*xold | xnew < 0.2*xold);
                x(iat) = x0(iat);
            end
        end
        % bt
        if (par.opt_bt == on)
            ibt = par.pindx.lbt;
            xnew = exp(x(ibt));
            xold = exp(x0(ibt));
            if (xnew > 15000 | xnew < 10000);
                x(ibt) = x0(ibt);
            end
        end
        % aa
        if (par.opt_aa == on)
            iaa = par.pindx.aa;
            xnew = x(iaa);
            xold = x0(iaa);
            if (xnew > 50 | xnew < -50);
                x(iaa) = x0(iaa);
            end
        end
        % bb
        if (par.opt_bb == on)
            ibb = par.pindx.lbb;
            xnew = exp(x(ibb));
            xold = exp(x0(ibb));
            if (xnew > 5*xold | xnew < 0.2*xold);
                x(ibb) = x0(ibb);
            end
        end    
    end
    % reset par.p0;
    par.p0 = x;
    %+++++++++restore the parameter values back to their original ones.
end

