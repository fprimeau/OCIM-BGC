function  x = reset_par(x, par);
    global iter
    on = true; off = false;
    load(par.fxhat,'x0')
    % check if the optimization routine suggests strange
    fb = 2 ;  fs = 0.5 ;
    % gradually decrease the constrains    
    fb = fb * 1.1^iter ;
    fs = fs * 0.9^iter ; 
    % parameter values
    if (par.opt_sigma == on)
        isigma = par.pindx.lsigma;
        xnew   = exp(x(isigma));
        xold   = exp(x0(isigma));
        if (xnew > 1 | xnew < 0)
            x(isigma) = x0(isigma);
        end
    end

    if (par.opt_kP_T == on & par.opt_kdP == on)
        ikP_T    = par.pindx.kP_T;
        kP_T_new = x(ikP_T);
        kP_T_old = x0(ikP_T);
        %
        ikdP    = par.pindx.lkdP;
        kdP_new = exp(x(ikdP));
        kdP_old = exp(x0(ikdP));
        %
        kP = kP_T_new*par.Tz + kdP_new ;
        if any(kP < 0)
            x(ikP_T) = x0(ikP_T);
            % x(ikdP) = x0(ikdP);
        end
    end
    
    if (par.opt_kP_T == off &par.opt_kdP == on)
        ikdP = par.pindx.lkdP;
        xnew = exp(x(ikdP));
        xold = exp(x0(ikdP));
        if (xnew > fb*xold | xnew < fs*xold);
            x(ikdP) = x0(ikdP);
        end
    end

    if (par.opt_bP_T == on & par.opt_bP == on)
        ibP_T = par.pindx.bP_T;
        bm1   = x(ibP_T);
        bm0   = x0(ibP_T);
        %
        ibP  = par.pindx.lbP;
        bb1  = exp(x(ibP));
        bb0  = exp(x0(ibP));
        %
        bP   = bm1*par.aveT + bb1 ;
        if any(bP(:) < 0)
            x(ibP_T) = x0(ibP_T);
            x(ibP) = x0(ibP);
        end
    end

    if (par.opt_bP_T == off & par.opt_bP == on)
        ibP  = par.pindx.lbP;
        xnew = exp(x(ibP));
        xold = exp(x0(ibP));
        if (xnew > fb*xold | xnew < fs*xold);
            x(ibP) = x0(ibP);
        end
    end

    if (par.opt_alpha == on)
        ialpha = par.pindx.lalpha;
        xnew   = exp(x(ialpha));
        xold   = exp(x0(ialpha));
        if (xnew > fb*xold | xnew < fs*xold);
            x(ialpha) = x0(ialpha);
        end
    end

    if (par.opt_beta == on)
        ibeta = par.pindx.lbeta;
        xnew  = exp(x(ibeta));
        xold  = exp(x0(ibeta));
        if (xnew > fb*xold | xnew < fs*xold);
            x(ibeta) = x0(ibeta);
        end
    end

    if (par.Cmodel == on) 
        if (par.opt_bC_T == on & par.opt_bC == on)
            ibC_T = par.pindx.bC_T;
            bm1   = x(ibC_T);
            bm0   = x0(ibC_T);
            %
            ibC  = par.pindx.lbC;
            bb1  = exp(x(ibC));
            bb0  = exp(x0(ibC));
            %
            bC   = bm1*par.aveT + bb1 ;
            if any(bC(:) < 0)
                x(ibC_T) = x0(ibC_T);
                x(ibC)   = x0(ibC);
            end
        end
        
        if (par.opt_bC_T == off & par.opt_bC == on)
            ibC  = par.pindx.lbC;
            xnew = exp(x(ibC));
            xold = exp(x0(ibC));
            if (xnew > 3 | xnew < 0.3);
                x(ibC) = x0(ibC);
            end
        end
        
        if (par.opt_d == on)
            id   = par.pindx.ld;
            xnew = exp(x(id));
            xold = exp(x0(id));
            if (xnew > fb*xold | xnew < fs*xold);
                x(id) = x0(id);
            end
        end

        if (par.opt_kC_T == on & par.opt_kdC == on)
            ikC_T    = par.pindx.kC_T;
            kC_T_new = x(ikC_T);
            kC_T_old = x0(ikC_T);
            %
            ikdC     = par.pindx.lkdC;
            kdC_new  = exp(x(ikdC));
            kdC_old  = exp(x0(ikdC));
            %
            kC = kC_T_new*par.Tz + kdC_new ;
            if any(kC < 0)
                x(ikC_T) = x0(ikC_T);
            end
        end

        if (par.opt_kC_T == off & par.opt_kdC == on)
            ikdC = par.pindx.lkdC;
            xnew = exp(x(ikdC));
            xold = exp(x0(ikdC));
            if (xnew > fb*xold | xnew < fs*xold);
                x(ikdC) = x0(ikdC);
            end
        end
        
        if (par.opt_RR == on)
            iRR  = par.pindx.lRR;
            xnew = exp(x(iRR));
            xold = exp(x0(iRR));
            if (xnew > fb*xold | xnew < fs*xold);
                x(iRR) = x0(iRR);
            end
        end
        
        if (par.opt_cc == on)
            icc  = par.pindx.lcc;
            xnew = exp(x(icc));
            xold = exp(x0(icc));
            if (xnew > fb*xold | xnew < fs*xold);
                x(icc) = x0(icc);
            end
        end
        
        if (par.opt_dd == on)
            idd  = par.pindx.ldd;
            xnew = exp(x(idd));
            xold = exp(x0(idd));
            if (xnew > fb*xold | xnew < fs*xold);
                x(idd) = x0(idd);
            end
        end
    end 
    % ------------------------------------------------------------
    if (par.Omodel == on)
        if (par.opt_slopeo == on)
            islopeo = par.pindx.slopeo;
            xnew    = abs(x(islopeo));
            xold    = abs(x0(islopeo));
            if (xnew > 50 | xnew < -50);
                x(islopeo) = x0(islopeo);
            end
        end
        
        if (par.opt_interpo == on)
            iinterpo = par.pindx.linterpo;
            xnew     = exp(x(iinterpo));
            xold     = exp(x0(iinterpo));
            if (xnew > fb*xold | xnew < fs*xold);
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
            if (xnew > 5000 | xnew < 1000);
                x(idsi) = x0(idsi);
            end
        end
        % at
        if (par.opt_at == on)
            iat  = par.pindx.lat;
            xnew = exp(x(iat));
            xold = exp(x0(iat));
            if (xnew > fb*xold | xnew < fs*xold);
                x(iat) = x0(iat);
            end
        end
        % bt
        if (par.opt_bt == on)
            ibt  = par.pindx.lbt;
            xnew = exp(x(ibt));
            xold = exp(x0(ibt));
            if (xnew > 15000 | xnew < 10000);
                x(ibt) = x0(ibt);
            end
        end
        % aa
        if (par.opt_aa == on)
            iaa  = par.pindx.aa;
            xnew = x(iaa);
            xold = x0(iaa);
            if (xnew > 50 | xnew < -50);
                x(iaa) = x0(iaa);
            end
        end
        % bb
        if (par.opt_bb == on)
            ibb  = par.pindx.lbb;
            xnew = exp(x(ibb));
            xold = exp(x0(ibb));
            if (xnew > fb*xold | xnew < fs*xold);
                x(ibb) = x0(ibb);
            end
        end    
    end
end

