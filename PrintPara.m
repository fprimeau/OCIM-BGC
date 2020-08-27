function PrintPara(x, par);
    on = true; off = false;
    %++++++++++ print out parameters to the log file
    if (par.opt_sigma == on)
        isigma = par.pindx.lsigma;
        fprintf('current sigma is    % 3.2e \n', exp(x(isigma)));
        xhat.sigma = exp(x(isigma));
    end

    if (par.opt_kappa_dp == on)
        ikappa_dp = par.pindx.lkappa_dp;
        fprintf('current kappa_dp is % 3.2e \n', exp(x(ikappa_dp)));
        xhat.kappa_dp = exp(x(ikappa_dp));
    end

    if (par.opt_bP_T == on)
        ibP_T = par.pindx.bP_T;
        fprintf('current bP_T is     % 3.2e \n', x(ibP_T));
        xhat.bP_T = x(ibP_T);
    end

    if (par.opt_bP == on)
        ibP = par.pindx.lbP;
        fprintf('current bP is       % 3.2e \n', exp(x(ibP)));
        xhat.bP = exp(x(ibP));
    end

    if (par.opt_alpha == on)
        ialpha = par.pindx.lalpha;
        fprintf('current alpha is    % 3.2e \n', exp(x(ialpha)));
        xhat.alpha = exp(x(ialpha));
    end

    if (par.opt_beta == on)
        ibeta = par.pindx.lbeta;
        fprintf('current beta is     % 3.2e \n', exp(x(ibeta)));
        xhat.beta = exp(x(ibeta));
    end

    if par.Cmodel == on 
        if (par.opt_bC_T == on)
            ibC_T = par.pindx.bC_T;
            fprintf('current bC_T is     % 3.2e \n', x(ibC_T));
            xhat.bC_T = x(ibC_T);
        end
        
        if (par.opt_bC == on)
            ibC = par.pindx.lbC;
            fprintf('current bC is       % 3.2e \n', exp(x(ibC)));
            xhat.bC = exp(x(ibC));
        end
        
        if (par.opt_d == on)
            id = par.pindx.ld;
            fprintf('current d is        % 3.2e \n', exp(x(id)));
            xhat.d = exp(x(id));
        end
        
        if (par.opt_kappa_dc == on)
            ikappa_dc = par.pindx.lkappa_dc;
            fprintf('current kappa_dc is % 3.2e \n', exp(x(ikappa_dc)));
            xhat.kappa_dc = exp(x(ikappa_dc));
        end
        
        if (par.opt_RR == on)
            iRR = par.pindx.lRR;
            fprintf('current RR is       % 3.2e \n', exp(x(iRR)));
            xhat.RR = exp(x(iRR));
        end
        
        if (par.opt_cc == on)
            icc = par.pindx.lcc;
            fprintf('current cc is       % 3.2e \n', exp(x(icc)));
            xhat.cc = exp(x(icc));
        end
        
        if (par.opt_dd == on)
            idd = par.pindx.ldd;
            fprintf('current dd is       % 3.2e \n', exp(x(idd)));
            xhat.dd = exp(x(idd));
        end
    end 
    % ------------------------------------------------------------
    if (par.Omodel == on)
        if (par.opt_slopeo == on)
            islopeo = par.pindx.slopeo;
            fprintf('current slopeo is   % 3.2e \n', x(islopeo));
            xhat.slopeo = x(islopeo);
        end
        
        if (par.opt_interpo == on)
            iinterpo = par.pindx.linterpo;
            fprintf('current interpo is  % 3.2e \n', exp(x(iinterpo)));
            xhat.interpo = exp(x(iinterpo));
        end
    end
    % ------------------------------------------------------------
    % dsi
    if (par.Simodel==on)
        if (par.opt_dsi == on)
            idsi = par.pindx.ldsi;
            fprintf('current dsi is      % 3.2e \n', exp(x(idsi)));
            xhat.dsi = exp(x(idsi));
        end
        
        % at
        if (par.opt_at == on)
            iat = par.pindx.lat;
            fprintf('current at is       % 3.2e \n', exp(x(iat)));
            xhat.at = exp(x(iat));
        end
        
        % bt
        if (par.opt_bt == on)
            ibt = par.pindx.lbt;
            fprintf('current bt is       % 3.2e \n', exp(x(ibt)));
            xhat.bt = exp(x(ibt));
        end
        
        % aa
        if (par.opt_aa == on)
            iaa = par.pindx.aa;
            fprintf('current aa is       % 3.2e \n', x(iaa));
            xhat.aa = x(iaa);
        end
        
        % bb
        if (par.opt_bb == on)
            ibb = par.pindx.lbb;
            fprintf('current bb is       % 3.2e \n', exp(x(ibb)));
            xhat.bb = exp(x(ibb));
        end
    end
    %
    if (par.optim == on)
        save(par.fxhat, 'xhat')
    end 
end

