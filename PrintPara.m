function PrintPara(x, par);
%
on = true; off = false;
%++++++++++ print out parameters in the log file
if (par.opt_sigma == on)
    isigma = par.pindx.lsigma;
    fprintf('current sigma is %3.2e \n', exp(x(isigma)));
end

if (par.opt_kappa_dp == on)
    ikappa_dp = par.pindx.lkappa_dp;
    fprintf('current kappa_dp is %3.2e \n', exp(x(ikappa_dp)));
end

if (par.opt_slopep == on)
    islopep = par.pindx.slopep;
    fprintf('current slopep is %3.2e \n', x(islopep));
end

if (par.opt_interpp == on)
    iinterpp = par.pindx.linterpp;
    fprintf('current interpp is %3.2e \n', exp(x(iinterpp)));
end

if (par.opt_alpha == on)
    ialpha = par.pindx.lalpha;
    fprintf('current alpha is %3.2e \n', exp(x(ialpha)));
end

if (par.opt_beta == on)
    ibeta = par.pindx.lbeta;
    fprintf('current beta is %3.2e \n', exp(x(ibeta)));
end

if par.Cmodel == on 
    if (par.opt_slopec == on)
        islopec = par.pindx.slopec;
        fprintf('current slopec is %3.2e \n', x(islopec));
    end
    
    if (par.opt_interpc == on)
        iinterpc = par.pindx.linterpc;
        fprintf('current interpc is %3.2e \n', exp(x(iinterpc)));
    end
    
    if (par.opt_d == on)
        id = par.pindx.ld;
        fprintf('current d is %3.2e \n', exp(x(id)));
    end
    
    if (par.opt_kappa_dc == on)
        ikappa_dc = par.pindx.lkappa_dc;
        fprintf('current kappa_dc is %3.2e \n', exp(x(ikappa_dc)));
    end
    
    if (par.opt_RR == on)
        iRR = par.pindx.lRR;
        fprintf('current RR is %3.2e \n', exp(x(iRR)));
    end
    
    if (par.opt_cc == on)
        icc = par.pindx.lcc;
        fprintf('current cc is %3.2e \n', exp(x(icc)));
    end
    
    if (par.opt_dd == on)
        idd = par.pindx.ldd;
        fprintf('current dd is %3.2e \n', exp(x(idd)));
    end
end 
% ------------------------------------------------------------
if (par.Omodel == on)
    if (par.opt_slopeo == on)
        islopeo = par.pindx.slopeo;
        fprintf('current slopeo is %3.2e \n', x(islopeo));
    end
    
    if (par.opt_interpo == on)
        iinterpo = par.pindx.linterpo;
        fprintf('current interpo is %3.2e \n', exp(x(iinterpo)));
    end
end
% ------------------------------------------------------------
% bsi
if (par.Simodel==on)
    if (par.opt_bsi == on)
        ibsi = par.pindx.lbsi;
        fprintf('current bsi is %3.2e \n', exp(x(ibsi)));
    end
    
    % at
    if (par.opt_at == on)
        iat = par.pindx.lat;
        fprintf('current at is %3.2e \n', exp(x(iat)));
    end
    
    % bt
    if (par.opt_bt == on)
        ibt = par.pindx.lbt;
        fprintf('current bt is %3.2e \n', exp(x(ibt)));
    end
    
    % aa
    if (par.opt_aa == on)
        iaa = par.pindx.aa;
        fprintf('current aa is %3.2e \n', x(iaa));
    end
    
    % bb
    if (par.opt_bb == on)
        ibb = par.pindx.lbb;
        fprintf('current bb is %3.2e \n', exp(x(ibb)));
    end
    
end

