function par = PackPar(par)
    on = true ; off = false ;
    npx = 0; ncx = 0;
    nox = 0; nsx = 0;
	nbx = 0; %nbx for cell model (BIO)
    p0 = [];
    % sigma
    if (par.opt_sigma == on)
        npx    = npx + 1;
        sigma  = par.sigma      ;
        lsigma = log(sigma)     ;
        strt   = length(p0) + 1 ;
        p0     = [p0; lsigma]   ;
        par.pindx.lsigma = strt : length(p0);
    end
    % kP_T
    if (par.opt_kP_T == on)
        npx  = npx + 1        ;
        kP_T = par.kP_T       ;
        strt = length(p0) + 1 ;
        p0   = [p0; kP_T]     ;
        par.pindx.kP_T = strt : length(p0);
    end
    % kdP
    if (par.opt_kdP == on)
        npx  = npx + 1        ;
        kdP  = par.kdP        ;
        lkdP = log(kdP)       ;
        strt = length(p0) + 1 ;
        p0   = [p0; lkdP]     ;
        par.pindx.lkdP = strt : length(p0);
    end
    % bP_T
    if (par.opt_bP_T == on)
        npx  = npx + 1        ;
        bP_T = par.bP_T       ;
        strt = length(p0) + 1 ;
        p0   = [p0; bP_T]     ;
        par.pindx.bP_T = strt : length(p0);
    end
    % bP
    if (par.opt_bP == on)
        npx  = npx + 1        ;
        bP   = par.bP         ;
        lbP  = log(bP)        ;
        strt = length(p0) + 1 ;
        p0   = [p0; lbP]      ;
        par.pindx.lbP = strt  : length(p0);
    end
    % alpha
    if (par.opt_alpha == on)
        npx    = npx + 1        ;
        alpha  = par.alpha      ;
        lalpha = log(alpha)     ;
        strt   = length(p0) + 1 ;
        p0     = [p0; lalpha]   ;
        par.pindx.lalpha = strt : length(p0);
    end
    % beta
    if (par.opt_beta == on)
        npx   = npx + 1        ;
        beta  = par.beta       ;
        lbeta = log(beta)      ;
        strt  = length(p0) + 1 ;
        p0    = [p0; lbeta]    ;
        par.pindx.lbeta = strt : length(p0);
    end
    %
    if (par.Cmodel == on)
        % bC_T
        if (par.opt_bC_T == on)
            ncx  = ncx + 1        ;
            bC_T = par.bC_T       ;
            strt = length(p0) + 1 ;
            p0   = [p0; bC_T]     ;
            par.pindx.bC_T = strt : length(p0);
        end
        % bC
        if (par.opt_bC == on)
            ncx  = ncx + 1        ;
            bC   = par.bC         ;
            lbC  = log(bC)        ;
            strt = length(p0) + 1 ;
            p0   = [p0; lbC]      ;
            par.pindx.lbC = strt  : length(p0);
        end
        % d
        if (par.opt_d == on)
            ncx  = ncx + 1        ;
            d    = par.d          ;
            ld   = log(d)         ;
            strt = length(p0) + 1 ;
            p0   = [p0; ld]       ;
            par.pindx.ld = strt   : length(p0);
        end
        % kC_T
        if (par.opt_kC_T == on)
            ncx  = ncx + 1        ;
            kC_T = par.kC_T       ;
            strt = length(p0) + 1 ;
            p0   = [p0; kC_T]     ;
            par.pindx.kC_T = strt : length(p0);
        end
        % kdC
        if (par.opt_kdC == on)
            ncx  = ncx + 1        ;
            kdC  = par.kdC        ;
            lkdC = log(kdC)       ;
            strt = length(p0) + 1 ;
            p0   = [p0; lkdC]     ;
            par.pindx.lkdC = strt : length(p0);
        end
        % R_Si
        if (par.opt_R_Si == on)
            ncx   = ncx + 1        ;
            R_Si  = par.R_Si       ;
            strt  = length(p0) + 1 ;
            p0    = [p0; R_Si]    ;
            par.pindx.R_Si = strt : length(p0);
        end
        % rR
        if (par.opt_rR == on)
            ncx  = ncx + 1        ;
            rR   = par.rR         ;
            lrR  = log(rR)        ;
            strt = length(p0) + 1 ;
            p0   = [p0; lrR]      ;
            par.pindx.lrR = strt  : length(p0);
        end
        % cc
        if (par.opt_cc == on)
            ncx  = ncx + 1        ;
            cc   = par.cc         ;
            lcc  = log(cc)        ;
            strt = length(p0) + 1 ;
            p0   = [p0; lcc]      ;
            par.pindx.lcc = strt  : length(p0);
        end
        % dd
        if (par.opt_dd == on)
            ncx  = ncx + 1        ;
            dd   = par.dd         ;
            ldd  = log(dd)        ;
            strt = length(p0) + 1 ;
            p0   = [p0; ldd]      ;
            par.pindx.ldd = strt  : length(p0);
        end
        par.ncx = ncx;
    end
    %
    if par.Omodel == on
        % O2C_T
        if (par.opt_O2C_T == on)
            nox   = nox + 1        ;
            O2C_T = par.O2C_T      ;
            strt  = length(p0)+1   ;
            p0    = [p0; O2C_T]    ;
            par.pindx.O2C_T = strt : length(p0);
        end
        % rO2C
        if (par.opt_rO2C == on)
            nox   = nox + 1        ;
            rO2C  = par.rO2C       ;
            lrO2C = log(rO2C)      ;
            strt  = length(p0) + 1 ;
            p0    = [p0; lrO2C]    ;
            par.pindx.lrO2C = strt : length(p0);
        end
        % O2P_T
        if (par.opt_O2P_T == on)
            nox   = nox + 1        ;
            O2P_T = par.O2P_T      ;
            strt  = length(p0) + 1 ;
            p0    = [p0; O2P_T]    ;
            par.pindx.O2P_T = strt : length(p0);
        end
        % rO2P
        if (par.opt_rO2P == on)
            nox   = nox + 1        ;
            rO2P  = par.rO2P       ;
            lrO2P = log(rO2P)      ;
            strt  = length(p0) + 1 ;
            p0    = [p0; lrO2P]    ;
            par.pindx.lrO2P = strt : length(p0);
        end
    end
    if par.Simodel == on
        % dsi
        if (par.opt_dsi == on)
            nsx  = nsx + 1        ;
            dsi  = par.dsi        ;
            ldsi = log(dsi)       ;
            strt = length(p0) + 1 ;
            p0   = [p0; ldsi]     ;
            par.pindx.ldsi = strt : length(p0);
        end
        % at
        if (par.opt_at == on)
            nsx  = nsx + 1        ;
            at   = par.at         ;
            lat  = log(at)        ;
            strt = length(p0) + 1 ;
            p0   = [p0; lat]      ;
            par.pindx.lat = strt  : length(p0);
        end
        % bt
        if (par.opt_bt == on)
            nsx  = nsx + 1        ;
            bt   = par.bt         ;
            lbt  = log(bt)        ;
            strt = length(p0) + 1 ;
            p0   = [p0; lbt]      ;
            par.pindx.lbt = strt  : length(p0);
        end
        % aa
        if (par.opt_aa == on)
            nsx  = nsx + 1        ;
            aa   = par.aa         ;
            strt = length(p0) + 1 ;
            p0   = [p0; aa]       ;
            par.pindx.aa = strt   : length(p0);
        end
        % bb
        if (par.opt_bb == on)
            nsx  = nsx + 1        ;
            bb   = par.bb         ;
            lbb  = log(bb)        ;
            strt = length(p0) + 1 ;
            p0   = [p0; lbb]      ;
            par.pindx.lbb = strt  : length(p0);
        end
    end
	if par.Cellmodel == on
		% Q10Photo
		if (par.opt_Q10Photo == on)
			nbx  = nbx + 1        		;
			Q10Photo = par.BIO.Q10Photo ;
			lQ10Photo = log(par.BIO.Q10Photo);
			strt = length(p0) + 1		;
			p0 = [p0; lQ10Photo]		;
			par.pindx.lQ10Photo = strt : length(p0);
		end
		if (par.opt_fStorage == on)
			nbx  = nbx + 1        		;
			lfStorage =log(par.BIO.fStorage);
			strt = length(p0) + 1		;
			p0 = [p0; lfStorage]		;
			par.pindx.lfStorage = strt : length(p0);
		end
		if (par.opt_fRibE == on)
			nbx  = nbx + 1        		;
			%lfRibE = log(par.BIO.fRibE) ;
			tfRibE = atanh(2*par.BIO.fRibE-1);
			%fRibE = 0.5*(1+tanh(x(par.pindx.tfRibE))); %change to tanh so fRibE must be between 0 and 1
			strt = length(p0) + 1		;
			p0 = [p0; tfRibE]			;
			par.pindx.tfRibE = strt : length(p0);
		end
		if (par.opt_kST0 == on);
			nbx  = nbx + 1        		;
			lkST0 = log(par.BIO.kST0)	;
			strt = length(p0) + 1		;
			p0 = [p0; lkST0]			;
			par.pindx.lkST0 = strt : length(p0);
		end
		if (par.opt_PLip_PCutoff == on);
			nbx  = nbx + 1        		;
			lPLip_PCutoff = log(par.BIO.PLip_PCutoff)	;
			strt = length(p0) + 1		;
			p0 = [p0; lPLip_PCutoff]			;
			par.pindx.lPLip_PCutoff = strt : length(p0);
		end
		if (par.opt_PLip_scale == on);
			nbx  = nbx + 1        		;
			lPLip_scale = log(par.BIO.PLip_scale)	;
			strt = length(p0) + 1		;
			p0 = [p0; lPLip_scale]			;
			par.pindx.lPLip_scale = strt : length(p0);
		end
		if (par.opt_PStor_rCutoff == on);
			nbx  = nbx + 1        		;
			lPStor_rCutoff = log(par.BIO.PStor_rCutoff)	;
			strt = length(p0) + 1		;
			p0 = [p0; lPStor_rCutoff]			;
			par.pindx.lPStor_rCutoff = strt : length(p0);
		end
		if (par.opt_PStor_scale == on);
			nbx  = nbx + 1        		;
			lPStor_scale = log(par.BIO.PStor_scale)	;
			strt = length(p0) + 1		;
			p0 = [p0; lPStor_scale]			;
			par.pindx.lPStor_scale = strt : length(p0);
		end
		if (par.opt_alphaS == on)
			nbx  = nbx + 1        		;
			lalphaS = log(par.BIO.alphaS)	;
			strt = length(p0) + 1		;
			p0 = [p0; lalphaS]			;
			par.pindx.lalphaS = strt : length(p0);
		end
		%if isfield(par,'opt_gammaDNA')
			if (par.opt_gammaDNA == on)
				nbx  = nbx + 1        		;
				%lfRibE = log(par.BIO.fRibE) ;
				tgammaDNA = atanh(2*par.BIO.gammaDNA-1);
				%gammaDNA = 0.5*(1+tanh(x(par.pindx.tgammaDNA))); %change to tanh so fRibE must be between 0 and 1
				strt = length(p0) + 1		;
				p0 = [p0; tgammaDNA]			;
				par.pindx.tgammaDNA = strt : length(p0);
			end
		%end
	end
    par.p0  = p0  ;
    par.npx = npx ; par.ncx = ncx ;
    par.nox = nox ; par.nsx = nsx ;
	par.nbx = nbx ;
end
