function PrintPara(x, par);
    global iter
    on = true; off = false;
    %++++++++++ print out all parameters to the log file
    if iter == 0
        fprintf('All parameters \n')
        fprintf('-----------------------------------------------\n')
        fprintf('kappa_p is  % 3.2e \n', par.kappa_p)   ;
        fprintf('kappa_g is  % 3.2e \n', par.kappa_g)   ;
        fprintf('sigma   is  % 3.2e \n', par.sigma)     ;
        fprintf('kP_T    is  % 3.2e \n', par.kP_T)      ;
        fprintf('kdP     is  % 3.2e \n', par.kdP)       ;
        fprintf('bP_T    is  % 3.2e \n', par.bP_T)      ;
        fprintf('bP      is  % 3.2e \n', par.bP)        ;
        fprintf('alpha   is  % 3.2e \n', par.alpha)     ;
        fprintf('beta    is  % 3.2e \n', par.beta)      ;
        if (par.Cmodel == on)
            fprintf('kPIC    is  % 3.2e \n', par.kPIC)  ;
            fprintf('bC_T    is  % 3.2e \n', par.bC_T)  ;
            fprintf('bC      is  % 3.2e \n', par.bC)    ;
            fprintf('d       is  % 3.2e \n', par.d)     ;
            fprintf('kC_T    is  % 3.2e \n', par.kC_T)  ;
            fprintf('kdC     is  % 3.2e \n', par.kdC)   ;
            fprintf('R_Si    is  % 3.2e \n', par.R_Si)  ;
            fprintf('rR      is  % 3.2e \n', par.rR)    ;
            fprintf('cc      is  % 3.2e \n', par.cc)    ;
            fprintf('dd      is  % 3.2e \n', par.dd)    ;
        end
        if (par.Omodel == on)
            fprintf('O2C_T   is  % 3.2e \n', par.O2C_T) ;
            fprintf('rO2C    is  % 3.2e \n', par.rO2C)  ;
            fprintf('O2P_T   is  % 3.2e \n', par.O2P_T) ;
            fprintf('rO2P    is  % 3.2e \n', par.rO2P)  ;
        end
        if (par.Simodel==on)
            fprintf('dsi     is  % 3.2e \n', par.dsi)   ;
            fprintf('at      is  % 3.2e \n', par.at)    ;
            fprintf('bt      is  % 3.2e \n', par.bt)    ;
            fprintf('aa      is  % 3.2e \n', par.iaa)   ;
            fprintf('bb      is  % 3.2e \n\n', par.bb)  ;
        end
		if (par.Cellmodel==on)
            fprintf('alphaS         is  % 3.2e \n', par.BIO.alphaS)   ;
            fprintf('gammDNA        is  % 3.2e \n', par.BIO.gammaDNA)    ;
            fprintf('gammaLipid     is  % 3.2e \n', par.BIO.gammaLipid)    ;
            fprintf('PCutoff        is  % 3.2e \n', exp(par.BIO.lPCutoff))   ;
            fprintf('r0Cutoff       is  % 3.2e \n', par.BIO.r0Cutoff)  ;
			fprintf('DNT0           is  % 3.2e \n', par.BIO.DNT0)  ;
			fprintf('DPT0           is  % 3.2e \n', par.BIO.DPT0)  ;
			fprintf('Q10Diffusivity is  % 3.2e \n', par.BIO.Q10Diffusivity)  ;
			fprintf('AMin           is  % 3.2e \n', par.BIO.AMin)  ;
			fprintf('CStor          is  % 3.2e \n', par.BIO.CStor)  ;
			fprintf('alphaPLip      is  % 3.2e \n', par.BIO.alphaPLip)  ;
			fprintf('PhiS           is  % 3.2e \n', par.BIO.PhiS)  ;
			fprintf('pDry           is  % 3.2e \n', par.BIO.pDry)  ;
			fprintf('rho      	    is  % 3.2e \n', par.BIO.rho)  ;
			fprintf('fProtM         is  % 3.2e \n', par.BIO.fProtM)  ;
			fprintf('fProtL         is  % 3.2e \n', par.BIO.fProtL)  ;
			fprintf('PDNA      	    is  % 3.2e \n', par.BIO.PDNA)  ;
			fprintf('PRib      	    is  % 3.2e \n', par.BIO.PRib)  ;
			fprintf('PPhospholipid  is  % 3.2e \n', par.BIO.PPhospholipid)  ;
			fprintf('NProt          is  % 3.2e \n', par.BIO.NProt)  ;
			fprintf('NDNA      	    is  % 3.2e \n', par.BIO.NDNA)  ;
			fprintf('NRib      	    is  % 3.2e \n', par.BIO.NRib)  ;
			fprintf('CProt          is  % 3.2e \n', par.BIO.CProt)  ;
			fprintf('CDNA      	    is  % 3.2e \n', par.BIO.CDNA)  ;
			fprintf('CPhospholipid  is  % 3.2e \n', par.BIO.CPhospholipid )  ;
			fprintf('CLipid         is  % 3.2e \n', par.BIO.CLipid)  ;
			fprintf('CRib      	    is  % 3.2e \n\n', par.BIO.CRib)  ;
        end
        fprintf('-----------------------------------------------\n\n')
    end
    %++++++++++ print out parameters to the log file
    fprintf('Tunable parameters \n')
    fprintf('-----------------------------------------------\n')
    if (par.opt_sigma == on)
        isigma = par.pindx.lsigma   ;
        fprintf('current sigma   is  % 3.2e \n', exp(x(isigma)));
        xhat.sigma = exp(x(isigma)) ;
    end

    if (par.opt_kP_T == on)
        ikP_T = par.pindx.kP_T  ;
        fprintf('current kP_T    is  % 3.2e \n', x(ikP_T));
        xhat.kP_T = x(ikP_T)    ;
    end

    if (par.opt_kdP == on)
        ikdP = par.pindx.lkdP   ;
        fprintf('current kdP     is  % 3.2e \n', exp(x(ikdP)));
        xhat.kdP = exp(x(ikdP)) ;
    end

    if (par.opt_bP_T == on)
        ibP_T = par.pindx.bP_T  ;
        fprintf('current bP_T    is  % 3.2e \n', x(ibP_T));
        xhat.bP_T = x(ibP_T)    ;
    end

    if (par.opt_bP == on)
        ibP = par.pindx.lbP     ;
        fprintf('current bP      is  % 3.2e \n', exp(x(ibP)));
        xhat.bP = exp(x(ibP))   ;
    end

    if (par.opt_alpha == on)
        ialpha = par.pindx.lalpha;
        fprintf('current alpha   is  % 3.2e \n', exp(x(ialpha)));
        xhat.alpha = exp(x(ialpha));
    end

    if (par.opt_beta == on)
        ibeta = par.pindx.lbeta;
        fprintf('current beta    is  % 3.2e \n', exp(x(ibeta)));
        xhat.beta = exp(x(ibeta));
    end

    if par.Cmodel == on
        if (par.opt_bC_T == on)
            ibC_T = par.pindx.bC_T;
            fprintf('current bC_T    is  % 3.2e \n', x(ibC_T));
            xhat.bC_T = x(ibC_T);
        end

        if (par.opt_bC == on)
            ibC = par.pindx.lbC;
            fprintf('current bC      is  % 3.2e \n', exp(x(ibC)));
            xhat.bC = exp(x(ibC));
        end

        if (par.opt_d == on)
            id = par.pindx.ld;
            fprintf('current d       is  % 3.2e \n', exp(x(id)));
            xhat.d = exp(x(id));
        end

        if (par.opt_kC_T == on)
            ikC_T = par.pindx.kC_T;
            fprintf('current kC_T    is  % 3.2e \n', x(ikC_T));
            xhat.kC_T = x(ikC_T);
        end

        if (par.opt_kdC == on)
            ikdC = par.pindx.lkdC   ;
            fprintf('current kdC     is  % 3.2e \n', exp(x(ikdC)));
            xhat.kdC = exp(x(ikdC)) ;
        end

        if (par.opt_R_Si == on)
            iR_Si = par.pindx.R_Si ;
            fprintf('current R_Si    is  % 3.2e \n', x(iR_Si)) ;
            xhat.R_Si = x(iR_Si)   ;
        end

        if (par.opt_rR == on)
            irR = par.pindx.lrR;
            fprintf('current rR      is  % 3.2e \n', exp(x(irR)));
            xhat.rR = exp(x(irR));
        end

        if (par.opt_cc == on)
            icc = par.pindx.lcc;
            fprintf('current cc      is  % 3.2e \n', exp(x(icc)));
            xhat.cc = exp(x(icc));
        end

        if (par.opt_dd == on)
            idd = par.pindx.ldd;
            fprintf('current dd      is  % 3.2e \n', exp(x(idd)));
            xhat.dd = exp(x(idd));
        end
    end
    % ------------------------------------------------------------
    if (par.Omodel == on)
        if (par.opt_O2C_T == on)
            iO2C_T = par.pindx.O2C_T;
            fprintf('current O2C_T   is  % 3.2e \n', x(iO2C_T));
            xhat.O2C_T = x(iO2C_T);
        end

        if (par.opt_rO2C == on)
            irO2C = par.pindx.lrO2C;
            fprintf('current rO2C    is  % 3.2e \n', exp(x(irO2C)));
            xhat.rO2C = exp(x(irO2C));
        end

        if (par.opt_O2P_T == on)
            iO2P_T = par.pindx.O2P_T;
            fprintf('current O2P_T   is  % 3.2e \n', x(iO2P_T));
            xhat.O2P_T = x(iO2P_T);
        end

        if (par.opt_rO2P == on)
            irO2P = par.pindx.lrO2P;
            fprintf('current rO2P    is  % 3.2e \n', exp(x(irO2P)));
            xhat.rO2P = exp(x(irO2P));
        end
    end
    % ------------------------------------------------------------
    % dsi
    if (par.Simodel==on)
        if (par.opt_dsi == on)
            idsi = par.pindx.ldsi;
            fprintf('current dsi     is  % 3.2e \n', exp(x(idsi)));
            xhat.dsi = exp(x(idsi));
        end

        % at
        if (par.opt_at == on)
            iat = par.pindx.lat;
            fprintf('current at      is  % 3.2e \n', exp(x(iat)));
            xhat.at = exp(x(iat));
        end

        % bt
        if (par.opt_bt == on)
            ibt = par.pindx.lbt;
            fprintf('current bt      is  % 3.2e \n', exp(x(ibt)));
            xhat.bt = exp(x(ibt));
        end

        % aa
        if (par.opt_aa == on)
            iaa = par.pindx.aa;
            fprintf('current aa      is  % 3.2e \n', x(iaa));
            xhat.aa = x(iaa);
        end

        % bb
        if (par.opt_bb == on)
            ibb = par.pindx.lbb;
            fprintf('current bb      is  % 3.2e \n\n', exp(x(ibb)));
            xhat.bb = exp(x(ibb));
        end
    end
	% ------------------------------------------------------------
    % Cell model
    if (par.Cellmodel==on)
        if (par.opt_Q10Photo == on)
            iQ10Photo = par.pindx.lQ10Photo;
            fprintf('current Q10Photo       is  % 3.2e \n', exp(x(iQ10Photo)));
            xhat.Q10Photo = exp(x(iQ10Photo));
        end
		if (par.opt_fStorage == on)
            ifStorage = par.pindx.lfStorage;
            fprintf('current fStorage       is  % 3.2e \n', exp(x(ifStorage)));
            xhat.fStorage = exp(x(ifStorage));
        end
		if (par.opt_PLip_PCutoff == on)
            iPLip_PCutoff = par.pindx.lPLip_PCutoff;
            fprintf('current PLip_PCutoff   is  % 3.2e \n', exp(x(iPLip_PCutoff)));
            xhat.PLip_PCutoff = exp(x(iPLip_PCutoff));
        end
		if (par.opt_PLip_scale == on)
            iPLip_Pscale = par.pindx.lPLip_scale;
            fprintf('current PLip_scale     is  % 3.2e \n', exp(x(iPLip_scale)));
            xhat.PLip_scale = exp(x(iPLip_scale));
        end
		if (par.opt_PStor_rCutoff == on)
            iPStor_rCutoff = par.pindx.lPStor_rCutoff;
            fprintf('current PStor_rCutoff  is  % 3.2e \n', exp(x(iPStor_rCutoff)));
            xhat.PStor_rCutoff = exp(x(iPStor_rCutoff));
        end
		if (par.opt_PStor_scale == on)
            iPStor_scale = par.pindx.lPStor_scale;
            fprintf('current PStor_scale    is  % 3.2e \n', exp(x(iPStor_scale)));
            xhat.PStor_scale = exp(x(iPStor_scale));
        end

	end
	% ------------------------------------------------------------
    x0 = x ;
    if (par.optim == on)
        save(par.fxhat, 'x0','xhat')
    end
end
