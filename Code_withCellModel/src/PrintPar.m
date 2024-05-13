function PrintPar(x, par);
    global iter
    on = true; off = false;
    %++++++++++ print out all parameters to the log file
    if iter == 0
        fprintf('All parameters \n')
        fprintf('-----------------------------------------------\n')
        fprintf('kappa_p is  % 3.2e \n', par.kappa_p)   ;
        fprintf('kappa_g is  % 3.2e \n', par.kappa_g)   ;
        fprintf('kappa_l is  % 3.2e \n', par.kappa_l)   ;
        fprintf('sigP    is  % 3.2e \n', par.sigP)      ;
        fprintf('gamma   is  % 3.2e \n', par.gamma)     ;
        fprintf('Q10P    is  % 3.2e \n', par.Q10P)      ;
        fprintf('kdP     is  % 3.2e \n', par.kdP)       ;
        fprintf('bP_T    is  % 3.2e \n', par.bP_T)      ;
        fprintf('bP      is  % 3.2e \n', par.bP)        ;
        fprintf('alpha   is  % 3.2e \n', par.alpha)     ;
        fprintf('beta    is  % 3.2e \n', par.beta)      ;
        if (par.Cmodel == on)
            fprintf('sigC    is  % 3.2e \n', par.sigC)  ;
            fprintf('kru     is  % 3.2e \n', par.kru)   ;
            fprintf('krd     is  % 3.2e \n', par.krd)   ;
            fprintf('etau    is  % 3.2e \n', par.etau)  ;
            fprintf('etad    is  % 3.2e \n', par.etad)  ;
            fprintf('kPIC    is  % 3.2e \n', par.kPIC)  ;
            fprintf('bC_T    is  % 3.2e \n', par.bC_T)  ;
            fprintf('bC      is  % 3.2e \n', par.bC)    ;
            fprintf('d       is  % 3.2e \n', par.d)     ;
            fprintf('Q10C    is  % 3.2e \n', par.Q10C)  ;
            fprintf('kdC     is  % 3.2e \n', par.kdC)   ;
            fprintf('R_Si    is  % 3.2e \n', par.R_Si)  ;
            fprintf('rR      is  % 3.2e \n', par.rR)    ;
            fprintf('cc      is  % 3.2e \n', par.cc)    ;
            fprintf('dd      is  % 3.2e \n', par.dd)    ;
        end 
        if (par.Cellmodel==on)
            fprintf('Q10Photo       is  % 3.2e \n', par.BIO.Q10Photo)   ;
			fprintf('fStorage       is  % 3.2e \n', par.BIO.fStorage)   ;
			fprintf('fRibE          is  % 3.2e \n', par.BIO.fRibE)      ;
			fprintf('PLip_PCutoff   is  % 3.2e \n', par.BIO.PLip_PCutoff) ;
			fprintf('PLip_scale     is  % 3.2e \n', par.BIO.PLip_scale)   ;
            fprintf('PStor_r0Cutoff is  % 3.2e \n', par.BIO.PStor_rCutoff);
			fprintf('PStor_scale    is  % 3.2e \n', par.BIO.PStor_scale)  ;
			fprintf('alphaS         is  % 3.2e \n', par.BIO.alphaS)     ;
			fprintf('alphaPLip      is  % 3.2e \n', par.BIO.alphaPLip)  ;
            fprintf('Q10Diffusivity is  % 3.2e \n', par.BIO.Q10Diffusivity) ;
			fprintf('Q10Bio         is  % 3.2e \n', par.BIO.Q10Bio)  ;
			fprintf('DNT0           is  % 3.2e \n', par.BIO.DNT0)  ;
			fprintf('DPT0           is  % 3.2e \n', par.BIO.DPT0)  ;
			fprintf('AMin           is  % 3.2e \n', par.BIO.AMin)  ;
			fprintf('PhiS           is  % 3.2e \n', par.BIO.PhiS)  ;
			fprintf('pDry           is  % 3.2e \n', par.BIO.pDry)  ;
			fprintf('rho            is  % 3.2e \n', par.BIO.rho)   ;
			fprintf('gammaDNA       is  % 3.2e \n', par.BIO.gammaDNA)   ;
            fprintf('gammaLipid     is  % 3.2e \n', par.BIO.gammaLipid) ;
			fprintf('fProtM         is  % 3.2e \n', par.BIO.fProtM)  ;
			fprintf('fProtL         is  % 3.2e \n', par.BIO.fProtL)  ;
			fprintf('PDNA           is  % 3.2e \n', par.BIO.PDNA)    ;
			fprintf('PRib           is  % 3.2e \n', par.BIO.PRib)    ;
			fprintf('PPhospholipid  is  % 3.2e \n', par.BIO.PPhospholipid) ;
			fprintf('NProt          is  % 3.2e \n', par.BIO.NProt)  ;
			fprintf('NDNA           is  % 3.2e \n', par.BIO.NDNA)  ;
			fprintf('NRib           is  % 3.2e \n', par.BIO.NRib)  ;
			fprintf('CProt          is  % 3.2e \n', par.BIO.CProt)  ;
			fprintf('CDNA           is  % 3.2e \n', par.BIO.CDNA)  ;
			fprintf('CPhospholipid  is  % 3.2e \n', par.BIO.CPhospholipid )  ;
			fprintf('CLipid         is  % 3.2e \n', par.BIO.CLipid)  ;
			fprintf('CRib           is  % 3.2e \n\n', par.BIO.CRib)  ;
		elseif (par.C2P_Tzmodel ==on)
			fprintf('ccT     is  % 3.2e \n', par.ccT)   ;
			fprintf('ddT     is  % 3.2e \n', par.ddT)   ;
		else
			fprintf('cc      is  % 3.2e \n', par.cc)    ;
            fprintf('dd      is  % 3.2e \n', par.dd)    ;
        end
        if (par.Omodel == on)
            fprintf('O2C_T   is  % 3.2e \n', par.O2C_T) ;
            fprintf('rO2C    is  % 3.2e \n', par.rO2C)  ;
        end
        if (par.Simodel==on)
            fprintf('dsi     is  % 3.2e \n', par.dsi)   ;
            fprintf('at      is  % 3.2e \n', par.at)    ;
            fprintf('bt      is  % 3.2e \n', par.bt)    ;
            fprintf('aa      is  % 3.2e \n', par.iaa)   ;
            fprintf('bb      is  % 3.2e \n\n', par.bb)  ;
        end
        fprintf('-----------------------------------------------\n\n')
    end 
    %++++++++++ print out parameters to the log file
    fprintf('Tunable parameters \n')
    fprintf('-----------------------------------------------\n')
    if (par.opt_sigP == on)
        isigP = par.pindx.lsigP   ;
        fprintf('current sigP    is  % 3.2e \n', exp(x(isigP)));
        xhat.sigP = exp(x(isigP)) ;
    end

    if (par.opt_Q10P == on)
        iQ10P = par.pindx.lQ10P   ;
        fprintf('current Q10P    is  % 3.2e \n', exp(x(iQ10P)));
        xhat.Q10P = exp(x(iQ10P)) ;
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
        if (par.opt_sigC == on)
            isigC = par.pindx.lsigC   ;
            fprintf('current sigC    is  % 3.2e \n', exp(x(isigC)));
            xhat.sigC = exp(x(isigC)) ;
        end
        
        if (par.opt_kru == on)
            ikru = par.pindx.lkru;
            fprintf('current kru     is  % 3.2e \n', exp(x(ikru)));
            xhat.kru = exp(x(ikru));
        end

        if (par.opt_krd == on)
            ikrd = par.pindx.lkrd;
            fprintf('current krd     is  % 3.2e \n', exp(x(ikrd)));
            xhat.krd = exp(x(ikrd));
        end

        if (par.opt_etau == on)
            ietau = par.pindx.letau;
            fprintf('current etau    is  % 3.2e \n', exp(x(ietau)));
            xhat.etau = exp(x(ietau));
        end

        if (par.opt_etad == on)
            ietad = par.pindx.letad;
            fprintf('current etad    is  % 3.2e \n', exp(x(ietad)));
            xhat.etad = exp(x(ietad));
        end

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

        if (par.opt_Q10C == on)
            iQ10C = par.pindx.lQ10C;
            fprintf('current Q10C    is  % 3.2e \n', exp(x(iQ10C)));
            xhat.Q10C = exp(x(iQ10C));
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

        % Temperature dependent C:P model
        if (par.opt_ccT == on)
            icc = par.pindx.ccT;
            fprintf('current ccT     is  % 3.2e \n', x(icc));
            xhat.ccT = x(icc);
        end
        if (par.opt_ddT == on)
            idd = par.pindx.lddT;
            fprintf('current ddT     is  % 3.2e \n', exp(x(idd)));
            xhat.ddT = exp(x(idd));
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
		if (par.opt_fRibE == on)
            ifRibE = par.pindx.tfRibE;
            fprintf('current fRibE          is  % 3.2e \n', 0.5*(1+tanh(x(ifRibE))) ); %fRibE = 0.5*(1+tanh(x(ifRibE)))
            xhat.fRibE = 0.5*(1+tanh(x(ifRibE)));
        end
		if (par.opt_kST0 == on)
            ikST0 = par.pindx.lkST0;
            fprintf('current kST0           is  % 3.2e \n', exp(x(ikST0)));
            xhat.kST0 = exp(x(ikST0));
        end
		if (par.opt_PLip_PCutoff == on)
            iPLip_PCutoff = par.pindx.lPLip_PCutoff;
            fprintf('current PLip_PCutoff   is  % 3.2e \n', exp(x(iPLip_PCutoff)));
            xhat.PLip_PCutoff = exp(x(iPLip_PCutoff));
        end
		if (par.opt_PLip_scale == on)
            iPLip_scale = par.pindx.lPLip_scale;
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
		if (par.opt_alphaS == on)
            ialphaS = par.pindx.lalphaS;
            fprintf('current alphaS         is  % 3.2e \n', exp(x(ialphaS)));
            xhat.alphaS = exp(x(ialphaS));
        end
		if (par.opt_gammaDNA == on)
            igammaDNA = par.pindx.tgammaDNA;
            fprintf('current gammaDNA       is  % 3.2e \n', 0.5*(1+tanh(x(igammaDNA))));
            xhat.gammaDNA = 0.5*(1+tanh(x(igammaDNA)));
        end
	end %end cell model params
	% ------------------------------------------------------------

    % Ensure only real part of parameter gets saved to avoid issues with
	% complex step derivatives
	nx = length(x);
	fnames = fieldnames(xhat);
	for ii = 1:nx
		xhat.(fnames{ii}) = real(xhat.(fnames{ii}));
	end

    % Save in neglogpost instead?
    x0 = x ;
    if (par.optim == on)
        save(par.fxhat, 'x0','xhat')
    end 
end
