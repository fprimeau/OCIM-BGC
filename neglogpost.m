function [f, fx, fxx, data] = neglogpost(x, par)
    global iter
    on = true; off = false;
	current_time = string(datetime('now')) ; %for runtime diagnostic purposes
	fprintf('current time: %s \n',current_time) ;

    % print and save current parameter values to
    % a file that is used to reset parameters ;
    if iter == 0
        PrintPara(x, par) ;
    end
    % reset parameters if optimization routine
    % suggests strange parameter values ;
    if iter < 3
        x = ResetPara(x, par) ;
    end
    % print current parameters
    if iter > 0
        PrintPara(x, par) ;
    end
    fprintf('current iteration is %d \n',iter) ;
    iter = iter + 1  ;

    nx   = length(x) ; % number of parameters
    dVt  = par.dVt   ;
    M3d  = par.M3d   ;
    iwet = par.iwet  ;
    nwet = par.nwet  ;
    %
    f    = 0 ;
    %%%%%%%%%%%%%%%%%%   Solve P    %%%%%%%%%%%%%%%%%%%%%%%%
	tic
    idip = find(par.po4raw(iwet) > 0.02) ;
    Wp   = d0(dVt(iwet(idip))/sum(dVt(iwet(idip)))) ;
    mu   = sum(Wp*par.po4raw(iwet(idip)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.po4raw(iwet(idip))-mu).^2)/sum(diag(Wp)) ;
    Wip  = Wp/var ;

    idop = find(par.dopraw(iwet) > 0.0) ;
    Wp   = d0(dVt(iwet(idop))/sum(dVt(iwet(idop)))) ;
    mu   = sum(Wp*par.dopraw(iwet(idop)))/sum(diag(Wp)) ;
    var  = sum(Wp*(par.dopraw(iwet(idop))-mu).^2)/sum(diag(Wp)) ;
    Wop  = par.pscale*Wp/var ;
    %
    [par, P, Px, Pxx] = eqPcycle(x, par) ;
    DIP = M3d+nan  ;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
    POP = M3d+nan  ;  POP(iwet) = P(1+1*nwet:2*nwet) ;
    DOP = M3d+nan  ;  DOP(iwet) = P(1+2*nwet:3*nwet) ;

    par.Px   = Px  ;
    par.Pxx  = Pxx ;
    par.DIP  = DIP(iwet) ;
    data.DIP = DIP ; data.POP = POP ; data.DOP = DOP ;
    % DIP error
    eip = DIP(iwet(idip)) - par.po4raw(iwet(idip)) ;
    eop = DOP(iwet(idop)) - par.dopraw(iwet(idop)) ;
    f  = f + 0.5*(eip.'*Wip*eip) + 0.5*(eop.'*Wop*eop);
	toc
    %%%%%%%%%%%%%%%%%%   End Solve P    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%   Solve Si       %%%%%%%%%%%%%%%%%%%%
    if (par.Simodel == on)
        isil = find(par.sio4raw(iwet)>0) ;
        Ws   = d0(dVt(iwet(isil))/sum(dVt(iwet(isil)))) ;
        mu   = sum(Ws*par.sDSi(iwet(isil)))/sum(diag(Ws)) ;
        var  = sum(Ws*(par.sDSi(iwet(isil))-mu).^2)/sum(diag(Ws));
        Ws   = Ws/var ;
        %
        [par,Si,Six,Sixx] = eqSicycle(x, par)   ;
        DSi = M3d+nan ;  DSi(iwet) = Si(1:nwet) ;
        bSi = M3d+nan ;  bSi(iwet) = Si(nwet+1:end) ;
        data.DSi = SI ;  data.bSi  = bSi ;
        % SiO error
        es = DSi(iwet(isil)) - par.sio4raw(iwet(isil)) ;
        f  = f + 0.5*(es.'*Ws*es) ;
    end
    %%%%%%%%%%%%%%%%%%   End Solve Si    %%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%     Solve for C2P with Cell model   %%%%%%%%%%%%%%%%%%%%%
	if (par.Cellmodel == on)
		tic
		[par, C2P, C2Px, C2Pxx] = eqC2Puptake(x, par, data); % this line replaces the rest of this section
		par.C2Px = C2Px;
		par.C2Pxx = C2Pxx;
%{
		iprod = find(M3d(:,:,1:2)); %production in top two layers
		P0 = data.DIP(iprod)./10^6; 			% convert ug/m^3 to g/m^3  data.DIP:[mmol/m^3] convert to mol/L
		N0 = par.no3obs(iprod)./10^6;   % convert ug/m^3 to g/m^3 <- NO [mmol/m^3 --> mol/L]
		T0 = par.Temp(iprod);
		Irr0 = par.PARobs(iprod);

   % keyboard;
		%parBIO=par.BIO;
        %M3dsurf=par.M3d(:,:,1:2);
		%dVtsurf=par.dVt(:,:,1:2);
		%pindx=par.pindx;
		%P0alt = par.po4raw(iprod)./10^6;
		%N0alt = par.no3raw(iprod)./10^6;
		%x0=x;
		%save('inputs_surf_CellCNP.mat','P0','N0','T0','Irr0','M3dsurf','dVtsurf','x0','parBIO','iprod','pindx','P0alt',N0alt');

 		[CellOut, parBIO] = CellCNP(par,x, P0,N0,T0,Irr0);
		par.BIO = parBIO;
		clear parBIO;
		par.CellOut.C2P = M3d*0;
		par.CellOut.N2P = M3d*0;
		par.CellOut.C2N = M3d*0;
		par.CellOut.LimType = NaN(size(M3d));
		par.CellOut.r = M3d*0;

		par.CellOut.C2P(iprod) = CellOut.CP;
		par.CellOut.N2P(iprod) = CellOut.NP;
		par.CellOut.C2N(iprod) = CellOut.CN;
		par.CellOut.LimType(iprod) = CellOut.LimType;
		par.CellOut.r(iprod) = CellOut.r;

		par.CellOut.C2P(isnan(par.CellOut.C2P)) = 0; %remove NaNs

		% Cell model derivatives
		% Q10Photo Derivatives
		if (par.opt_Q10Photo)
			par.CellOut.dC2P_dQ10Photo = M3d*0;
			par.CellOut.dC2P_dQ10Photo(iprod) = CellOut.dC2P_dQ10Photo;

			%second Derivatives w.r.t. Q10Photo
			xim = zeros(size(x));
			xim(par.pindx.lQ10Photo) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_dQ10Photo2 = M3d*0;
			par.CellOut.d2C2P_dQ10Photo2(iprod) = imag(CellOut.dC2P_dQ10Photo)./eps^3;

			if (par.opt_fStorage)
				par.CellOut.d2C2P_dfStorage_dQ10Photo = M3d*0;
				par.CellOut.d2C2P_dfStorage_dQ10Photo(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3;
			end
			if (par.opt_fRibE)
                par.CellOut.d2C2P_dfRibE_dQ10Photo = M3d*0;
				par.CellOut.d2C2P_dfRibE_dQ10Photo(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
            end
			if (par.opt_PLip_PCutoff)
				par.CellOut.d2C2P_dPCutoff_dQ10Photo = M3d*0;
				par.CellOut.d2C2P_dPCutoff_dQ10Photo(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
			end
			if (par.opt_PStor_rCutoff)
				par.CellOut.d2C2P_drCutoff_dQ10Photo = M3d*0;
				par.CellOut.d2C2P_drCutoff_dQ10Photo(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
			end
			if (par.opt_PStor_scale)
				par.CellOut.d2C2P_dPStorscale_dQ10Photo = M3d*0;
				par.CellOut.d2C2P_dPStorscale_dQ10Photo(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			end
			if (par.opt_PLip_scale)
				par.CellOut.d2C2P_dPLipscale_dQ10Photo = M3d*0;
				par.CellOut.d2C2P_dPLipscale_dQ10Photo(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
			end
			if (par.opt_alphaS)
				par.CellOut.d2C2P_dalphaS_dQ10Photo = M3d*0;
				par.CellOut.d2C2P_dalphaS_dQ10Photo(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			end
		end

		% fStorage Derivatives
		if (par.opt_fStorage)
			par.CellOut.dC2P_dfStorage = M3d*0;
			par.CellOut.dC2P_dfStorage(iprod) = real(CellOut.dC2P_dfStorage);

			%second Derivatives w.r.t. fStorage
			xim = zeros(size(x));
			xim(par.pindx.lfStorage) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_dfStorage2 = M3d*0;
			par.CellOut.d2C2P_dfStorage2(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3;

			%if (par.opt_Q10Photo)
			%	par.CellOut.d2C2P_dQ10Photo_dfStorage = M3d*0;
			%	par.CellOut.d2C2P_dQ10Photo_dfStorage(iprod) = imag(CellOut.dC2P_dQ10Photo)./eps^3;
			%end
			if (par.opt_fRibE)
                par.CellOut.d2C2P_dfRibE_dfStorage = M3d*0;
				par.CellOut.d2C2P_dfRibE_dfStorage(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
            end
			if (par.opt_PLip_PCutoff)
				par.CellOut.d2C2P_dPCutoff_dfStorage = M3d*0;
				par.CellOut.d2C2P_dPCutoff_dfStorage(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
			end
			if (par.opt_PLip_scale)
				par.CellOut.d2C2P_dPLipscale_dfStorage = M3d*0;
				par.CellOut.d2C2P_dPLipscale_dfStorage(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
			end
			if (par.opt_PStor_rCutoff)
				par.CellOut.d2C2P_drCutoff_dfStorage = M3d*0;
				par.CellOut.d2C2P_drCutoff_dfStorage(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
			end
			if (par.opt_PStor_scale)
				par.CellOut.d2C2P_dPStorscale_dfStorage = M3d*0;
				par.CellOut.d2C2P_dPStorscale_dfStorage(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			end
			if (par.opt_alphaS)
				par.CellOut.d2C2P_dalphaS_dfStorage = M3d*0;
				par.CellOut.d2C2P_dalphaS_dfStorage(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			end
		end

		%fRibE derivatives
		if (par.opt_fRibE)
			par.CellOut.dC2P_dfRibE  = M3d*0;
			par.CellOut.dC2P_dfRibE(iprod)  = real(CellOut.dC2P_dfRibE);

			%second Derivatives w.r.t. fRibE
			xim = zeros(size(x));
			xim(par.pindx.tfRibE) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_dfRibE2 = M3d*0;
			par.CellOut.d2C2P_dfRibE2(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;

            if (par.opt_PLip_PCutoff)
				par.CellOut.d2C2P_dPCutoff_dfRibE = M3d*0;
				par.CellOut.d2C2P_dPCutoff_dfRibE(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
			end
			if (par.opt_PLip_scale)
				par.CellOut.d2C2P_dPLipscale_dfRibE = M3d*0;
				par.CellOut.d2C2P_dPLipscale_dfRibE(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
			end
			if (par.opt_PStor_rCutoff)
				par.CellOut.d2C2P_drCutoff_dfRibE = M3d*0;
				par.CellOut.d2C2P_drCutoff_dfRibE(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
			end
			if (par.opt_PStor_scale)
				par.CellOut.d2C2P_dPStorscale_dfRibE = M3d*0;
				par.CellOut.d2C2P_dPStorscale_dfRibE(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			end
			if (par.opt_alphaS)
				par.CellOut.d2C2P_dalphaS_dfRibE = M3d*0;
				par.CellOut.d2C2P_dalphaS_dfRibE(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			end
        end

		% PLip_PCutoff derivatives
		if (par.opt_PLip_PCutoff)
			par.CellOut.dC2P_dPCutoff  = M3d*0;
			par.CellOut.dC2P_dPCutoff(iprod)  = real(CellOut.dC2P_dPCutoff);

			%second Derivatives w.r.t. PLip_PCutoff
			xim = zeros(size(x));
			xim(par.pindx.lPLip_PCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_dPCutoff2 = M3d*0;
			par.CellOut.d2C2P_dPCutoff2(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;

			if (par.opt_PLip_scale)
				par.CellOut.d2C2P_dPLipscale_dPCutoff = M3d*0;
				par.CellOut.d2C2P_dPLipscale_dPCutoff(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
			end
			if (par.opt_PStor_rCutoff)
				par.CellOut.d2C2P_drCutoff_dPCutoff = M3d*0;
				par.CellOut.d2C2P_drCutoff_dPCutoff(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
			end
			if (par.opt_PStor_scale)
				par.CellOut.d2C2P_dPStorscale_dPCutoff = M3d*0;
				par.CellOut.d2C2P_dPStorscale_dPCutoff(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			end
			if (par.opt_alphaS)
				par.CellOut.d2C2P_dalphaS_dPCutoff = M3d*0;
				par.CellOut.d2C2P_dalphaS_dPCutoff(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			end
		end

		% PLip_scale Derivatives
		if (par.opt_PLip_scale)
			par.CellOut.dC2P_dPLip_scale = M3d*0;
			par.CellOut.dC2P_dPLip_scale(iprod) = real(CellOut.dC2P_dPLipscale);

			%second Derivatives w.r.t. PLip_scale
			xim = zeros(size(x));
			xim(par.pindx.lPLip_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_dPLipscale2 = M3d*0;
			par.CellOut.d2C2P_dPLipscale2(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;

			if (par.opt_PStor_rCutoff)
				par.CellOut.d2C2P_drCutoff_dPLipscale = M3d*0;
				par.CellOut.d2C2P_drCutoff_dPLipscale(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
			end
			if (par.opt_PStor_scale)
				par.CellOut.d2C2P_dPStorscale_dPLipscale = M3d*0;
				par.CellOut.d2C2P_dPStorscale_dPLipscale(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			end
			if (par.opt_alphaS)
				par.CellOut.d2C2P_dalphaS_dPLipscale = M3d*0;
				par.CellOut.d2C2P_dalphaS_dPLipscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			end
		end

		% PStor_rCutoff derivatives
		if (par.opt_PStor_rCutoff)
			par.CellOut.dC2P_drCutoff  = M3d*0;
			par.CellOut.dC2P_drCutoff(iprod)  = real(CellOut.dC2P_drCutoff);

			%second Derivatives w.r.t. PStor_rCutoff
			xim = zeros(size(x));
			xim(par.pindx.lPStor_rCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_drCutoff2 = M3d*0;
			par.CellOut.d2C2P_drCutoff2(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;

			if (par.opt_PStor_scale)
				par.CellOut.d2C2P_dPStorscale_drCutoff = M3d*0;
				par.CellOut.d2C2P_dPStorscale_drCutoff(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			end
			if (par.opt_alphaS)
				par.CellOut.d2C2P_dalphaS_drCutoff = M3d*0;
				par.CellOut.d2C2P_dalphaS_drCutoff(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			end
		end

		% PStor_scale derivatives
		if (par.opt_PStor_scale)
			par.CellOut.dC2P_dPStor_scale = M3d*0;
			par.CellOut.dC2P_dPStor_scale(iprod) = real(CellOut.dC2P_dPStorscale);

			%second Derivatives w.r.t. PStor_scale
			xim = zeros(size(x));
			xim(par.pindx.lPStor_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_dPStorscale2 = M3d*0;
			par.CellOut.d2C2P_dPStorscale2(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;

			if (par.opt_alphaS)
				par.CellOut.d2C2P_dalphaS_dPStorscale = M3d*0;
				par.CellOut.d2C2P_dalphaS_dPStorscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			end
		end

		% alphaS derivatives
		if (par.opt_alphaS)
			par.CellOut.dC2P_dalphaS = M3d*0;
			par.CellOut.dC2P_dalphaS(iprod) = real(CellOut.dC2P_dalphaS);

			%second Derivatives w.r.t. alphaS
			xim = zeros(size(x));
			xim(par.pindx.lalphaS) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_dalphaS2 = M3d*0;
			par.CellOut.d2C2P_dalphaS2(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
		end

		% par.CellOut.dC2P_dQ10Photo = M3d*0;
		% par.CellOut.dC2P_dfStorage = M3d*0;
		% par.CellOut.dC2P_dPCutoff  = M3d*0;
		% par.CellOut.dC2P_drCutoff  = M3d*0;
		% par.CellOut.dC2P_dPStor_scale = M3d*0;
		% par.CellOut.dC2P_dPLip_scale = M3d*0;
		% %par.CellOut.dC2P_dalphaS = M3d*0;
		%
        % par.CellOut.dC2P_dQ10Photo(iprod) = CellOut.dC2P_dQ10Photo;
		% par.CellOut.dC2P_dfStorage(iprod) = CellOut.dC2P_dfStorage;
		% par.CellOut.dC2P_dPCutoff(iprod)  = Cellout.dC2P_dPCutoff;
		% par.CellOut.dC2P_drCutoff(iprod)  = CellOut.dC2P_drCutoff;
		% par.CellOut.dC2P_dPStor_scale(iprod) = CellOut.dC2P_dPStorscale;
		% par.CellOut.dC2P_dPLip_scale(iprod) = CellOut.dC2P_dPLipscale;
		% %par.CellOut.dC2P_dalphaS(iprod) = CellOut.dC2P_dalphaS;
%}
		%data.CellOut = par.CellOut;
		data.CellOut.C2P =  par.CellOut.C2P;
		data.CellOut.N2P = par.CellOut.N2P;
		data.CellOut.C2N = par.CellOut.C2N;
		data.CellOut.LimType = par.CellOut.LimType;
		data.CellOut.r = par.CellOut.r;
		toc
	end

    %%%%%%%%%%%%%%%%%%     Solve C   %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Cmodel == on)
		tic
        idic = find(par.dicraw(iwet) > 0) ;
        Wic  = d0(dVt(iwet(idic))/sum(dVt(iwet(idic)))) ;
        mu   = sum(Wic*par.dicraw(iwet(idic)))/sum(diag(Wic)) ;
        var  = sum(Wic*(par.dicraw(iwet(idic))-mu).^2)/sum(diag(Wic));
        Wic  = Wic/var  ;

        ialk = find(par.alkraw(iwet)>0) ;
        Wlk  = d0(dVt(iwet(ialk))/sum(dVt(iwet(ialk)))) ;
        mu   = sum(Wlk*par.alkraw(iwet(ialk)))/sum(diag(Wlk)) ;
        var  = sum(Wlk*(par.alkraw(iwet(ialk))-mu).^2)/sum(diag(Wlk));
        Wlk  = Wlk/var  ;

        idoc = find(par.docraw(iwet)>0) ;
        Woc  = d0(dVt(iwet(idoc))/sum(dVt(iwet(idoc)))) ;
        mu   = sum(Woc*par.docraw(iwet(idoc)))/sum(diag(Woc)) ;
        var  = sum(Woc*(par.docraw(iwet(idoc))-mu).^2)/sum(diag(Woc));
        Woc  = par.cscale*Woc/var ;

        [par, C, Cx, Cxx] = eqCcycle(x, par) ;
        DIC = M3d+nan ;  DIC(iwet) = C(0*nwet+1:1*nwet) ;
        POC = M3d+nan ;  POC(iwet) = C(1*nwet+1:2*nwet) ;
        DOC = M3d+nan ;  DOC(iwet) = C(2*nwet+1:3*nwet) ;
        PIC = M3d+nan ;  PIC(iwet) = C(3*nwet+1:4*nwet) ;
        ALK = M3d+nan ;  ALK(iwet) = C(4*nwet+1:5*nwet) ;

        par.DIC = DIC(iwet) ;
        par.DOC = DOC(iwet) ;
        DIC = DIC + par.dicant  ;
        par.Cx   = Cx  ;  par.Cxx  = Cxx ;
        data.DIC = DIC ;  data.POC = POC ;
        data.DOC = DOC ;  data.PIC = PIC ;
        data.ALK = ALK ;
        % DIC error
        eic = DIC(iwet(idic)) - par.dicraw(iwet(idic)) ;
        eoc = DOC(iwet(idoc)) - par.docraw(iwet(idoc)) ;
        elk = ALK(iwet(ialk)) - par.alkraw(iwet(ialk)) ;
        f   = f + 0.5*(eic.'*Wic*eic) + 0.5*(eoc.'*Woc*eoc) + ...
              0.5*(elk.'*Wlk*elk);
		toc
    end

% debugging gradient --> fixed!
%	[ibadx,ibady] = find(Cx ~= real(Cx));
%	ibadparam = unique(ibady)
%	keyboard;
    %%%%%%%%%%%%%%%%%%   End Solve C    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%   Solve O    %%%%%%%%%%%%%%%%%%%%%%%%
    if (par.Omodel == on)
		tic
        io2 = find(par.o2raw(iwet)>0) ;
        Wo  = d0(dVt(iwet(io2))/sum(dVt(iwet(io2)))) ;
        mu  = sum(Wo*par.o2raw(iwet(io2)))/sum(diag(Wo)) ;
        var = sum(Wo*(par.o2raw(iwet(io2))-mu).^2)/sum(diag(Wo)) ;
        Wo  = Wo/var ;
        %
        [par, O, Ox, Oxx] = eqOcycle(x, par) ;
        O2 = M3d+nan ;  O2(iwet) = O ;
        data.O2 = O2 ;
        eo = O2(iwet(io2)) - par.o2raw(iwet(io2)) ;
        f  = f + 0.5*(eo.'*Wo*eo)   ;
		toc
    end
    %%%%%%%%%%%%%%%%%%   End Solve O    %%%%%%%%%%%%%%%%%%%%
    fprintf('current objective function value is %3.3e \n\n',f)
    if mod(iter, 10) == 0
        save(par.fname, 'data')
    end
    %% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % calculate gradient
    if (nargout > 1)
        fx = zeros(length(x), 1)   ;
        ipx  = Px(0*nwet+1:nwet,:) ;
        opx  = Px(2*nwet+1:end ,:) ;
        npx = par.npx              ;
        % ---------------------------------
        for ji = 1 : npx
            fx(ji) = eip.'*Wip*ipx(idip,ji) + eop.'*Wop*opx(idop,ji);
        end
        % ---------------------------------
        if (par.Simodel == on & par.Omodel == off & par.Cmodel == off)
            sx  = Six(1:nwet,:) ;
            nsx = par.nsx       ;
            % --------------------------
            for ji = 1 : npx+nsx
                fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji) ;
            end
        end
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == off)
            icx = Cx(0*nwet+1 : 1*nwet, :) ;
            ocx = Cx(2*nwet+1 : 3*nwet, :) ;
            lkx = Cx(4*nwet+1 : 5*nwet, :) ;
            ncx = par.ncx   ;
			nbx = par.nbx   ; % treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
            % --------------------------
            for ji = 1 : npx+ncx+nbx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         fx(ji);
            end
        end
        % ----------------------------------
        % ----------------------------------
        if (par.Omodel == on & par.Simodel == off)
            ox  = Ox(1:nwet, :) ;
            nox = par.nox       ;
            % --------------------------
            for ji = 1 : npx+ncx+nbx+nox	%treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
                fx(ji) = fx(ji) + eo.'*Wo*ox(io2, ji) ;
            end
        end
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == on & par.Omodel == off)
            ncx = par.ncx ;
			nbx = par.nbx ; % treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
            nsx = par.nsx ;
            icx = Cx(0*nwet+1 : 1*nwet, :) ;
            ocx = Cx(2*nwet+1 : 3*nwet, :) ;
            lkx = Cx(4*nwet+1 : 5*nwet, :) ;
            sx  = Six(1:nwet, :) ;
            % --------------------------
            for ji = 1 : npx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         es.'*Ws*sx(isil, ji) + fx(ji) ;
            end
            % --------------------------
            for ji = npx+1 : npx+ncx+nbx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         fx(ji);
            end
            % --------------------------
            for ji = npx+ncx+nbx+1 : npx+ncx+nbx+nsx
                fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji) ;
            end
        end
        % ----------------------------------
        % ----------------------------------
        if (par.Cmodel == on & par.Simodel == on & par.Omodel == on)
            ncx = par.ncx ;
			nbx = par.nbx ; % treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
            nox = par.nox ;
            nsx = par.nsx ;
            icx = Cx(0*nwet+1 : 1*nwet, :) ;
            ocx = Cx(2*nwet+1 : 3*nwet, :) ;
            lkx = Cx(4*nwet+1 : 5*nwet, :) ;
            ox  = Ox(1:nwet,:) ;
            sx  = Six(1:nwet,:);
            % --------------------------
            for ji = 1 : npx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         eo.'*Wo*ox(io2 , ji) + ...
                         es.'*Ws*sx(isil, ji) + fx(ji) ;
            end
            % --------------------------
            for ji = npx+1 : npx+ncx+nbx
                fx(ji) = eic.'*Wic*icx(idic, ji) + ...
                         eoc.'*Woc*ocx(idoc, ji) + ...
                         elk.'*Wlk*lkx(ialk, ji) + ...
                         fx(ji);
            end
            % --------------------------
            for ji = npx+ncx+nbx+1 : npx+ncx+nbx+nox       % BUG: needs npx+ncx starting indx
                fx(ji) = fx(ji) + eo.'*Wo*ox(io2, ji) ;
            end
            % --------------------------
            for ji = npx+ncx+nbx+nox+1 : npx+ncx+nbx+nox+nsx
                fx(ji) = fx(ji) + es.'*Ws*sx(isil, ji) ;
            end
        end
        % ----------------------------------
	end
    %
    if (nargout>2)
        fxx = sparse(npx, npx)  ;
        ipxx = Pxx(0*nwet+1 : 1*nwet, :) ;
        opxx = Pxx(2*nwet+1 : 3*nwet, :) ;
        % ----------------------------------------------------------------
        kk = 0;
        for ju = 1:npx
            for jo = ju:npx
                kk = kk + 1 ;
                fxx(ju,jo) = ...
                    ipx(idip,ju).'*Wip*ipx(idip,jo) + eip.'*Wip*ipxx(idip,kk) + ...
                    opx(idop,ju).'*Wop*opx(idop,jo) + eop.'*Wop*opxx(idop,kk);
                % Simodel
                if par.Simodel == on
                    sxx = Sixx(1:nwet,:);
                    fxx(ju,jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);
                end
                % Cmodel
                if par.Cmodel == on
                    icxx = Cxx(0*nwet+1 : 1*nwet, :) ;
                    ocxx = Cxx(2*nwet+1 : 3*nwet, :) ;
                    lkxx = Cxx(4*nwet+1 : 5*nwet, :) ;
                    fxx(ju,jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;
                end
                % Omodel
                if par.Omodel == on
                    oxx = Oxx(1:nwet,:);
                    fxx(ju,jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);
                end
                % make Hessian symetric;
                fxx(jo, ju) = fxx(ju, jo);
            end
        end
        kpp = kk;
        % ----------------------------------------------------------------
        % C model
			% treating cell model as C-cycle parameters (replace ncx with ncx+nbx)
        if (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
            %
            tmp = sparse(npx+ncx+nbx, npx+ncx+nbx); %
            tmp(1:npx,1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P model parameters and C parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo,ju) = fxx(ju, jo);
                end
            end
			% P model parameters and Cell parameters
            for ju = 1:npx
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo,ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % Only C parameters
            for ju = (npx+1):(npx+ncx)  %to npx+ncx+nbx to include cell params (check order in Cxx)
                for jo = ju:(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
			% -------------------------------
            % C model parameters and Cell model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo,ju) = fxx(ju, jo);
                end
            end
			% -------------------------------
            % Only Cell parameters
            for ju = (npx+ncx+1):(npx+ncx+nbx)
                for jo = ju:(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
        end
        % ----------------------------------------------------------------
        % C and O model
        if (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
            %
            tmp = sparse(npx+ncx+nbx+nox, npx+ncx+nbx+nox);
            tmp(1:npx,1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P and C model parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
			% P and Cell model parameters
            for ju = 1:npx
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % C model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = ju:(npx+ncx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
			% C and cell model parameters
            for ju = (npx+1):(npx+ncx)
                for jo = (npx+ncx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
			% Cell model parameters
            for ju = (npx+ncx+1):(npx+ncx+nbx)
                for jo = ju:(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % P and O model parameters
            for ju = 1:npx
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % C and O model parameters     % need to add cell model with O
            for ju = (npx+1):(npx+ncx+nbx)
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % O model parameters
            for ju = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                for jo = ju:(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2, ju).'*Wo*ox(io2, jo) + eo.'*Wo*oxx(io2, kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
        end
        % ----------------------------------------------------------------
        % Si model on; Cmodel off; Omodel off;
        if (par.Cmodel == off & par.Omodel == off & par.Simodel == on)
            % --------------------------
            tmp = sparse(npx+nsx, npx+nsx);
            tmp(1:npx,1:npx) = fxx;
            fxx = tmp;
            % --------------------------
            % P model parameters and Si parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -----------------------------
            % Only Si parameters
            for ju = (npx+1):(npx+nsx)
                for jo = ju:(npx+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
        end
        % ----------------------------------------------------------------
        % Cmodel on; O model off; and Simodel on;
        if (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
            %
            tmp = sparse(npx+ncx+nbx+nsx, npx+ncx+nbx+nsx);
            tmp(1:npx, 1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P and C model parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) ;

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % C model parameters
            for ju = (npx+1):(npx+ncx+nbx)
                for jo = ju:(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % P model parameters and Si parameters
            kk = kpp; % starting fro P-P parameters
            for ju = 1:npx
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -----------------------------
            % Only Si parameters
            for ju = (npx+ncx+nbx+1):(npx+ncx+nbx+nsx)
                for jo = ju:(npx+ncx+nbx+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
        end
        % ----------------------------------------------------------------
        % Cmodel on; O model on; and Simodel on;
        if (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
            %
            tmp = sparse(npx+ncx+nbx+nox+nsx, npx+ncx+nbx+nox+nsx);
            tmp(1:npx, 1:npx) = fxx;
            fxx = tmp;
            % -------------------------------
            % P and C model parameters
            for ju = 1:npx
                for jo = (npx+1):(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % C model parameters
            for ju = (npx+1):(npx+ncx+nbx)
                for jo = ju:(npx+ncx+nbx)
                    kk = kk + 1;
                    fxx(ju, jo) = ...
                        fxx(ju, jo) + ...
                        icx(idic,ju).'*Wic*icx(idic,jo) + eic.'*Wic*icxx(idic,kk) + ...
                        ocx(idoc,ju).'*Woc*ocx(idoc,jo) + eoc.'*Woc*ocxx(idoc,kk) + ...
                        lkx(ialk,ju).'*Wlk*lkx(ialk,jo) + elk.'*Wlk*lkxx(ialk,kk) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % P and O model parameters
            for ju = 1:npx
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % C and O model parameters
            for ju = (npx+1):(npx+ncx+nbx)
                for jo = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2,ju).'*Wo*ox(io2,jo) + eo.'*Wo*oxx(io2,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -------------------------------
            % O model parameters
            for ju = (npx+ncx+nbx+1):(npx+ncx+nbx+nox)
                for jo = ju:(npx+ncx+nbx+nox)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        ox(io2, ju).'*Wo*ox(io2, jo) + eo.'*Wo*oxx(io2, kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % --------------------------
            % P model parameters and Si parameters
            kk = kpp; % starting fro P-P parameters
            for ju = 1:npx
                for jo = (npx+ncx+nbx+nox+1):(npx+ncx+nbx+nox+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
            % -----------------------------
            % Only Si parameters
            for ju = (npx+ncx+nbx+nox+1):(npx+ncx+nbx+nox+nsx)
                for jo = ju:(npx+ncx+nbx+nox+nsx)
                    kk = kk + 1;
                    fxx(ju, jo) = fxx(ju, jo) + ...
                        sx(isil,ju).'*Ws*sx(isil,jo) + es.'*Ws*sxx(isil,kk);

                    fxx(jo, ju) = fxx(ju, jo);
                end
            end
        end
    end

end % end function
