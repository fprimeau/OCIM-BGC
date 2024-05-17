%CellUptake
function [par, C2P, C2Px, C2Pxx, C2Ppxx] = eqC2Puptake(x, par, data)
    on = true; off = false;
    % temporary variable for testing methods
    %TestVer = 1 %structure with all derivs as fields
    TestVer = 2; %matrix derivative matrix

    M3d =par.M3d;
    pindx = par.pindx;
    iwet = par.iwet;
    nwet = par.nwet;

	% if par.Cellmodel == off
	% 	DIP = par.DIP ;
	% 	C2P = 1./(cc*DIP + dd);
	% end
		% data.DIP units = mmol/m^3

		iprod = find(M3d(:,:,1:2)); %production in top two layers
		N0 = par.no3obs(iprod)./10^6;   % *par.permil ; convert [umol/kg --> mmol/m^3 --> mol/L]
		T0 = par.Temp(iprod);
		Irr0 = par.PARobs(iprod);

		if par.dynamicP == on
			P0 = data.DIP(iprod)./10^6;		% data.DIP:[mmol/m^3] convert to mol/L
		else
			P0 = par.po4obs(iprod)./10^6;	% *par.permil ; convert [umol/kg --> mmol/m^3 --> mol/L]
		end

		if isfield(par,'future_case_options')
			tmp = par.future_case_options;
			if par.future_case_options.Pproj_in_C2P
				P0 = par.po4proj(iprod)./10^6;	% convert [mmol/m^3 --> mol/L]
				if par.dynamicP == on
					fprintf('Warning: curently, there is no option to run future case using modelled P instead of observed field in the cell model. Using projected PO4obs instead. \n')
				end
			end
			if par.future_case_options.Nproj_in_C2P
				N0 = par.no3proj(iprod)./10^6;	% convert [mmol/m^3 --> mol/L]
			end
			if par.future_case_options.Tproj_in_C2P
				T0 = par.Temp_proj(iprod);
			end
		end

		fprintf('Solving Cell model ...\n') ;
		%set negative phosphate values to smallest positive concentration.
	    % (negative values break CellCHNOP code)
	    fprintf('replacing %d negative Phosphate concentrations with the minimum positive concentration \n',length(P0(P0<0)))
		negDIPindx = (P0<0);
	    P0(P0<0)= real(min(P0(P0>=0)));

		fprintf('replacing %d negative DIN concentrations with the minimum positive concentration \n',length(N0(N0<0)))
		negDINindx = (N0<0);
	    N0(N0<0)= real(min(N0(N0>=0)));

		[CellOut, parBIO] = CellCHNOP(par,x, P0,N0,T0,Irr0);
		par.BIO = parBIO;
		clear parBIO;
		par.CellOut.C2P = M3d*0;
		par.CellOut.N2P = M3d*0;
		par.CellOut.C2N = M3d*0;
		par.CellOut.O2C = M3d*0;

		par.CellOut.C2P(iprod) = CellOut.CP;
		par.CellOut.N2P(iprod) = CellOut.NP;
		par.CellOut.C2N(iprod) = CellOut.CN;
		par.CellOut.O2C(iprod) = CellOut.RQtotalO2toC;

		par.CellOut.C2P(isnan(par.CellOut.C2P)) = 0; %remove NaNs

		% save cell model allocations for analysis
		par.CellOut.LimType = M3d(:,:,1:2)*NaN;
		par.CellOut.r       = M3d(:,:,1:2)*NaN;
		par.CellOut.mu      = M3d(:,:,1:2)*NaN;
		par.CellOut.E       = M3d(:,:,1:2)*NaN;
		par.CellOut.L       = M3d(:,:,1:2)*NaN;
		par.CellOut.A       = M3d(:,:,1:2)*NaN;
		par.CellOut.PLip    = M3d(:,:,1:2)*NaN;
		par.CellOut.PStor   = M3d(:,:,1:2)*NaN;
		par.CellOut.QP      = M3d(:,:,1:2)*NaN;
		par.CellOut.QC      = M3d(:,:,1:2)*NaN;

		par.CellOut.LimType(iprod) = CellOut.LimType;
		par.CellOut.r(iprod)       = CellOut.r;
		par.CellOut.mu(iprod)      = CellOut.mu;
		par.CellOut.E(iprod)       = CellOut.E;
		par.CellOut.L(iprod)       = CellOut.L;
		par.CellOut.A(iprod)       = CellOut.A;
		par.CellOut.PLip(iprod)    = CellOut.PLip;
		par.CellOut.PStor(iprod)   = CellOut.PStor;
		par.CellOut.QP(iprod)      = CellOut.QP;
		par.CellOut.QC(iprod)      = CellOut.QC;


        C2P = par.CellOut.C2P(iwet);
		%data.CellOut = par.CellOut;


%% ------------------------------------------------
%%%%% Gradient %%%%%

		% for points where DIP was reset from a negative value, set the derivative to 0
		CellOut.dC2P_dDIP(negDIPindx) = 0;

		%% always save this for sensitivity test
		% cell model uses DIP in units of mol/L
		% need to convert output that is function of DIP back to mmol/m^3
		% DIP [mol/L] = 1e-6 * DIP [mmol/m^3]
		dDIPmolperL_dDIPmmolperm3 = 1e-6;

		% if par.dynamicP == on
		% 	dC2P_dDIPmolperL = M3d*NaN;
		% 	dC2P_dDIPmolperL(iprod) = CellOut.dC2P_dDIP;
		% 	dC2P_dDIP = dC2P_dDIPmolperL * dDIPmolperL_dDIPmmolperm3;
		% else
		% 	dC2P_dDIPmolperL = M3d*NaN;
		% 	dC2P_dDIPmolperL(iprod) = 0;
		% 	dC2P_dDIP = dC2P_dDIPmolperL * dDIPmolperL_dDIPmmolperm3;
		% end

		dC2P_dDIPmolperL = M3d*NaN;
		dC2P_dDIPmolperL(iprod) = CellOut.dC2P_dDIP;
		dC2P_dDIP = dC2P_dDIPmolperL * dDIPmolperL_dDIPmmolperm3;

		par.CellOut.dC2P_dDIP = dC2P_dDIP(:,:,1:2);

%%%%% Gradient v2
if TestVer == 2
C2Px = [];
C2Pxx = [];
    if par.optim == off
        C2Px = [];
    elseif (par.optim & nargout > 1)
        % gradient of uptake operator
        nbx  = par.nbx; ncx=par.ncx; npx =par.npx;
        C2Px  = zeros(nwet,npx+ncx+nbx);

		if par.npx > 0
	        DIPx = par.Px(1:nwet,:);
			% define C2Px(:,pindx.sigP) [for all P model parameters]
			%Then in eqCcycle: +dDIC_dC2P* C2Px(:,pindx.sigP)  ;where dDIC_dC2P = -(I+(1-sigP)*RR)*G)

			%--------- P Model Parameters ----------------
			if par.dynamicP == on
				% C2P is a function of p model parameters through the DIP term
				if (par.opt_sigP == on)
					C2Px(:,par.pindx.lsigP) = d0(dC2P_dDIP(iwet))*DIPx(:,par.pindx.lsigP) ;
				end
				if (par.opt_Q10P == on)
			        C2Px(:,par.pindx.lQ10P) = d0(dC2P_dDIP(iwet))*DIPx(:,par.pindx.lQ10P) ;
				end
				if (par.opt_kdP   == on)
					C2Px(:,par.pindx.lkdP) = dC2P_dDIP(iwet).*DIPx(:,par.pindx.lkdP) ;
				end
				if (par.opt_bP_T  == on)
					C2Px(:,par.pindx.bP_T) = dC2P_dDIP(iwet).*DIPx(:,par.pindx.bP_T) ;
				end
				if (par.opt_bP    == on)
					C2Px(:,par.pindx.lbP) = dC2P_dDIP(iwet).*DIPx(:,par.pindx.lbP) ;
				end
				if (par.opt_alpha == on)
					C2Px(:,par.pindx.lalpha) = dC2P_dDIP(iwet).*DIPx(:,par.pindx.lalpha) ;
				end
				if (par.opt_beta  == on)
					C2Px(:,par.pindx.lbeta) = dC2P_dDIP(iwet).*DIPx(:,par.pindx.lbeta) ;
				end
			else
				% cell model uses observed PO4, so DIP is not a function of the P model parameters.
				dC2P_dPparam = M3d*0;
				if (par.opt_sigP == on)
					C2Px(:,par.pindx.lsigP) = dC2P_dPparam(iwet) ;
				end
				if (par.opt_Q10P == on)
			        C2Px(:,par.pindx.lQ10P) = dC2P_dPparam(iwet) ;
				end
				if (par.opt_kdP   == on)
					C2Px(:,par.pindx.lkdP) = dC2P_dPparam(iwet) ;
				end
				if (par.opt_bP_T  == on)
					C2Px(:,par.pindx.bP_T) = dC2P_dPparam(iwet) ;
				end
				if (par.opt_bP    == on)
					C2Px(:,par.pindx.lbP) = dC2P_dPparam(iwet) ;
				end
				if (par.opt_alpha == on)
					C2Px(:,par.pindx.lalpha) = dC2P_dPparam(iwet) ;
				end
				if (par.opt_beta  == on)
					C2Px(:,par.pindx.lbeta) = dC2P_dPparam(iwet) ;
				end
			end % end if par.dynamicP ==on ; else
		end
		%--------------------------------------------

		%--------- Cell Model Parameters ----------------
		dQ10_lQ10Photo = par.BIO.Q10Photo;
		dfStor_lfStorage = par.BIO.fStorage;
		tfRibE = atanh(2*par.BIO.fRibE-1);
		dfRibE_tfRibE = 0.5*sech(tfRibE)^2;
		dkST0_lkST0 = par.BIO.kST0;
		dPCutoff_lPCutoff = par.BIO.PLip_PCutoff;
		dPLipscale_lPLipscale = par.BIO.PLip_scale;
		drCutoff_lrCutoff = par.BIO.PStor_rCutoff;
		dPStorscale_lPStorscale = par.BIO.PStor_scale;
		dalphaS_lalphaS = par.BIO.alphaS;
		tgammaDNA = atanh(2*par.BIO.gammaDNA-1);
		dgammaDNA_tgammaDNA = 0.5*sech(tgammaDNA)^2;

        if (par.opt_Q10Photo)
            dC2P_dQ10Photo = M3d*0;
			dC2P_dQ10Photo(iprod) = CellOut.dC2P_dQ10Photo;
            C2Px(:,par.pindx.lQ10Photo) = dC2P_dQ10Photo(iwet) *dQ10_lQ10Photo;
			par.CellOut.dC2P_dQ10Photo = dC2P_dQ10Photo;
        end
        if (par.opt_fStorage)
            dC2P_dfStorage = M3d*0;
			dC2P_dfStorage(iprod) = CellOut.dC2P_dfStorage;
            C2Px(:,par.pindx.lfStorage) = dC2P_dfStorage(iwet) * dfStor_lfStorage;
			par.CellOut.dC2P_dfStorage = dC2P_dfStorage;
        end
        if par.opt_fRibE
            dC2P_dfRibE  = M3d*0;
			dC2P_dfRibE(iprod)  = CellOut.dC2P_dfRibE;
            C2Px(:,par.pindx.tfRibE) = dC2P_dfRibE(iwet) * dfRibE_tfRibE;
			par.CellOut.dC2P_dfRibE = dC2P_dfRibE;
        end
		if par.opt_kST0
            dC2P_dkST0  = M3d*0;
			dC2P_dkST0(iprod)  = CellOut.dC2P_dkST0;
            C2Px(:,par.pindx.lkST0) = dC2P_dkST0(iwet) * dkST0_lkST0;
			par.CellOut.dC2P_dkST0 = dC2P_dkST0;
        end
		if (par.opt_PLip_PCutoff)
            dC2P_dPCutoff  = M3d*0;
			dC2P_dPCutoff(iprod)  = CellOut.dC2P_dPCutoff;
			C2Px(:,par.pindx.lPLip_PCutoff) = dC2P_dPCutoff(iwet) * dPCutoff_lPCutoff;
			par.CellOut.dC2P_dPCutoff = dC2P_dPCutoff;
        end
        if (par.opt_PLip_scale)
            dC2P_dPLip_scale = M3d*0;
			dC2P_dPLip_scale(iprod) = CellOut.dC2P_dPLipscale;
			C2Px(:,par.pindx.lPLip_scale) = dC2P_dPLip_scale(iwet) * dPLipscale_lPLipscale;
			par.CellOut.dC2P_dPLip_scale = dC2P_dPLip_scale;
		end
		if (par.opt_PStor_rCutoff)
            dC2P_drCutoff  = M3d*0;
			dC2P_drCutoff(iprod)  = CellOut.dC2P_drCutoff;
			C2Px(:,par.pindx.lPStor_rCutoff) = dC2P_drCutoff(iwet) * drCutoff_lrCutoff;
			par.CellOut.dC2P_drCutoff = dC2P_drCutoff;
        end
        if par.opt_PStor_scale
            dC2P_dPStor_scale = M3d*0;
			dC2P_dPStor_scale(iprod) = CellOut.dC2P_dPStorscale;
			C2Px(:,par.pindx.lPStor_scale) = dC2P_dPStor_scale(iwet) * dPStorscale_lPStorscale;
			par.CellOut.dC2P_dPStor_scale = dC2P_dPStor_scale;
        end
        if par.opt_alphaS
            dC2P_dalphaS = M3d*0;
			dC2P_dalphaS(iprod) = CellOut.dC2P_dalphaS;
            C2Px(:,par.pindx.lalphaS) =  dC2P_dalphaS(iwet) * dalphaS_lalphaS;
			par.CellOut.dC2P_dalphaS = dC2P_dalphaS;
        end
		if par.opt_gammaDNA
            dC2P_dgammaDNA  = M3d*0;
			dC2P_dgammaDNA(iprod)  = CellOut.dC2P_dgammaDNA;
            C2Px(:,par.pindx.tgammaDNA) = dC2P_dgammaDNA(iwet) * dgammaDNA_tgammaDNA;
			par.CellOut.dC2P_dgammaDNA = dC2P_dgammaDNA;
        end
    end


%iwet = par.iwet;

 %%%%%% hessian v2 %%%%%%%%%
 % derivative of cell model w.r.t x (w.r.t. log(param) for most parameters)
    if par.optim == off
        C2Pxx = [] ;
		C2Ppxx = [] ;
    elseif (par.optim & nargout > 2)
		if par.npx > 0
			DIPxx = par.Pxx(1:nwet,:);
		else
			C2Ppxx = [] ;
		end

		% add C2Pxx for P Parameters
	    	%C2Pxx  = zeros(nwet,nbx*nbx);
	    	%%%% recall: C2Px(:,pindx.sigP) = dC2P_dDIP*DIPx(:,pindx.lQ10P)
	   	% DIPxx = par.Pxx(1:nwet,:);
	   	% C2Pxx(:,pindx.lQ10P) = dC2P_dDIP*DIPxx(:,pindx.lQ10P) + d2C2P_dDIP2*DIPx(:,pindx.lQ10P).^2;
		xim = sqrt(-1)*eps^3;
		[CellOut, ~] = CellCHNOP(par,x,P0+xim,N0,T0,Irr0);
		d2C2P_dDIPmolperL = M3d*0;
		d2C2P_dDIPmolperL(iprod) = imag(CellOut.dC2P_dDIP)./eps^3 ; %% this might be a problem when doing complex step test because complex perturbation on DIP and on parameter
		d2C2P_dDIPmolperL(iprod(negDIPindx)) = 0;   % seting deriv of reset (neg) P values to zero. (CellOut.dC2P_dDIP(negDIPindx) = 0;)
		d2C2P_dDIP2 = d2C2P_dDIPmolperL * dDIPmolperL_dDIPmmolperm3^2 ;

		if par.dynamicP == on
			kk = 0;
			for jj = 1:par.npx
				for jk = jj:par.npx
					kk = kk + 1;
					% C2Px = dC2P_dDIP * dDIP_dx1
					% C2Pxx = dC2P_dDIP * d2DIP_dx1_dx2 + d2C2P_dDIP2*dDIP_dx1*dDIP_dx2;
					% C2Pxx(:,jj,jk)
					C2Ppxx(:,kk) = dC2P_dDIP(iwet).*DIPxx(:,kk) + d2C2P_dDIP2(iwet).*DIPx(:,jj).*DIPx(:,jk) ;
				end
			end
		else
			kk = 0;
			d2C2P_dPparam = M3d*0;
			for jj = 1:par.npx
				for jk = jj:par.npx
					kk = kk + 1;
					C2Ppxx(:,kk) = d2C2P_dPparam(iwet) ;
				end
			end
		end % end if par.dynamicP == on

		% add P model Parametrs x Cell model parameters
			% need to extract from CellOut dC2P_dcellmodelparameter (C2Px)
			% d2C2P_dcellparam_dPparam
			% eg  dC2P_dQ10photo * dQ10photo_dlQ10Photo
			% d2C2P_dlQ10Photo_dlDIP = imag(CellOut.dC2P_dQ10Photo)./eps^3 * dQ10photo_dlQ10Photo
			% d2C2P_dlQ10Photo_dPparam = d2C2P_dlQ10Photo_dDIP*dDIPmolperL_dDIPmmolperm3 * DIPx()
			% loop through phosphorus params; inside the loop, will have to do if statements for each cell model parameter
		ck = kk; % ck = sum(1:npx)
		kk = 0;
		if par.dynamicP == on
			for jj=1:par.npx
				if (par.opt_Q10Photo)
					kk = kk + 1;
					d2C2P_dlQ10Photo_dDIP = M3d*0;
					CellOut.dC2P_dQ10Photo(negDIPindx) = 0;   % for reset P values
					d2C2P_dlQ10Photo_dDIP(iprod) = imag(CellOut.dC2P_dQ10Photo)./eps^3 .* dQ10_lQ10Photo;
					C2Ppcellxx(:,kk) = d2C2P_dlQ10Photo_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dlQ10Photo_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj);
				end
				if (par.opt_fStorage)
	                kk = kk + 1;
					d2C2P_dlfStorage_dDIP = M3d*0;
					CellOut.dC2P_dfStorage(negDIPindx) = 0;
					d2C2P_dlfStorage_dDIP(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3 *dfStor_lfStorage ;
					C2Ppcellxx(:,kk) = d2C2P_dlfStorage_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dlfStorage_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj);
				end
				if (par.opt_fRibE)
	                kk = kk + 1;
	                d2C2P_dtfRibE_dDIP = M3d*0;
					CellOut.dC2P_dfRibE(negDIPindx) = 0;
					d2C2P_dtfRibE_dDIP(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3 *dfRibE_tfRibE;
					C2Ppcellxx(:,kk) = d2C2P_dtfRibE_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dtfRibE_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
				end
				if (par.opt_kST0)
	                kk = kk + 1;
	                d2C2P_dlkST0_dDIP = M3d*0;
					CellOut.dC2P_dkST0(negDIPindx) = 0;
					d2C2P_dlkST0_dDIP(iprod) = imag(CellOut.dC2P_dkST0)./eps^3 * dkST0_lkST0;
					C2Ppcellxx(:,kk) = d2C2P_dlkST0_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dlkST0_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
				end
				if (par.opt_PLip_PCutoff)
	                kk = kk + 1;
					d2C2P_dlPCutoff_dDIP = M3d*0;
					CellOut.dC2P_dPCutoff(negDIPindx) = 0;
					d2C2P_dlPCutoff_dDIP(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3 *dPCutoff_lPCutoff;
					C2Ppcellxx(:,kk) = d2C2P_dlPCutoff_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dlPCutoff_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
				end
				if (par.opt_PLip_scale)
				    kk = kk + 1;
					d2C2P_dlPLipscale_dDIP = M3d*0;
					CellOut.dC2P_dPLipscale(negDIPindx) = 0;
					d2C2P_dlPLipscale_dDIP(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3 * dPLipscale_lPLipscale;
					C2Ppcellxx(:,kk) = d2C2P_dlPLipscale_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dlPLipscale_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
				end
				if (par.opt_PStor_rCutoff)
	                kk = kk + 1;
					d2C2P_dlrCutoff_dDIP = M3d*0;
					CellOut.dC2P_drCutoff(negDIPindx) = 0;
					d2C2P_dlrCutoff_dDIP(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3 * drCutoff_lrCutoff;
					C2Ppcellxx(:,kk) = d2C2P_dlrCutoff_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dlrCutoff_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
				end
				if (par.opt_PStor_scale)
	                kk = kk + 1;
					d2C2P_dlPStorscale_dDIP = M3d*0;
					CellOut.dC2P_dPStorscale(negDIPindx) = 0;
					d2C2P_dlPStorscale_dDIP(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3 * dPStorscale_lPStorscale;
					C2Ppcellxx(:,kk) = d2C2P_dlPStorscale_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dlPStorscale_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
				end
				if (par.opt_alphaS)
	                kk = kk + 1;
					d2C2P_dlalphaS_dDIP = M3d*0;
					CellOut.dC2P_dalphaS(negDIPindx) = 0;
					d2C2P_dlalphaS_dDIP(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3 * dalphaS_lalphaS;
					C2Ppcellxx(:,kk) = d2C2P_dlalphaS_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dlalphaS_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
				end
				if (par.opt_gammaDNA)
	                kk = kk + 1;
	                d2C2P_dtgammaDNA_dDIP = M3d*0;
					CellOut.dC2P_dgammaDNA(negDIPindx) = 0;
					d2C2P_dtgammaDNA_dDIP(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3 *dgammaDNA_tgammaDNA;
					C2Ppcellxx(:,kk) = d2C2P_dtgammaDNA_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dtgammaDNA_dDIP(iwet).*dDIPmolperL_dDIPmmolperm3.*DIPx(:,jj) ;
				end
			end % end for jj=1:par.npx
	    else
		% cell model doeos not depend on modelled DIP
			% loop through phosphorus params; set all to zero
			d2C2P_dcell_dPparam = M3d*0;
			for jj=1:par.npx
				if (par.opt_Q10Photo)
					kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_fStorage)
	                kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_fRibE)
	                kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_kST0)
	                kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_PLip_PCutoff)
	                kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_PLip_scale)
					kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_PStor_rCutoff)
					kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_PStor_scale)
					kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_alphaS)
					kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
				if (par.opt_gammaDNA)
					kk = kk + 1;
					C2Ppcellxx(:,kk) = d2C2P_dcell_dPparam(iwet);
					ck = ck + 1;
					C2Ppxx(:,ck) = d2C2P_dcell_dPparam(iwet);
				end
			end % end for jj=1:par.npx
	    end % end if par.dynamicP == on; else

	%----------Cell model parameter----------------
		kk = 0;
        % Q10Photo Derivatives
		if (par.opt_Q10Photo)
            kk = kk + 1;

			%second Derivatives w.r.t. Q10Photo          % change to dlQ10Photo?
			xim = zeros(size(x));
			xim(par.pindx.lQ10Photo) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dQ10Photo2 = M3d*0;
			d2C2P_dQ10Photo2(iprod) = imag(CellOut.dC2P_dQ10Photo)./eps^3;
			%C2Pxx(:,kk) = d2C2P_dQ10Photo2(iwet);
			C2P_Q10_lQ10 = d2C2P_dQ10Photo2(iwet); %complex step takes second deriv wrt log(Q10)
			d2Q10_lQ10Photo = par.BIO.Q10Photo;
			%C2Pxx(:,kk) = (C2Px(:,par.pindx.lQ10Photo)*d2Q10_lQ10Photo + C2P_Q10_lQ10*dQ10_lQ10Photo);
			% need to remove d2Q10_lQ10Photo after updating C2Px ;
			C2Pxx(:,kk) = (C2Px(:,par.pindx.lQ10Photo) + C2P_Q10_lQ10*dQ10_lQ10Photo);

			if (par.opt_fStorage)
                kk = kk + 1;
				d2C2P_dfStorage_dQ10Photo = M3d*0;
				d2C2P_dfStorage_dQ10Photo(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3;
				%C2Pxx(:,kk) = d2C2P_dfStorage_dQ10Photo(iwet);
				C2P_fStor_lQ10 = d2C2P_dfStorage_dQ10Photo(iwet); %d2C2P/(dfStorage dlogQ10photo)
	            C2Pxx(:,kk) = C2P_fStor_lQ10*dfStor_lfStorage;
            end

            if (par.opt_fRibE)
                kk = kk + 1;
                d2C2P_dfRibE_dQ10Photo = M3d*0;
				d2C2P_dfRibE_dQ10Photo(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dfRibE_dQ10Photo(iwet);
				C2P_fRibE_lQ10 = d2C2P_dfRibE_dQ10Photo(iwet);
	            C2Pxx(:,kk) = C2P_fRibE_lQ10 *dfRibE_tfRibE;
            end

			if (par.opt_kST0)
                kk = kk + 1;
                d2C2P_dkST0_dQ10Photo = M3d*0;
				d2C2P_dkST0_dQ10Photo(iprod) = imag(CellOut.dC2P_dkST0)./eps^3;
				%C2Pxx(:,kk) = d2C2P_dkST0_dQ10Photo(iwet);
				C2P_kST0_lQ10 = d2C2P_dkST0_dQ10Photo(iwet);
                C2Pxx(:,kk) = C2P_kST0_lQ10*dkST0_lkST0;
            end

			if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dQ10Photo = M3d*0;
				d2C2P_dPCutoff_dQ10Photo(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dPCutoff_dQ10Photo(iwet);
				C2P_PCutoff_lQ10 = d2C2P_dPCutoff_dQ10Photo(iwet);
				C2Pxx(:,kk) = C2P_PCutoff_lQ10 *dPCutoff_lPCutoff;
            end

			if (par.opt_PLip_scale)
			    kk = kk + 1;
				d2C2P_dPLipscale_dQ10Photo = M3d*0;
				d2C2P_dPLipscale_dQ10Photo(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
			    %C2Pxx(:,kk) = d2C2P_dPLipscale_dQ10Photo(iwet);
				C2P_PLipscale_lQ10 = d2C2P_dPLipscale_dQ10Photo(iwet);
				C2Pxx(:,kk) = C2P_PLipscale_lQ10 *dPLipscale_lPLipscale;
			end

			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dQ10Photo = M3d*0;
				d2C2P_drCutoff_dQ10Photo(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                %C2Pxx(:,kk) = d2C2P_drCutoff_dQ10Photo(iwet);
				C2P_rCutoff_lQ10 = d2C2P_drCutoff_dQ10Photo(iwet);
				C2Pxx(:,kk) = C2P_rCutoff_lQ10 * drCutoff_lrCutoff;
            end

			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dQ10Photo = M3d*0;
				d2C2P_dPStorscale_dQ10Photo(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dPStorscale_dQ10Photo(iwet);
				C2P_PStorscale_lQ10 = d2C2P_dPStorscale_dQ10Photo(iwet);
				C2Pxx(:,kk) = C2P_PStorscale_lQ10 * dPStorscale_lPStorscale;
            end

			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dQ10Photo = M3d*0;
				d2C2P_dalphaS_dQ10Photo(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dalphaS_dQ10Photo(iwet);
				C2P_alphaS_lQ10 = d2C2P_dalphaS_dQ10Photo(iwet);
				C2Pxx(:,kk) = C2P_alphaS_lQ10 * dalphaS_lalphaS;
			end

			if (par.opt_gammaDNA)
                kk = kk + 1;
                d2C2P_dgammaDNA_dQ10Photo = M3d*0;
				d2C2P_dgammaDNA_dQ10Photo(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dgammaDNA_dQ10Photo(iwet);
				C2P_gammaDNA_lQ10 = d2C2P_dgammaDNA_dQ10Photo(iwet);
	            C2Pxx(:,kk) = C2P_gammaDNA_lQ10 *dgammaDNA_tgammaDNA;
            end
        end

        % fStorage
		if (par.opt_fStorage)
			kk = kk + 1;

			%second Derivatives w.r.t. fStorage
			xim = zeros(size(x));
			xim(par.pindx.lfStorage) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dfStorage2 = M3d*0;
			d2C2P_dfStorage2(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3;
            %C2Pxx(:,kk) = d2C2P_dfStorage2(iwet);
			C2P_fStor_lfStor = d2C2P_dfStorage2(iwet);
			d2fStor_lfStorage = par.BIO.fStorage;
			%C2Pxx(:,kk) = C2Px(:,par.pindx.lfStorage)*d2fStor_lfStorage + C2P_fStor_lfStor*dfStor_lfStorage;
			% need to remove d2fStor_lfStorage after updating C2Px ;
			C2Pxx(:,kk) = C2Px(:,par.pindx.lfStorage) + C2P_fStor_lfStor*dfStor_lfStorage;

            if (par.opt_fRibE)
                kk = kk + 1;
                d2C2P_dfRibE_dfStorage = M3d*0;
				d2C2P_dfRibE_dfStorage(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dfRibE_dfStorage(iwet);
				C2P_fRibE_lfStor = d2C2P_dfRibE_dfStorage(iwet);
				C2Pxx(:,kk) = C2P_fRibE_lfStor * dfRibE_tfRibE;
            end
			if (par.opt_kST0)
                kk = kk + 1;
                d2C2P_dkST0_dfStorage = M3d*0;
				d2C2P_dkST0_dfStorage(iprod) = imag(CellOut.dC2P_dkST0)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dkST0_dfStorage(iwet);
				C2P_kST0_lfStor = d2C2P_dkST0_dfStorage(iwet);
				C2Pxx(:,kk) = C2P_kST0_lfStor* dkST0_lkST0;
            end
			if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dfStorage = M3d*0;
				d2C2P_dPCutoff_dfStorage(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dPCutoff_dfStorage(iwet);
				C2P_PCutoff_lfStor = d2C2P_dPCutoff_dfStorage(iwet);
				C2Pxx(:,kk) = C2P_PCutoff_lfStor * dPCutoff_lPCutoff;
			end
			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dfStorage = M3d*0;
				d2C2P_dPLipscale_dfStorage(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dPLipscale_dfStorage(iwet);
				C2P_PLipscale_lfStor = d2C2P_dPLipscale_dfStorage(iwet);
				C2Pxx(:,kk) = C2P_PLipscale_lfStor * dPLipscale_lPLipscale;
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dfStorage = M3d*0;
				d2C2P_drCutoff_dfStorage(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                %C2Pxx(:,kk) = d2C2P_drCutoff_dfStorage(iwet);
				C2P_rCutoff_lfStor = d2C2P_drCutoff_dfStorage(iwet);
				C2Pxx(:,kk) = C2P_rCutoff_lfStor * drCutoff_lrCutoff;
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dfStorage = M3d*0;
				d2C2P_dPStorscale_dfStorage(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dPStorscale_dfStorage(iwet);
				C2P_PStorscale_lfStor = d2C2P_dPStorscale_dfStorage(iwet);
				C2Pxx(:,kk) = C2P_PStorscale_lfStor * dPStorscale_lPStorscale;
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dfStorage = M3d*0;
				d2C2P_dalphaS_dfStorage(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dalphaS_dfStorage(iwet);
				C2P_alphaS_lfStor = d2C2P_dalphaS_dfStorage(iwet);
				C2Pxx(:,kk) = C2P_alphaS_lfStor * dalphaS_lalphaS;
			end
			if (par.opt_gammaDNA)
                kk = kk + 1;
                d2C2P_dgammaDNA_dfStorage = M3d*0;
				d2C2P_dgammaDNA_dfStorage(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dgammaDNA_dfStorage(iwet);
				C2P_gammaDNA_lfStor = d2C2P_dgammaDNA_dfStorage(iwet);
				C2Pxx(:,kk) = C2P_gammaDNA_lfStor * dgammaDNA_tgammaDNA;
            end
        end

        %fRibE derivatives
		if (par.opt_fRibE)
            kk = kk + 1;

			%second Derivatives w.r.t. PLip_PCutoff
			xim = zeros(size(x));
			xim(par.pindx.tfRibE) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dfRibE2 = M3d*0;
			d2C2P_dfRibE2(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
            %C2Pxx(:,kk) = d2C2P_dfRibE2(iwet);
			C2P_fRibE_tfRibE = d2C2P_dfRibE2(iwet);
			% tfRibE = atanh(2*par.BIO.fRibE-1);
			% dfRibE_tfRibE = 0.5*sech(tfRibE)^2;
			d2fRibE_tfRibE = -sech(tfRibE)^2 * tanh(tfRibE);
			% dC2Ptmp = d/dtfRibE[dC2P/dfRibE *dfRibE/dtfRibE]
			% = dC2P/dfRibE * d/dtfRibE[dfRibE/dtfRibE] + d/dtfRibE[dC2P/dfRibE] *dfRibE/dtfRibE
			% C2Pxx(:,kk) = C2Px(:,par.pindx.tfRibE) * d2fRibE_tfRibE + C2P_fRibE_tfRibE * dfRibE_tfRibE;
			% need to remove d2fRibE_tfRibE after updating C2Px? ;
			C2Pxx(:,kk) = C2P_fRibE_tfRibE * dfRibE_tfRibE + C2Px(:,par.pindx.tfRibE) * d2fRibE_tfRibE / dfRibE_tfRibE;

			if (par.opt_kST0)
                kk = kk + 1;
                d2C2P_dkST0_dfRibE = M3d*0;
				d2C2P_dkST0_dfRibE(iprod) = imag(CellOut.dC2P_dkST0)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dkST0_dfRibE(iwet);
				C2Pxx(:,kk) = d2C2P_dkST0_dfRibE(iwet) * dkST0_lkST0;
            end
            if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dfRibE = M3d*0;
				d2C2P_dPCutoff_dfRibE(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dPCutoff_dfRibE(iwet);
				C2P_PCutoff_tfRibE = d2C2P_dPCutoff_dfRibE(iwet);
				C2Pxx(:,kk) = C2P_PCutoff_tfRibE * dPCutoff_lPCutoff;
			end
			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dfRibE = M3d*0;
				d2C2P_dPLipscale_dfRibE(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dPLipscale_dfRibE(iwet);
				C2P_PLipscale_tfRibE = d2C2P_dPLipscale_dfRibE(iwet);
				C2Pxx(:,kk) = C2P_PLipscale_tfRibE * dPLipscale_lPLipscale;
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dfRibE = M3d*0;
				d2C2P_drCutoff_dfRibE(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                %C2Pxx(:,kk) = d2C2P_drCutoff_dfRibE(iwet);
				C2P_rCutoff_tfRibE = d2C2P_drCutoff_dfRibE(iwet);
				C2Pxx(:,kk) = C2P_rCutoff_tfRibE * drCutoff_lrCutoff;
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dfRibE = M3d*0;
				d2C2P_dPStorscale_dfRibE(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dPStorscale_dfRibE(iwet);
				C2P_PStorscale_tfRibE = d2C2P_dPStorscale_dfRibE(iwet);
				C2Pxx(:,kk) = C2P_PStorscale_tfRibE * dPStorscale_lPStorscale;
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dfRibE = M3d*0;
				d2C2P_dalphaS_dfRibE(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dalphaS_dfRibE(iwet);
				C2P_alphaS_tfRibE = d2C2P_dalphaS_dfRibE(iwet);
				C2Pxx(:,kk) = C2P_alphaS_tfRibE * dalphaS_lalphaS;
			end
			if (par.opt_gammaDNA)
                kk = kk + 1;
                d2C2P_dgammaDNA_dfRibE = M3d*0;
				d2C2P_dgammaDNA_dfRibE(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
                %C2Pxx(:,kk) = d2C2P_dgammaDNA_dfRibE(iwet);
				C2P_gammaDNA_tfRibE = d2C2P_dgammaDNA_dfRibE(iwet);
				C2Pxx(:,kk) = C2P_gammaDNA_tfRibE * dgammaDNA_tgammaDNA;
            end
        end

		%kST0 derivatives
		if (par.opt_kST0)
            kk = kk + 1;

			%second Derivatives w.r.t. kST0
			xim = zeros(size(x));
			xim(par.pindx.lkST0) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			%d2C2P_dkST02 = M3d*0;
			%d2C2P_dkST02(iprod) = imag(CellOut.dC2P_dkST0)./eps^3;
            %C2Pxx(:,kk) = d2C2P_dkST02(iwet);

			C2P_kST0_lkST0 = M3d*0;
			C2P_kST0_lkST0(iprod) = imag(CellOut.dC2P_dkST0)./eps^3;
			d2kST0_lkST0 = par.BIO.kST0;
			% C2Pxx(:,kk) = C2Px(:,par.pindx.lkST0) * d2kST0_lkST0 + C2P_kST0_lkST0(iwet) * dkST0_lkST0;
			% need to remove d2kST0_lkST0 after updating C2Px? ;
			C2Pxx(:,kk) = C2Px(:,par.pindx.lkST0) + C2P_kST0_lkST0(iwet) * dkST0_lkST0;

            if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				C2P_PCutoff_lkST0 = M3d*0;
				C2P_PCutoff_lkST0(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                C2Pxx(:,kk) = C2P_PCutoff_lkST0(iwet) * dPCutoff_lPCutoff;
			end
			if (par.opt_PLip_scale)
                kk = kk + 1;
				C2P_PLipscale_lkST0 = M3d*0;
				C2P_PLipscale_lkST0(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = C2P_PLipscale_lkST0(iwet) * dPLipscale_lPLipscale;
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				C2P_rCutoff_lkST0 = M3d*0;
				C2P_rCutoff_lkST0(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = C2P_rCutoff_lkST0(iwet) * drCutoff_lrCutoff;
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				C2P_PStorscale_lkST0 = M3d*0;
				C2P_PStorscale_lkST0(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = C2P_PStorscale_lkST0(iwet) * dPStorscale_lPStorscale;
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				C2P_alphaS_lkST0 = M3d*0;
				C2P_alphaS_lkST0(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = C2P_alphaS_lkST0(iwet) * dalphaS_lalphaS;
			end
			if (par.opt_gammaDNA)
                kk = kk + 1;
				C2P_gammaDNA_lkST0 = M3d*0;
				C2P_gammaDNA_lkST0(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
	            C2Pxx(:,kk) = C2P_gammaDNA_lkST0(iwet) * dgammaDNA_tgammaDNA;
            end
        end

        %PLip_PCutoff
        % PLip_PCutoff derivatives
        if (par.opt_PLip_PCutoff)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPLip_PCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			%d2C2P_dPCutoff2 = M3d*0;
			%d2C2P_dPCutoff2(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
            %C2Pxx(:,kk) = d2C2P_dPCutoff2(iwet);
			C2P_PCutoff_lPCutoff = M3d*0;
			C2P_PCutoff_lPCutoff(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
			d2PCutoff_lPCutoff = par.BIO.PLip_PCutoff;
            % C2Pxx(:,kk) = C2P_PCutoff_lPCutoff(iwet) * dPCutoff_lPCutoff + C2Px(:,par.pindx.lPLip_PCutoff) * d2PCutoff_lPCutoff;
			% need to remove d2PCutoff_lPCutoff after updating C2Px? ;
			C2Pxx(:,kk) = C2P_PCutoff_lPCutoff(iwet) * dPCutoff_lPCutoff + C2Px(:,par.pindx.lPLip_PCutoff);

			if (par.opt_PLip_scale)
                kk = kk + 1;
				C2P_PLipscale_lPCutoff = M3d*0;
				C2P_PLipscale_lPCutoff(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = C2P_PLipscale_lPCutoff(iwet) * dPLipscale_lPLipscale;
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				C2P_rCutoff_lPCutoff = M3d*0;
				C2P_rCutoff_lPCutoff(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = C2P_rCutoff_lPCutoff(iwet) * drCutoff_lrCutoff;
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				C2P_PStorscale_lPCutoff = M3d*0;
				C2P_PStorscale_lPCutoff(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = C2P_PStorscale_lPCutoff(iwet) * dPStorscale_lPStorscale;
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				C2P_alphaS_lPCutoff = M3d*0;
				C2P_alphaS_lPCutoff(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = C2P_alphaS_lPCutoff(iwet) * dalphaS_lalphaS;
            end
			if (par.opt_gammaDNA)
                kk = kk + 1;
				C2P_gammaDNA_lPCutoff = M3d*0;
				C2P_gammaDNA_lPCutoff(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
	            C2Pxx(:,kk) = C2P_gammaDNA_lPCutoff(iwet) * dgammaDNA_tgammaDNA;
            end
        end

        %PLip_scale
        if (par.opt_PLip_scale)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPLip_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			C2P_PLipscale_lPLipscale = M3d*0;
			C2P_PLipscale_lPLipscale(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
			d2PLipscale_lPLipscale = par.BIO.PLip_scale;
            % C2Pxx(:,kk) = C2P_PLipscale_lPLipscale(iwet) * dPLipscale_lPLipscale + C2Px(:,par.pindx.lPLip_scale)*d2PLipscale_lPLipscale;
			% need to remove d2PLipscale_lPLipscale after updating C2Px? ;
			C2Pxx(:,kk) = C2P_PLipscale_lPLipscale(iwet) * dPLipscale_lPLipscale + C2Px(:,par.pindx.lPLip_scale);

			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				C2P_rCutoff_lPLipscale = M3d*0;
				C2P_rCutoff_lPLipscale(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = C2P_rCutoff_lPLipscale(iwet) * drCutoff_lrCutoff;
			end
			if (par.opt_PStor_scale)
				kk = kk + 1;
                C2P_PStorscale_lPLipscale = M3d*0;
				C2P_PStorscale_lPLipscale(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = C2P_PStorscale_lPLipscale(iwet) * dPStorscale_lPStorscale;
			end
			if (par.opt_alphaS)
				kk = kk + 1;
                C2P_alphaS_lPLipscale = M3d*0;
				C2P_alphaS_lPLipscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = C2P_alphaS_lPLipscale(iwet) * dalphaS_lalphaS;
            end
			if (par.opt_gammaDNA)
                kk = kk + 1;
				C2P_gammaDNA_lPLipscale = M3d*0;
				C2P_gammaDNA_lPLipscale(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
	            C2Pxx(:,kk) = C2P_gammaDNA_lPLipscale(iwet) * dgammaDNA_tgammaDNA;
            end
        end

        %PStor_rCutoff
        if (par.opt_PStor_rCutoff)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPStor_rCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			C2P_rCutoff_lrCutoff = M3d*0;
			C2P_rCutoff_lrCutoff(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
			d2rCutoff_lrCutoff = par.BIO.PStor_rCutoff;
            % C2Pxx(:,kk) = C2P_rCutoff_lrCutoff(iwet) * drCutoff_lrCutoff + C2Px(:,par.pindx.lPStor_rCutoff)*d2rCutoff_lrCutoff;
			% need to remove d2rCutoff_lrCutoff after updating C2Px ;
			C2Pxx(:,kk) = C2P_rCutoff_lrCutoff(iwet) * drCutoff_lrCutoff + C2Px(:,par.pindx.lPStor_rCutoff);

			if (par.opt_PStor_scale)
                kk = kk + 1;
				C2P_PStorscale_lrCutoff = M3d*0;
				C2P_PStorscale_lrCutoff(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = C2P_PStorscale_lrCutoff(iwet) * dPStorscale_lPStorscale;
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				C2P_alphaS_lrCutoff = M3d*0;
				C2P_alphaS_lrCutoff(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = C2P_alphaS_lrCutoff(iwet) * dalphaS_lalphaS;
            end
			if (par.opt_gammaDNA)
                kk = kk + 1;
				C2P_gammaDNA_lrCutoff = M3d*0;
				C2P_gammaDNA_lrCutoff(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
	            C2Pxx(:,kk) = C2P_gammaDNA_lrCutoff(iwet) * dgammaDNA_tgammaDNA;
            end
        end

        %PStor_scale
        if (par.opt_PStor_scale)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPStor_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			C2P_PStorscale_lPStorscale = M3d*0;
			C2P_PStorscale_lPStorscale(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			d2PStorscale_lPStorscale = par.BIO.PStor_scale;
            % C2Pxx(:,kk) = C2P_PStorscale_lPStorscale(iwet)*dPStorscale_lPStorscale + C2Px(:,par.pindx.lPStor_scale)*d2PStorscale_lPStorscale;
			% need to remove d2PLipscale_lPLipscale after updating C2Px? ;
			C2Pxx(:,kk) = C2P_PStorscale_lPStorscale(iwet)*dPStorscale_lPStorscale + C2Px(:,par.pindx.lPStor_scale);

			if (par.opt_alphaS)
                kk = kk + 1;
				C2P_alphaS_lPStorscale = M3d*0;
				C2P_alphaS_lPStorscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = C2P_alphaS_lPStorscale(iwet) * dalphaS_lalphaS;
			end
			if (par.opt_gammaDNA)
                kk = kk + 1;
				C2P_gammaDNA_lPStorscale = M3d*0;
				C2P_gammaDNA_lPStorscale(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
	            C2Pxx(:,kk) = C2P_gammaDNA_lPStorscale(iwet) * dgammaDNA_tgammaDNA;
            end
        end

        % alphaS
        if (par.opt_alphaS)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lalphaS) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			C2P_alphaS_lalphaS = M3d*0;
			C2P_alphaS_lalphaS(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			d2alphaS_lalphaS = par.BIO.alphaS;
            % C2Pxx(:,kk) = C2P_alphaS_lalphaS(iwet)*dalphaS_lalphaS + C2Px(:,par.pindx.lalphaS)*d2alphaS_lalphaS;
			% need to remove d2alphaS_lalphaS after updating C2Px? ;
			C2Pxx(:,kk) = C2P_alphaS_lalphaS(iwet)*dalphaS_lalphaS + C2Px(:,par.pindx.lalphaS);

			if (par.opt_gammaDNA)
                kk = kk + 1;
				C2P_gammaDNA_lalphaS = M3d*0;
				C2P_gammaDNA_lalphaS(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
	            C2Pxx(:,kk) = C2P_gammaDNA_lalphaS(iwet) * dgammaDNA_tgammaDNA;
            end
        end

		if (par.opt_gammaDNA)
            kk = kk + 1;

			%second Derivatives w.r.t. gammaDNA
			xim = zeros(size(x));
			xim(par.pindx.tgammaDNA) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

			C2P_gammaDNA_tgammaDNA = M3d*0;
			C2P_gammaDNA_tgammaDNA(iprod) = imag(CellOut.dC2P_dgammaDNA)./eps^3;
			%dgammaDNA_dtgammaDNA  = 0.5*sech(tgammaDNA)^2;
			d2gammaDNA_dtgammaDNA = -sech(tgammaDNA)^2 * tanh(tgammaDNA);
			C2Pxx(:,kk) = C2Px(:,par.pindx.tgammaDNA)./dgammaDNA_tgammaDNA.*d2gammaDNA_dtgammaDNA + C2P_gammaDNA_tgammaDNA(iwet)*dgammaDNA_tgammaDNA;
		end

    end %if par.optim
end % if TestVer == 2

if TestVer == 3 %trying to index explicitly in 2 dimensions instead of using kk
	C2Pxx = zeros(nwet,npx+ncx+nbx,npx+ncx+nbx);

	%-----P model Derivatives------------
	DIPxx = par.Pxx(1:nwet,:);

	xim = sqrt(-1)*eps^3;
	[CellOut, ~] = CellCHNOP(par,x,P0+xim,N0,T0,Irr0);
	d2C2P_dDIPmolperL = M3d*0;
	d2C2P_dDIPmolperL(iprod) = imag(CellOut.dC2P_dDIP)./eps^3 ;
	d2C2P_dDIP2 = dC2P_dDIPmolperL * dDIPmolperL_dDIPmmolperm3^2 ;

	kk = 0;
	for jj = 1:npx
		for jk = jj:npx
			kk = kk + 1;
			%C2Pxx = dC2P_dDIP * d2DIP_dx1_dx2 + d2C2P_dDIP2*dDIP_dx1*dDIP_dx2;
			% C2Pxx(:,jj,jk)
			C2Ppxx(:,jj,jk) = dC2P_dDIP(iwet).*DIPxx(:,kk) + d2C2P_dDIP2(iwet).*DIPx(:,jj).*DIPx(:,jk) ;
		end
	end

	% Q10Photo Derivatives
	kk=0;
	if (par.opt_Q10Photo)
		kk = kk + 1;

		%second Derivatives w.r.t. Q10Photo          % change to dlQ10Photo?
		xim = zeros(size(x));
		xim(par.pindx.lQ10Photo) = sqrt(-1)*eps^3;
		[CellOut, ~] = CellCHNOP(par,x+xim, P0,N0,T0,Irr0);

		d2C2P_dQ10Photo2 = M3d*0;
		d2C2P_dQ10Photo2(iprod) = imag(CellOut.dC2P_dQ10Photo)./eps^3;
		%C2Pxx(:,kk) = d2C2P_dQ10Photo2(iwet);
		C2P_Q10_lQ10 = d2C2P_dQ10Photo2(iwet); %complex step takes second deriv wrt log(Q10)
		d2Q10_lQ10Photo = par.BIO.Q10Photo;
		C2Pxx(:,pindx.lQ10Photo,pindx.lQ10Photo) = (C2Px(:,pindx.lQ10Photo)*d2Q10_lQ10Photo + C2P_Q10_lQ10*dQ10_lQ10Photo);
	end %if (par.opt_Q10Photo)


end %if TestVer == 3

end % end function
