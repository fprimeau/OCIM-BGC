%CellUptake
function [C2P, C2Px, C2Pxx] = eqC2Puptake(par)
    on = true; off = false;
    % testing versions
    TestVer = 1 %structure with all derivs as fields
    %TestVer = 2; %matrix derivative matrix

    M3d =par.M3d;
    pindx = par.pindx;
    % unpack the parameters to be optimized
    alpha = par.alpha;
    beta  = par.beta;

    iwet = par.iwet;
    nwet = par.nwet;


if (par.Cellmodel == on)
		iprod = find(M3d(:,:,1:2)); %production in top two layers
		P0 = par.DIPgrd(iprod)./10^6; 			% convert ug/m^3 to g/m^3
		N0 = par.no3obs(iprod)./10^6;   % convert ug/m^3 to g/m^3
		T0 = par.Tobs(iprod);
		Irr0 = par.PARobs(iprod);

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

        C2P=par.CellOut.C2P;

if TestVer == 1
		% Q10Photo Derivatives
		if (par.opt_Q10Photo)
			par.CellOut.dC2P_dQ10Photo = M3d*0;
			par.CellOut.dC2P_dQ10Photo(iprod) = real(CellOut.dC2P_dQ10Photo);

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

			if (par.opt_Q10Photo)
				par.CellOut.d2C2P_dQ10Photo_dfStorage = M3d*0;
				par.CellOut.d2C2P_dQ10Photo_dfStorage(iprod) = imag(CellOut.dC2P_dQ10Photo)./eps^3;
            end
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

		data.CellOut = par.CellOut;

end
end


if par.Cellmodel == off
	DIP = par.DIP ;
	C2P = 1./(cc*DIP + dd);
end

    % of uptake operator
    po4obs = par.po4obs(iwet);
    % P uptake operator
    L = par.L;

    DIP = par.DIP ;
    G   = d0(alpha*L*DIP) ;
    Gp  = alpha*L ;


%% %%% Gradient v1
if TestVer == 1
    if par.optim == off
        C2Px = [];
    elseif (par.optim & nargout > 1)
        % gradient of uptake operator
        nbx  = par.nbx; ncx=par.ncx; npx =par.npx;
        C2Px  = zeros(nwet,npx+ncx+nbx);
        %DIPx = par.Px(1:nwet,:);

        if (par.opt_Q10Photo)
            C2Px(:,par.pindx.lQ10Photo) = par.CellOut.dC2P_dQ10Photo(iwet);
        end
        if (par.opt_fStorage)
            C2Px(:,par.pindx.lfStorage) = par.CellOut.dC2P_dfStorage(iwet);
        end
        if par.opt_fRibE
            C2Px(:,par.pindx.tfRibE) = par.CellOut.dC2P_dfRibE(iwet);
        end
		if (par.opt_PLip_PCutoff)
			C2Px(:,par.pindx.lPLip_PCutoff) = par.CellOut.dC2P_dPCutoff(iwet);
        end
        if (par.opt_PLip_scale)
			C2Px(:,par.pindx.lPLip_scale) = par.CellOut.dC2P_dPLip_scale(iwet);
		end
		if (par.opt_PStor_rCutoff)
			C2Px(:,par.pindx.lPStor_rCutoff) = par.CellOut.dC2P_drCutoff(iwet);
        end
        if par.opt_PStor_scale
			C2Px(:,par.pindx.lPStor_scale) = par.CellOut.dC2P_dPStor_scale(iwet);
        end
        if par.alphaS == on
            C2Px(:,par.pindx.lalphaS) =  par.CellOut.dC2P_dalphaS(iwet);
        end
    end
   %%% Hessian
    if par.optim == off
        C2Pxx = [];
    elseif (par.optim & nargout > 2)
        kk = 0;
        %DIPxx = par.Pxx(1:nwet,:);

        % Q10Photo Q10Photo
        if (par.opt_Q10Photo)
            kk = kk + 1;
            C2Pxx(:,kk) = par.CellOut.d2C2P_dQ10Photo2(iwet);
        end

        % Q10Photo fStorage
        if (par.opt_Q10Photo & par.opt_fStorage)
            kk = kk + 1;
            C2Pxx(:,kk) = par.CellOut.d2C2P_dfStorage_dQ10Photo(iwet);
        end
    end

end

%% ------------------------------------------------
%%%%% Gradient v2
if TestVer == 2

    if par.optim == off
        C2Px = [];
    elseif (par.optim & nargout > 1)
        % gradient of uptake operator
        nbx  = par.nbx; ncx=par.ncx; npx =par.npx;
        C2Px  = zeros(nwet,npx+ncx+nbx);
        %DIPx = par.Px(1:nwet,:);

        if (par.opt_Q10Photo)
            dC2P_dQ10Photo = M3d*0;
			dC2P_dQ10Photo(iprod) = CellOut.dC2P_dQ10Photo;
            C2Px(:,par.pindx.lQ10Photo) = dC2P_dQ10Photo(iwet);
        end
        if (par.opt_fStorage)
            dC2P_dfStorage = M3d*0;
			dC2P_dfStorage(iprod) = CellOut.dC2P_dfStorage;
            C2Px(:,par.pindx.lfStorage) = dC2P_dfStorage(iwet);
        end
        if par.opt_fRibE
            dC2P_dfRibE  = M3d*0;
			dC2P_dfRibE(iprod)  = CellOut.dC2P_dfRibE;
            C2Px(:,par.pindx.tfRibE) = dC2P_dfRibE(iwet);
        end
		if (par.opt_PLip_PCutoff)
            dC2P_dPCutoff  = M3d*0;
			dC2P_dPCutoff(iprod)  = CellOut.dC2P_dPCutoff;
			C2Px(:,par.pindx.lPLip_PCutoff) = dC2P_dPCutoff(iwet);
        end
        if (par.opt_PLip_scale)
            dC2P_dPLip_scale = M3d*0;
			dC2P_dPLip_scale(iprod) = CellOut.dC2P_dPLipscale;
			C2Px(:,par.pindx.lPLip_scale) = dC2P_dPLip_scale(iwet);
		end
		if (par.opt_PStor_rCutoff)
            dC2P_drCutoff  = M3d*0;
			dC2P_drCutoff(iprod)  = CellOut.dC2P_drCutoff;
			C2Px(:,par.pindx.lPStor_rCutoff) = dC2P_drCutoff(iwet);
        end
        if par.opt_PStor_scale
            dC2P_dPStor_scale = M3d*0;
			dC2P_dPStor_scale(iprod) = CellOut.dC2P_dPStorscale;
			C2Px(:,par.pindx.lPStor_scale) = dC2P_dPStor_scale(iwet);
        end
        if par.alphaS
            dC2P_dalphaS = M3d*0;
			dC2P_dalphaS(iprod) = CellOut.dC2P_dalphaS;
            C2Px(:,par.pindx.lalphaS) =  dC2P_dalphaS(iwet);
        end
    end


 %%%%%% hessian v2 %%%%%%%%%
    if par.optim == off
        C2Pxx = [];
    elseif (par.optim & nargout > 2)
        kk = 0;

        % Q10Photo Derivatives
		if (par.opt_Q10Photo)
            kk = kk + 1;

			%second Derivatives w.r.t. Q10Photo          % change to dlQ10Photo?
			xim = zeros(size(x));
			xim(par.pindx.lQ10Photo) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dQ10Photo2 = M3d*0;
			d2C2P_dQ10Photo2(iprod) = imag(CellOut.dC2P_dQ10Photo)./eps^3;
            C2Pxx(:,kk) = d2C2P_dQ10Photo2(iwet);

			if (par.opt_fStorage)
                kk = kk + 1;
				d2C2P_dfStorage_dQ10Photo = M3d*0;
				d2C2P_dfStorage_dQ10Photo(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3;
                C2Pxx(:,kk) = d2C2P_dfStorage_dQ10Photo(iwet);
            end

            if (par.opt_fRibE)
                kk = kk + 1;
                d2C2P_dfRibE_dQ10Photo = M3d*0;
				d2C2P_dfRibE_dQ10Photo(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
                C2Pxx(:,kk) = d2C2P_dfRibE_dQ10Photo(iwet);
            end

			if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dQ10Photo = M3d*0;
				d2C2P_dPCutoff_dQ10Photo(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPCutoff_dQ10Photo(iwet);
            end

			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dQ10Photo = M3d*0;
				d2C2P_drCutoff_dQ10Photo(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dQ10Photo(iwet);
            end

			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dQ10Photo = M3d*0;
				d2C2P_dPStorscale_dQ10Photo(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dQ10Photo(iwet);
            end

			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dQ10Photo = M3d*0;
				d2C2P_dPLipscale_dQ10Photo(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPLipscale_dQ10Photo(iwet);
            end

			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dQ10Photo = M3d*0;
				d2C2P_dalphaS_dQ10Photo(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dQ10Photo(iwet);
			end
        end

        % fStorage
		if (par.opt_fStorage)
			kk = kk + 1;

			%second Derivatives w.r.t. fStorage
			xim = zeros(size(x));
			xim(par.pindx.lfStorage) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dfStorage2 = M3d*0;
			d2C2P_dfStorage2(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3;
            C2Pxx(:,kk) = d2C2P_dfStorage2(iwet);

            if (par.opt_fRibE)
                kk = kk + 1;
                d2C2P_dfRibE_dfStorage = M3d*0;
				d2C2P_dfRibE_dfStorage(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
                C2Pxx(:,kk) = d2C2P_dfRibE_dfStorage(iwet);
            end
			if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dfStorage = M3d*0;
				d2C2P_dPCutoff_dfStorage(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPCutoff_dfStorage(iwet);
			end
			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dfStorage = M3d*0;
				d2C2P_dPLipscale_dfStorage(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPLipscale_dfStorage(iwet);
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dfStorage = M3d*0;
				d2C2P_drCutoff_dfStorage(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dfStorage(iwet);
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dfStorage = M3d*0;
				d2C2P_dPStorscale_dfStorage(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dfStorage(iwet);
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dfStorage = M3d*0;
				d2C2P_dalphaS_dfStorage(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dfStorage(iwet);
			end
        end

        %fRibE derivatives
		if (par.opt_fRibE)
            kk = kk + 1;

			%second Derivatives w.r.t. PLip_PCutoff
			xim = zeros(size(x));
			xim(par.pindx.tfRibE) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dfRibE2 = M3d*0;
			d2C2P_dfRibE2(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
            C2Pxx(:,kk) = d2C2P_dfRibE2(iwet);

            if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dfRibE = M3d*0;
				d2C2P_dPCutoff_dfRibE(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPCutoff_dfRibE(iwet);
			end
			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dfRibE = M3d*0;
				d2C2P_dPLipscale_dfRibE(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPLipscale_dfRibE(iwet);
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dfRibE = M3d*0;
				d2C2P_drCutoff_dfRibE(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dfRibE(iwet);
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dfRibE = M3d*0;
				d2C2P_dPStorscale_dfRibE(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dfRibE(iwet);
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dfRibE = M3d*0;
				d2C2P_dalphaS_dfRibE(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dfRibE(iwet);
			end
        end

        %PLip_PCutoff
        % PLip_PCutoff derivatives
        if (par.opt_PLip_PCutoff)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPLip_PCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dPCutoff2 = M3d*0;
			d2C2P_dPCutoff2(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
            C2Pxx(:,kk) = d2C2P_dPCutoff2(iwet);

			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dPCutoff = M3d*0;
				d2C2P_dPLipscale_dPCutoff(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPLipscale_dPCutoff(iwet);
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dPCutoff = M3d*0;
				d2C2P_drCutoff_dPCutoff(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dPCutoff(iwet);
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dPCutoff = M3d*0;
				d2C2P_dPStorscale_dPCutoff(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dPCutoff(iwet);
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dPCutoff = M3d*0;
				d2C2P_dalphaS_dPCutoff(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dPCutoff(iwet);
            end
        end

        %PLip_scale
        if (par.opt_PLip_scale)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPLip_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dPLipscale2 = M3d*0;
			d2C2P_dPLipscale2(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
            C2Pxx(:,kk) = d2C2P_dPLipscale2(iwet);

			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dPLipscale = M3d*0;
				d2C2P_drCutoff_dPLipscale(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dPLipscale(iwet);
			end
			if (par.opt_PStor_scale)
				kk = kk + 1;
                d2C2P_dPStorscale_dPLipscale = M3d*0;
				d2C2P_dPStorscale_dPLipscale(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dPLipscale(iwet);
			end
			if (par.opt_alphaS)
				kk = kk + 1;
                d2C2P_dalphaS_dPLipscale = M3d*0;
				d2C2P_dalphaS_dPLipscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dPLipscale(iwet);
            end
        end

        %PStor_rCutoff
        if (par.opt_PStor_rCutoff)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPStor_rCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_drCutoff2 = M3d*0;
			d2C2P_drCutoff2(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
            C2Pxx(:,kk) = d2C2P_drCutoff2(iwet);

			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_drCutoff = M3d*0;
				d2C2P_dPStorscale_drCutoff(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_drCutoff(iwet);
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_drCutoff = M3d*0;
				d2C2P_dalphaS_drCutoff(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_drCutoff(iwet);
            end
        end

        %PStor_scale
        if (par.opt_PStor_scale)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPStor_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dPStorscale2 = M3d*0;
			d2C2P_dPStorscale2(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
            C2Pxx(:,kk) = d2C2P_dPStorscale2(iwet);

			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dPStorscale = M3d*0;
				d2C2P_dalphaS_dPStorscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dPStorscale(iwet);
			end
        end

        % alphaS
        if (par.opt_alphaS)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lalphaS) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dalphaS2 = M3d*0;
			d2C2P_dalphaS2(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
            C2Pxx(:,kk) = d2C2P_dalphaS2(iwet);
        end

    end
end

end
