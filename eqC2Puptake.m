%CellUptake
function [par, C2P, C2Px, C2Pxx] = eqC2Puptake(x, par, data)
    on = true; off = false;
    % testing versions
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
	%     % of uptake operator
	%     po4obs = par.po4obs(iwet);
	%     % P uptake operator
	%     L = par.L;
	%
	%     DIP = par.DIP ;
	%     G   = d0(alpha*L*DIP) ;
	%     Gp  = alpha*L ;


		iprod = find(M3d(:,:,1:2)); %production in top two layers
		P0 = data.DIP(iprod)./10^6;		% convert ug/m^3 to g/m^3
		N0 = par.no3obs(iprod)./10^6;   % convert ug/m^3 to g/m^3
		T0 = par.Temp(iprod);
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

        C2P = par.CellOut.C2P;
		%data.CellOut = par.CellOut;

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
			if (par.opt_kST0)
				par.CellOut.d2C2P_dkST0_dQ10Photo = M3d*0;
				par.CellOut.d2C2P_dkST0_dQ10Photo(iprod) = imag(CellOut.dC2P_kST0)./eps^3;
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
			if (par.opt_kST0)
				par.CellOut.d2C2P_dkST0_dfStorage = M3d*0;
				par.CellOut.d2C2P_dkST0_dfStorage(iprod) = imag(CellOut.dC2P_kST0)./eps^3;
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

			if (par.opt_kST0)
				par.CellOut.d2C2P_dkST0_dfRibE = M3d*0;
				par.CellOut.d2C2P_dkST0_dfRibE(iprod) = imag(CellOut.dC2P_kST0)./eps^3;
			end
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

		%kST0 derivatives
		if (par.opt_kST0)
			par.CellOut.dC2P_dkST0  = M3d*0;
			par.CellOut.dC2P_dkST0(iprod)  = real(CellOut.dC2P_dkST0);

			%second Derivatives w.r.t. kST0
			xim = zeros(size(x));
			xim(par.pindx.lkST0) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			par.CellOut.d2C2P_dkST02 = M3d*0;
			par.CellOut.d2C2P_dkST02(iprod) = imag(CellOut.dC2P_dkST0)./eps^3;

            if (par.opt_PLip_PCutoff)
				par.CellOut.d2C2P_dPCutoff_dkST0 = M3d*0;
				par.CellOut.d2C2P_dPCutoff_dkST0(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
			end
			if (par.opt_PLip_scale)
				par.CellOut.d2C2P_dPLipscale_dkST0 = M3d*0;
				par.CellOut.d2C2P_dPLipscale_dkST0(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
			end
			if (par.opt_PStor_rCutoff)
				par.CellOut.d2C2P_drCutoff_dkST0 = M3d*0;
				par.CellOut.d2C2P_drCutoff_dkST0(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
			end
			if (par.opt_PStor_scale)
				par.CellOut.d2C2P_dPStorscale_dkST0 = M3d*0;
				par.CellOut.d2C2P_dPStorscale_dkST0(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			end
			if (par.opt_alphaS)
				par.CellOut.d2C2P_dalphaS_dkST0 = M3d*0;
				par.CellOut.d2C2P_dalphaS_dkST0(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
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
end


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
		if par.opt_kST0
            C2Px(:,par.pindx.lkST0) = par.CellOut.dC2P_dkST0(iwet);
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
			par.CellOut.dC2P_dQ10Photo = dC2P_dQ10Photo;
        end
        if (par.opt_fStorage)
            dC2P_dfStorage = M3d*0;
			dC2P_dfStorage(iprod) = CellOut.dC2P_dfStorage;
            C2Px(:,par.pindx.lfStorage) = dC2P_dfStorage(iwet);
			par.CellOut.dC2P_dfStorage = dC2P_dfStorage;
        end
        if par.opt_fRibE
            dC2P_dfRibE  = M3d*0;
			dC2P_dfRibE(iprod)  = CellOut.dC2P_dfRibE;
            C2Px(:,par.pindx.tfRibE) = dC2P_dfRibE(iwet);
			par.CellOut.dC2P_dfRibE = dC2P_dfRibE;
        end
		if par.opt_kST0
            dC2P_dkST0  = M3d*0;
			dC2P_dkST0(iprod)  = CellOut.dC2P_dkST0;
            C2Px(:,par.pindx.lkST0) = dC2P_dkST0(iwet);
			par.CellOut.dC2P_dkST0 = dC2P_kST0;
        end
		if (par.opt_PLip_PCutoff)
            dC2P_dPCutoff  = M3d*0;
			dC2P_dPCutoff(iprod)  = CellOut.dC2P_dPCutoff;
			C2Px(:,par.pindx.lPLip_PCutoff) = dC2P_dPCutoff(iwet);
			par.CellOut.dC2P_dPCutoff = dC2P_dPCutoff;
        end
        if (par.opt_PLip_scale)
            dC2P_dPLip_scale = M3d*0;
			dC2P_dPLip_scale(iprod) = CellOut.dC2P_dPLipscale;
			C2Px(:,par.pindx.lPLip_scale) = dC2P_dPLip_scale(iwet);
			par.CellOut.dC2P_dPLip_scale = dC2P_dPLip_scale;
		end
		if (par.opt_PStor_rCutoff)
            dC2P_drCutoff  = M3d*0;
			dC2P_drCutoff(iprod)  = CellOut.dC2P_drCutoff;
			C2Px(:,par.pindx.lPStor_rCutoff) = dC2P_drCutoff(iwet);
			par.CellOut.dC2P_drCutoff = dC2P_drCutoff;
        end
        if par.opt_PStor_scale
            dC2P_dPStor_scale = M3d*0;
			dC2P_dPStor_scale(iprod) = CellOut.dC2P_dPStorscale;
			C2Px(:,par.pindx.lPStor_scale) = dC2P_dPStor_scale(iwet);
			par.CellOut.dC2P_dPStor_scale = dC2P_dPStor_scale;
        end
        if par.opt_alphaS
            dC2P_dalphaS = M3d*0;
			dC2P_dalphaS(iprod) = CellOut.dC2P_dalphaS;
            C2Px(:,par.pindx.lalphaS) =  dC2P_dalphaS(iwet);
			par.CellOut.dC2P_dalphaS = dC2P_dalphaS;
        end
    end

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
%iwet = par.iwet;

 %%%%%% hessian v2 %%%%%%%%%
 % derivative of cell model w.r.t x (w.r.t. log(param) for most parameters)
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
			%C2Pxx(:,kk) = d2C2P_dQ10Photo2(iwet);
			C2P_Q10_lQ10 = d2C2P_dQ10Photo2(iwet); %complex step takes second deriv wrt log(Q10)
			d2Q10_lQ10Photo = par.BIO.Q10Photo;
			C2Pxx(:,kk) = (C2Px(:,par.pindx.lQ10Photo)*d2Q10_lQ10Photo + C2P_Q10_lQ10*dQ10_lQ10Photo);


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
				C2P_kST0_lQ10 = d2C2P_dkST0_dQ10Photo(iwet)
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
            %C2Pxx(:,kk) = d2C2P_dfStorage2(iwet);
			C2P_fStor_lfStor = d2C2P_dfStorage2(iwet);
			d2fStor_lfStorage = par.BIO.fStorage;
			C2Pxx(:,kk) = C2Px(:,par.pindx.lfStorage)*d2fStor_lfStorage + C2P_fStor_lfStor*dfStor_lfStorage;

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
            %C2Pxx(:,kk) = d2C2P_dfRibE2(iwet);
			C2P_fRibE_tfRibE = d2C2P_dfRibE2(iwet);
			d2fRibE_tfRibE = -sech(tfRibE)^2 * tanh(tfRibE);
			% dC2Ptmp = d/dtfRibE[dC2P/dfRibE *dfRibE/dtfRibE]
			% = dC2P/dfRibE * d/dtfRibE[dfRibE/dtfRibE] + d/dtfRibE[dC2P/dfRibE] *dfRibE/dtfRibE
			C2Pxx(:,kk) = C2Px(:,par.pindx.tfRibE) * d2fRibE_tfRibE + C2P_fRibE_tfRibE * dfRibE_tfRibE;

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
        end

		%kST0 derivatives
		if (par.opt_kST0)
            kk = kk + 1;

			%second Derivatives w.r.t. kST0
			xim = zeros(size(x));
			xim(par.pindx.lkST0) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			%d2C2P_dkST02 = M3d*0;
			%d2C2P_dkST02(iprod) = imag(CellOut.dC2P_dkST0)./eps^3;
            %C2Pxx(:,kk) = d2C2P_dkST02(iwet);

			C2P_kST0_lkST0 = M3d*0;
			C2P_kST0_lkST0(iprod) = imag(CellOut.dC2P_dkST0)./eps^3;
			d2kST0_lkST0 = par.BIO.kST0;
			C2Pxx(:,kk) = C2Px(:,par.pindx.lkST0) * d2kST0_lkST0 + C2P_kST0_lkST0(iwet) * dkST0_lkST0;

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
        end

        %PLip_PCutoff
        % PLip_PCutoff derivatives
        if (par.opt_PLip_PCutoff)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPLip_PCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			%d2C2P_dPCutoff2 = M3d*0;
			%d2C2P_dPCutoff2(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
            %C2Pxx(:,kk) = d2C2P_dPCutoff2(iwet);
			C2P_PCutoff_lPCutoff = M3d*0;
			C2P_PCutoff_lPCutoff(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
			d2PCutoff_lPCutoff = par.BIO.PLip_PCutoff;
            C2Pxx(:,kk) = C2P_PCutoff_lPCutoff(iwet) * dPCutoff_lPCutoff + C2Px(:,par.pindx.lPLip_PCutoff) * d2PCutoff_lPCutoff;

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
        end

        %PLip_scale
        if (par.opt_PLip_scale)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPLip_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			C2P_PLipscale_lPLipscale = M3d*0;
			C2P_PLipscale_lPLipscale(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
			d2PLipscale_lPLipscale = par.BIO.PLip_scale;
            C2Pxx(:,kk) = C2P_PLipscale_lPLipscale(iwet) * dPLipscale_lPLipscale + C2Px(:,par.pindx.lPLip_scale)*d2PLipscale_lPLipscale;

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
        end

        %PStor_rCutoff
        if (par.opt_PStor_rCutoff)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPStor_rCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			C2P_rCutoff_lrCutoff = M3d*0;
			C2P_rCutoff_lrCutoff(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
			d2rCutoff_lrCutoff = par.BIO.PStor_rCutoff;
            C2Pxx(:,kk) = C2P_rCutoff_lrCutoff(iwet) * drCutoff_lrCutoff + C2Px(:,par.pindx.lPStor_rCutoff)*d2rCutoff_lrCutoff;

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
        end

        %PStor_scale
        if (par.opt_PStor_scale)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lPStor_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			C2P_PStorscale_lPStorscale = M3d*0;
			C2P_PStorscale_lPStorscale(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
			d2PStorscale_lPStorscale = par.BIO.PStor_scale;
            C2Pxx(:,kk) = C2P_PStorscale_lPStorscale(iwet)*dPStorscale_lPStorscale + C2Px(:,par.pindx.lPStor_scale)*d2PStorscale_lPStorscale;

			if (par.opt_alphaS)
                kk = kk + 1;
				C2P_alphaS_lPStorscale = M3d*0;
				C2P_alphaS_lPStorscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = C2P_alphaS_lPStorscale(iwet) * dalphaS_lalphaS;
			end
        end

        % alphaS
        if (par.opt_alphaS)
			kk = kk + 1;

			xim = zeros(size(x));
			xim(par.pindx.lalphaS) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			C2P_alphaS_lalphaS = M3d*0;
			C2P_alphaS_lalphaS(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
			d2alphaS_lalphaS = par.BIO.alphaS;
            C2Pxx(:,kk) = C2P_alphaS_lalphaS(iwet)*dalphaS_lalphaS + C2Px(:,par.pindx.lalphaS)*d2alphaS_lalphaS;
        end

    end %if par.optim
end % if TestVer == 2

end % end function
