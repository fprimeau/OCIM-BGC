 function [out,M] = CellCHNOP(par,x,P,N,T,Irr)
  % DESCRIPTION:
  % This function calculates the C:P ratio of the maximally growing plankton
  % type in the environment specified by the function arguments
	% Currently: only relies on observed PO4,NO3,Temp,&PAR fields

  % CellCNOP adds respiration quotient to previous code CellCNP
  % this function additionally calculates cellular quotas of oxygen and hydrogen

  % INPUTS:
  % Irr is light level (PAR, not total irradiance) in microeinsteins m^-2 s^-1 (i.e. umol photons m^-2 s^-1)
  % T is temperature in Celsius
  % P is phosphate concentration in mol/L (note: OCIM-BGC in mmol/m^3. function input was converted to mol/L in eqC2Puptake)
  % N is nitrate concentration in mol/L (same note for OCIM-BGC)

  % Original Author: George Hagstrom.
  % Modified by: Megan Sullivan

  %%
    on = true; off = false;
	% unpack the parameters to be optimized
	M = par.BIO;

	%Q10Photo
	if (par.opt_Q10Photo == on)
	    lQ10Photo  = x(par.pindx.lQ10Photo);
		M.Q10Photo  = exp(lQ10Photo);
    end

	%fStorage
	if (par.opt_fStorage == on)
		lfStorage = x(par.pindx.lfStorage);
		M.fStorage = exp(lfStorage);
    end

    %fRibE
	if (par.opt_fRibE == on)
		tfRibE = x(par.pindx.tfRibE);
		M.fRibE = 0.5*(1+tanh(tfRibE));
	end

	%kST0
	if (par.opt_kST0 == on)
		lkST0 = x(par.pindx.lkST0);
		M.kST0 = exp(lkST0);
    end

    %PLip_PCutoff
    if (par.opt_PLip_PCutoff == on)
		lPLip_PCutoff = x(par.pindx.lPLip_PCutoff);
		M.PLip_PCutoff = exp(lPLip_PCutoff);
    end
	%PLip_scale
    if (par.opt_PLip_scale == on)
		lPLip_scale = x(par.pindx.lPLip_scale);
		M.PLip_scale = exp(lPLip_scale);
    end

	%PStor_rCutoff
	if (par.opt_PStor_rCutoff == on)
		lPStor_rCutoff = x(par.pindx.lPStor_rCutoff);
		M.PStor_rCutoff = exp(lPStor_rCutoff);
    end
	%PStor_scale
	if (par.opt_PStor_scale == on)
		lPStor_scale = x(par.pindx.lPStor_scale);
		M.PStor_scale = exp(lPStor_scale);
    end
    %alphaS
    if (par.opt_alphaS == on)
		lalphaS = x(par.pindx.lalphaS);
		M.alphaS = exp(lalphaS);
    end
    %gammaDNA
    if (par.opt_gammaDNA == on)
		tgammaDNA = x(par.pindx.tgammaDNA);
        M.gammaDNA = 0.5*(1+tanh(tgammaDNA)); 
    end

  % Read in parameter values from par
     Q10Photo 		= M.Q10Photo;
     fRibE 			= M.fRibE;
     fStorage 		= M.fStorage;
	 kST0 			= M.kST0;
	 alphaS 		= M.alphaS;
     gammaDNA 		= M.gammaDNA;
     gammaLipid 	= M.gammaLipid;
     PLip_PCutoff 	= M.PLip_PCutoff;
     PLip_scale 	= M.PLip_scale;
     PStor_rCutoff 	= M.PStor_rCutoff;
	 PStor_scale 	= M.PStor_scale;
     DNT0 			= M.DNT0;
     DPT0 			= M.DPT0;
     Q10Diffusivity = M.Q10Diffusivity;
	 Q10Bio 		= M.Q10Bio; 
     AMin 			= M.AMin;
     PhiS 			= M.PhiS;
	 alphaPLip 		= M.alphaPLip;
	% %%% everything below here should remain fixed %%%
     pDry 			= M.pDry;
     rho 			= M.rho;
     fProtM 		= M.fProtM;
     fProtL 		= M.fProtL;
     % Define the P, N, and C content of each cellular component (protein, lipid, ribosome, etc.)
     PDNA 			= M.PDNA;
     PRib 			= M.PRib;
     PPhospholipid 	= M.PPhospholipid;
     NProt 			= M.NProt;
	 NRib 			= M.NRib;
     NDNA 			= M.NDNA;
     CProt 			= M.CProt;
     CDNA 			= M.CDNA;
     CPhospholipid 	= M.CPhospholipid;
     CLipid 		= M.CLipid;
	 CRib 			= M.CRib;
    %
    % ----- added hydrogen and oxygen fractions for RQ code update ----
    HProt = 0.0704;
    OProt = 0.23;
    SulfProt = 0.006;

    HLipid = 0.118; 
    OLipid = 0.126;

    HDNA = 0.042;
    ODNA = 0.216; 

    % overwrite previous mass fractions for ribosomes
    HRib = 0.056;
    ORib = 0.304;
    CRib = 0.434; % updated value from 0.419
    NRib = 0.159; % updated value from 0.16
    PRib = 0.046; % updated value from 0.047
    SulfRib = 0.003; 

    % ------

    %
	par.BIO = M;

	nanindx = find(isnan(Irr) | isnan(T) | isnan(N) |isnan(P));
	if ~isempty(nanindx)
		fprintf('Warning: [ %i ] NaNs found in Cell Model input fields. \n',length(nanindx))
    end

%% Define constants
	    molarC = 12.0;              % molar mass of carbon [g/mol]
	    molarN = 14.0;              % molar mass of nitrogen [g/mol]
	    molarP = 31.0;              % molar mass of phosphorus [g/mol]
        molarO = 16.0;              % molar mass of oxygen [g/mol]
        molarS = 32.1;              % molar mass of sulfur [g/mol]
        molarH = 1.01;              % molar mass of hydrogen [g/mol]
		T0 	   = 25.0;              % standard temperature [degC]


	gammaS = gammaDNA+gammaLipid;   % mass fraction of non-membrane or periplasm structural components (gammaDNA+gammaLipid+gammaCarb)
    dgammaS_dgammaDNA = 1;

    rFullA = alphaS/(2.0*AMin);
						 % rFullA is the radius at which if the entire periplasm was filled with
						 % uptake proteins, the percent biomass would be AMin. For cells above this
						 % size, I treat the periplasmic space as being 100% filled with P-uptake
						 % proteins. When AMin = 1, it basically forces all cells to totally invest
						 % in periplasmic proteins, so it disables the part of the model where a
						 % cell can have less uptake proteins to reduce N-quota when P is not too
						 % low.
    

% temperature dependent rates
    DN = DNT0*Q10Diffusivity.^((T-T0)./10.0);    % diffusivity of nitrate [ m^2/hr]
    DP = DPT0*Q10Diffusivity.^((T-T0)./10.0);     % diffusivity of phosphate [ m^2/hr]

    kST = kST0*Q10Bio.^((T-T0)./10.0);  % specific synthesis rate of biosynthetic apparatus [1/hr]
	dkST_dkST0 = Q10Bio.^((T-T0)./10.0);

    alphaI = alphaPhoto(Irr).*Q10Photo.^((T-T0)./10.0); 
				% alphaI is the amount of carbon generated by 1 unit of
				% photosynthetic apparatus. We multiply by Q10Photo to 
                % adjust for the temperature dependence


	CI  = 1+(kST.*(1+PhiS))./alphaI;
                % CI is a constant defined for convenience that is used in finding optimal solutions.
				% It comes from balancing biosynthetic and photosynthetic rates, ie if we set:
				% kST*E = alphaI*L/(1+PhiS), then L = kST*E*(1+PhiS)/alphaI
				% Thus whenever E+L occurs in the equations, we can replace it with E*CI

    % derivatives of CI
	dCI_dQ10Photo = (kST.*(1+PhiS))./alphaPhoto(Irr) .* (-(T-T0)./10.0) .*Q10Photo.^((-(T-T0)./10.0)-1) ;
	dCI_dkST = (1+PhiS)./alphaI;

% | cell membrane is made up of lipid and protein
  %  fProtM = 0.25;               % protein fraction of cell membranes
    fLipidM = 1-fProtM;           % Lipid fraction of cell membranes 


% | light harvesting apparatus is 70% protein. The rest is Lipid
  %  fProtL = .7;            % protein fraction of light harvesting apparatus
	fLipidL = 1-fProtL; 	 % lipid fraction of light harvesting apparatus

% | Biosynthetic apparatus is ribosomes and proteins
	fProtE = 1-fRibE;



%% FUNCTIONAL POOLS
    %%% Calculate mass fractions of each element in each functional pool by
    %%% finding average of elemental mass fractions in all molecules within 
    %%% the pool, weighted by weight fraction of each molecule in the pool
    %
    % Composition of each pool:
    %%% _L = light/photosynthesis --> (proteins) % +lipids??
    %%% _M = cell membrane       --> (lipids + proteins)
    %%% _E = biosynthesis         --> (proteins,ribosomes,)
    %%% _S = other structure components (not membrane) --> (DNA +lipids)
    %%% Ctotal = CL*L + CM*M + CE*E + CS*S_other
    %%% 1 = L + M + E + S
    % periplasmic space (A) - prokaryotes have periplasmic space, eukaryotes dont - generally filled with proteins for nutrient uptake. can be a big fraction of mass for small bateria (maybe 40% or more)
      % fProtA - cells only invest in periplasmic proteins if they need too.
      % exact formula for how periplasm relates to p uptake is too complicated

    %%% Photosynthetic pool CNP
    NL = fProtL*NProt;        % N fraction of light harvesting apparatus (assumes all N is in proteins)
                              %   = fraction of light apparatus that is protein * fraction of protein that is N
    CL = CProt*fProtL+CLipid*(1-fProtL);     % Carbon fraction of light harvesting apparatus
    %
    HL    = HProt*fProtL+HLipid*(1-fProtL);     % Hydrogen weight fraction
    OL  = OProt*fProtL+OLipid*(1-fProtL);   % Oxygen weight fraction
    SulfL = SulfProt*fProtL;                 % Sulfur weight fraction (only in proteins)

    %%% Membrane pool CNP (membrane carbon from proteins + lipids) 
    NM = fProtM*NProt;        % N fraction of membrane
                              %   = fraction of cell membrane that is protein * fraction of protein that is N
    CM = fProtM*CProt+fLipidM*CLipid;     % carbon fraction of membrane 
	%PM     % for M&L,only P content is from phospholipids, which are added at the end with the luxury pool
    %
    HM    = fProtM*HProt + fLipidM*HLipid;     % Hydrogen weight fraction
    OM    = fProtM*OProt + fLipidM*OLipid;   % Oxygen weight fraction
    SulfM = fProtM*SulfProt;                 % Sulfur weight fraction (only in proteins)

    %%% Biosynthetic pool CNP
    NE = NProt*fProtE+NRib*fRibE;  % Nitrogen fraction of biosynthetic components
    CE = CProt*fProtE+CRib*fRibE;  % Carbon fraction of biosynthesis pool
    PE = fRibE*PRib;               % PRib is different for prokaryotes and eukaryotes. the input value is for prokaryotes
    %
    HE    = HProt*fProtE + HRib*fRibE;     % Hydrogen weight fraction
    OE    = OProt*fProtE + ORib*fRibE;   % Oxygen weight fraction
    SulfE = SulfProt*fProtE + SulfRib*fRibE;                 % Sulfur weight fraction (only in proteins)
    
    %%% Other structure pool
    NS = gammaDNA/gammaS*NDNA;                          % N fraction of other structural components
    CS = (gammaDNA*CDNA+gammaLipid*CLipid)/gammaS;      % Carbon fraction of other structural components
                                % fLipid stands for other lipids that are not phospholipids.
                                % cell has some energy reserves always as part of its structural pool
	PS = gammaDNA/gammaS*PDNA;
    %
    HS = (gammaDNA*HDNA+gammaLipid*HLipid)/gammaS;     % Hydrogen weight fraction
    OS = (gammaDNA*ODNA+gammaLipid*OLipid)/gammaS;  % Oxygen weight fraction
    SulfS = 0;                  % Sulfur weight fraction (only in proteins)

    %%% derivatives of pools wrt optimizable parameters
    dNE_dfRibE = -NProt+NRib;
    dCE_dfRibE = -CProt+CRib;  
    dPE_dfRibE = PRib;

    dPS_dgammaDNA = (gammaS*PDNA-(gammaDNA*PDNA))/(gammaS^2);
    dCS_dgammaDNA = gammaLipid*(CDNA-CLipid)/gammaS^2; 
    dNS_dgammaDNA = (gammaS*NDNA-(gammaDNA*NDNA))/(gammaS^2);
    %
    dHS_dgammaDNA = gammaLipid*(HDNA-HLipid)/gammaS^2; 
    dOS_dgammaDNA = gammaLipid*(ODNA-OLipid)/gammaS^2; 

%% Solving for PLim
    aP= (3*DP.*P*1000)*molarP/10^6 /(pDry*rho);
	% units = P:[mol/L]*[1000L/m^3]*DP:[m^2/hr]-> [mol/m/hr]*molarP:[g/mol] -> [g/m/hr]/[10^6um/m] -> [g/um/hr] % rho:[g/um^3] % alphaS:[um]
    % units of aP = [um^2/hr]*[gP/gDrymass]
    coefA = aP.*CI.^2/alphaS^2-kST*PE;	
    coefB = -2*(1-gammaS)*CI.*aP/alphaS^2-gammaDNA*PDNA*kST;
    coefC = (1-gammaS)^2.*aP/alphaS^2;

    EPLim = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);

    rPLim = alphaS./(1-gammaS - CI.*EPLim); %1.0./( (1-gammaS)/alphaS-CI.*EPLim/alphaS);
    APLim = (1-gammaS -CI.*EPLim)./2; %should equal AMin;  APLim = alphaS/(2*rPLim)
    fProtAPLim = APLim.*0+1; % =1 (if purely P limited = invest as much as you can in uptake proteins)


%%% derivatives of EPlim w.r.t. parameters
	%%%% dE_dDIP
	daP_dDIP = 3*DP*molarP/10^6 /(pDry*rho)*1000;
	dcoefA_daP = CI.^2/alphaS^2;
	dcoefB_daP = -2*(1-gammaS)*CI/alphaS^2;
	dcoefC_daP = (1-gammaS)^2/alphaS^2;

	dEPLim_daP = ((2*coefA.*(-dcoefB_daP - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_daP-4*dcoefA_daP.*coefC-4*coefA.*dcoefC_daP)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_daP)./(4*coefA.^2) );
	dEPLim_dDIP = dEPLim_daP.*daP_dDIP;

	%%%% dE_dCI
    dcoefA_dCI = 2.*aP.*CI/alphaS^2;
	dcoefB_dCI = -2*(1-gammaS).*aP/alphaS^2;

	dEPLim_dCI = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );

	%%%% dE_dkST
	dcoefA_dkST = -PE +dcoefA_dCI.*dCI_dkST;
	dcoefB_dkST = -gammaDNA*PDNA + dcoefB_dCI.*dCI_dkST;
	dEPLim_dkST = ((2*coefA.*(-dcoefB_dkST - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dkST-4*dcoefA_dkST.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dkST)./(4*coefA.^2) );

	%%%% dE_dalphaS
	dcoefA_dalphaS = -2*aP.*CI.^2/alphaS^3;
	dcoefB_dalphaS = -2*(-2*(1-gammaS)*CI.*aP)/alphaS^3;
	dcoefC_dalphaS = -2*(1-gammaS)^2.*aP/alphaS^3;

	dEPLim_dalphaS = ((2*coefA.*(-dcoefB_dalphaS - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dalphaS-4*dcoefA_dalphaS.*coefC-4*coefA.*dcoefC_dalphaS)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dalphaS)./(4*coefA.^2) );

    %%%% dE_dCI
    dcoefA_dgammaDNA = 0;
	dcoefB_dgammaDNA = 2*CI.*aP./alphaS^2 -PDNA*kST;
    dcoefC_dgammaDNA = -2*(1-gammaS).*aP./alphaS^2;

    dEPLim_dgammaDNA = ((2*coefA.*(-dcoefB_dgammaDNA - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dgammaDNA-4*dcoefA_dgammaDNA.*coefC-4*coefA.*dcoefC_dgammaDNA)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dgammaDNA)./(4*coefA.^2) );


%% Solving for NLim
    aN = 3*DN.*(N*1000)*molarN/10^6 /(pDry*rho); %converting top to g/um/hr
    coefA = CI.^2.*aN/alphaS^2+kST.*CI*NM/2.0-kST*NL.*(CI-1)-kST*NE;
    coefB = -2*aN*(1-gammaS).*CI/alphaS^2-kST*gammaDNA*NDNA-kST*NM*(1-gammaS)/2.0-kST*AMin*NProt;
    coefC = aN.*(1-gammaS)^2/alphaS^2;

	ENLim = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);

    rNLim = alphaS./(1-gammaS - CI.*ENLim);
	ANLim = (1-gammaS-CI.*ENLim)./2;    % mass fraction of periplasm in the cell. (= alphaS./(2*rNLim))
	fProtANLim = AMin./ANLim;           % protein fraction of periplasm

%%% derivatives of ENlim w.r.t. parameters
	%%%% dE_dDIP
	daN_dDIP = 0;
	dENLim_dDIP = aN.*0;

	%%%% dE_dCI
    dcoefA_dCI = 2*CI.*aN/alphaS^2 + kST*NM/2.0 -kST*NL;
    dcoefB_dCI = -2*aN*(1-gammaS)/alphaS^2;

    dENLim_dCI = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );

	%%%% dE_dkST   % check this!
	dcoefA_dkST = dcoefA_dCI.*dCI_dkST +CI*NM/2.0 -NL.*(CI-1) -NE;
	dcoefB_dkST = dcoefB_dCI.*dCI_dkST -gammaDNA*NDNA-NM*(1-gammaS)/2.0-AMin*NProt;
	dENLim_dkST = ((2*coefA.*(-dcoefB_dkST - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dkST-4*dcoefA_dkST.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dkST)./(4*coefA.^2) );

	%%%% dE_alphaS
	dcoefA_dalphaS = -2*CI.^2.*aN/alphaS^3;
	dcoefB_dalphaS = -2*(-2*aN*(1-gammaS).*CI)/alphaS^3;
	dcoefC_dalphaS = -2*aN.*(1-gammaS)^2/alphaS^3;

	dENLim_dalphaS = ((2*coefA.*(-dcoefB_dalphaS - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dalphaS-4*dcoefA_dalphaS.*coefC-4*coefA.*dcoefC_dalphaS)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dalphaS)./(4*coefA.^2) );

	%%%% dE_dfRibE
	dcoefA_dNE = -kST;
	dcoefB_dNE = 0;
	dENLim_dNE = ((2*coefA.*(-dcoefB_dNE - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dNE-4*dcoefA_dNE.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dNE)./(4*coefA.^2) );
	dENLim_dfRibE = dENLim_dNE.*dNE_dfRibE;

    %%%% dE_dgammaDNA
	dcoefA_dgammaDNA = 0;
    dcoefB_dgammaDNA = 2*aN.*CI./alphaS^2 - kST.*NDNA +kST.*NM./2;
    dcoefC_dgammaDNA = -2*aN.*(1-gammaS)./alphaS^2;
    dENLim_dgammaDNA = ((2*coefA.*(-dcoefB_dgammaDNA - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dgammaDNA-4*dcoefA_dgammaDNA.*coefC-4*coefA.*dcoefC_dgammaDNA)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dgammaDNA)./(4*coefA.^2) );


%% reset ANLim to AMin if the mass fraction of periplasm is too small (for large phytoplankton),
	rLG_ind = find(rNLim>rFullA);  % equivalent to find(ANLim<AMin) ;
		fProtANLim(rLG_ind) = 1.0; % fill periplasmic space with uptake proteins

        coefA = CI.^2.*aN/alphaS^2+kST.*CI*(NM+NProt)/2.0-kST*NL.*(CI-1)-kST*NE;
        coefB = -2*aN*(1-gammaS).*CI/alphaS^2-kST*gammaDNA*NDNA-kST*(NM+NProt)*(1-gammaS)/2.0;
        coefC = aN.*(1-gammaS)^2/alphaS^2;

        ENLim_LG = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);
        rNLim_LG = alphaS./(1-gammaS - CI.*ENLim_LG);
        ANLim_LG = (1-gammaS-CI.*ENLim_LG)./2; % mass fraction of periplasm in the cell. (= alphaS./(2*rNLim))

		ENLim(rLG_ind) = ENLim_LG(rLG_ind);
		rNLim(rLG_ind) = rNLim_LG(rLG_ind);
        ANLim(rLG_ind) = ANLim_LG(rLG_ind);

    %%% derivatives of ENlim for large phyto w.r.t. parameters
		%%%% dENLim_dDIP = aN.*0; %stays zero in both NLim cases
		%%%% dE_dCI
        dcoefA_dCI = 2*CI.*aN/alphaS^2 + kST*(NM+NProt)/2.0 -kST*NL;
        dcoefB_dCI = -2*aN*(1-gammaS)/alphaS^2;
        dENLim_dCI_LG = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );
        dENLim_dCI(rLG_ind) = dENLim_dCI_LG(rLG_ind);

		%%%% dE_dkST
		dcoefA_dkST = dcoefA_dCI.*dCI_dkST +CI*(NM+NProt)/2.0-NL.*(CI-1)-NE;
		dcoefB_dkST = dcoefB_dCI.*dCI_dkST -gammaDNA*NDNA-(NM+NProt)*(1-gammaS)/2.0;
		dENLim_dkST_LG = ((2*coefA.*(-dcoefB_dkST - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dkST-4*dcoefA_dkST.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dkST)./(4*coefA.^2) );
		dENLim_dkST(rLG_ind) =  dENLim_dkST_LG(rLG_ind);

		% dE_dalphaS doesnt change for large case, so this bbit of code is redundant
		dcoefA_dalphaS = -2*CI.^2.*aN/alphaS^3;
		dcoefB_dalphaS = -2*(-2*aN*(1-gammaS).*CI)/alphaS^3;
		dcoefC_dalphaS = -2*aN.*(1-gammaS)^2/alphaS^3;
		dENLim_dalphaS_LG = ((2*coefA.*(-dcoefB_dalphaS - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dalphaS-4*dcoefA_dalphaS.*coefC-4*coefA.*dcoefC_dalphaS)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dalphaS)./(4*coefA.^2) );
		dENLim_dalphaS(rLG_ind) = dENLim_dalphaS_LG(rLG_ind);

		% dfRibE
		dcoefA_dNE = -kST;
		dcoefB_dNE = 0;
		dENLim_dNE = ((2*coefA.*(-dcoefB_dNE - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dNE-4*dcoefA_dNE.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dNE)./(4*coefA.^2) );
		dENLim_dfRibE_LG = dENLim_dNE.*dNE_dfRibE;
		dENLim_dfRibE(rLG_ind) = dENLim_dfRibE_LG(rLG_ind);

        %%%% dE_dgammaDNA
	    dcoefA_dgammaDNA = 0;
        dcoefB_dgammaDNA = 2*aN.*CI./alphaS^2 - kST.*NDNA +kST.*(NM+NProt)./2;
        dcoefC_dgammaDNA = -2*aN.*(1-gammaS)./alphaS^2;
        dENLim_dgammaDNA_LG = ((2*coefA.*(-dcoefB_dgammaDNA - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dgammaDNA-4*dcoefA_dgammaDNA.*coefC-4*coefA.*dcoefC_dgammaDNA)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dgammaDNA)./(4*coefA.^2) );

        dENLim_dgammaDNA(rLG_ind) = dENLim_dgammaDNA_LG(rLG_ind);

%% Calculate growth rates (biosynthetic rate) under each condition 
    muNLim = kST.*ENLim;  % units: [1/hr]
    muPLim = kST.*EPLim;  % units: [1/hr]

    muNLim_P= aP.*fProtANLim./(rNLim.^2.*(ENLim*PE+gammaDNA*PDNA));
    muPLim_N= aN./(rPLim.^2.*((CI-1)*NL.*EPLim+EPLim*NE+gammaDNA*NDNA+alphaS./(2.0.*rPLim)*NProt+alphaS./(2.0.*rPLim)*NM));

        % The optimal strategy will either be N limited, P-limited, or
        % co-limited. The method for determining which starts by looking at how
        % much P the cell gets when it optimizes for N limitation, and how much N
        % it gets when it optimizes for P limitation.
        % 
		% muNLim_P is the "P-limited" growth rate at the strategy where
		% N-limitation, Biosynthesis, and Photosynthesis balance.
		% (P-limited growth rate equation, evaluated at the N limited solution for E)
		% Likewise, muPLim_N is the "N-Limited" growth rate at the strategy where
		% P-limitation, biosynthesis, and photosynthesis balance.
		% If, for example, muNLim_P < kST*ENLim, then we conclude that the NLim
		% strategy is suboptimal (thus, not N-limited), and vice versa for the P-lim strategy.

%% initialize vector fields
all_nan = NaN(size(P));

fProtAOpt = all_nan;
EOpt = all_nan;
LimType = all_nan;

dfProtAOpt_dCI = all_nan;
dfProtAOpt_dE = all_nan;
dE_dCI = all_nan;
dE_dalphaS = all_nan;
dE_dfRibE = all_nan;
pdfProtAOpt_fRibE =all_nan;
dE_dkST = all_nan;
dfProtAOpt_dkST =all_nan;
dE_dDIP = all_nan;
dfProtAOpt_dDIP = all_nan;
dE_dgammaDNA = all_nan;
dfProtAOpt_dgammaDNA =  all_nan;

errcount = 0;
%% loop through points to calculate cell quotas
% handle definitions needed to run loop in parallel
EColim_handle = @EColim;
AColim_handle = @AColim;
%parfor i =1:length(P)   % to run in parallel
for i =1:length(P)
    % for debugging: these NaN values should never appear in final vectors
	EOpt_i = NaN;
	fProtAOpti = NaN;
	LimState = -1;
    if muNLim(i)<muPLim(i)
        if muNLim(i)<muNLim_P(i) % N Limitation
            LimState = 0;
            EOpt_i = ENLim(i);
            fProtAOpti = fProtANLim(i); % = 2*AMin./(1-gammaS-CI.*ENLim)

            %%% derivatives
            dE_dCI(i) = dENLim_dCI(i);
			dE_dalphaS(i) = dENLim_dalphaS(i);
			dE_dfRibE(i) = dENLim_dfRibE(i);
			dE_dkST(i) = dENLim_dkST(i);
			dE_dDIP(i) = 0; 
            dE_dgammaDNA(i) = dENLim_dgammaDNA(i);

			if ANLim(i)>AMin
				dfProtAOpt_dE(i) = fProtAOpti^2*CI(i)/(2*AMin);
				dfProtAOpt_dCI(i) = fProtAOpti^2/(2*AMin)*(EOpt_i);
                pdfProtAOpt_fRibE(i) = 0;
				dfProtAOpt_dkST(i) = dfProtAOpt_dE(i)*dE_dkST(i)+dfProtAOpt_dCI(i)*dCI_dkST(i);
				dfProtAOpt_dDIP(i) = dfProtAOpt_dE(i)*dE_dDIP(i);
                dfProtAOpt_dgammaDNA(i) = dfProtAOpt_dE(i)*dE_dgammaDNA(i) + 2*AMin./(1-gammaS-CI(i)*EOpt_i)^2;
			else
				dfProtAOpt_dE(i) = 0;
				dfProtAOpt_dCI(i) = 0;
                pdfProtAOpt_fRibE(i) = 0;
				dfProtAOpt_dkST(i) = 0;
				dfProtAOpt_dDIP(i) = 0;
                dfProtAOpt_dgammaDNA(i) = 0;
            end
        else
            LimState = 2; % Co-limitation: weak N limitation

            EMin=1e-4;
            EMax=(1-gammaS)/CI(i);
			[EOpt_i, dE_dCI(i), dE_dalphaS(i), dE_dfRibE(i),pdE_kST,dE_daP,dE_dgammaDNA(i)] = feval(EColim_handle,EMin,EMax,i); %if using parfor loop
			dE_dkST(i) = pdE_kST + dE_dCI(i).*dCI_dkST(i);
			dE_dDIP(i) = dE_daP*daP_dDIP(i);

			[fProtAOpti, dfProtAOpt_dCI(i), dfProtAOpt_dE(i), pdfProtAOpt_fRibE(i),pdfProtAOpt_DIP,pdfProtAOpt_gammaS,dfProtAOpt_dNS,dfProtAOpt_dPS] = feval(AColim_handle,EOpt_i,CI(i),i);
			dfProtAOpt_dkST(i) = dfProtAOpt_dE(i)*dE_dkST(i)+dfProtAOpt_dCI(i)*dCI_dkST(i);
			dfProtAOpt_dDIP(i) = pdfProtAOpt_DIP + dfProtAOpt_dE(i)*dE_dDIP(i);
            dfProtAOpt_dgammaDNA(i) = dfProtAOpt_dE(i)*dE_dgammaDNA(i) + pdfProtAOpt_gammaS*1 + dfProtAOpt_dNS*dNS_dgammaDNA + dfProtAOpt_dPS*dPS_dgammaDNA;

        end
    elseif muPLim(i) <= muNLim(i)
        if muPLim(i)<muPLim_N(i) % P Limitation
            LimState = 1;
            EOpt_i = EPLim(i);
			fProtAOpti = fProtAPLim(i); % = 1

			%%% derivatives
            dE_dCI(i) = dEPLim_dCI(i);
			dE_dalphaS(i) = dEPLim_dalphaS(i);
            dfProtAOpt_dE(i) = 0;
			dfProtAOpt_dCI(i) = 0;
			dE_dfRibE(i) = 0;
            pdfProtAOpt_fRibE(i) = 0;
			dE_dkST(i) = dEPLim_dkST(i);
			dfProtAOpt_dkST(i) = 0;
			dE_dDIP(i) = dEPLim_dDIP(i);
			dfProtAOpt_dDIP(i) = 0;
            dE_dgammaDNA(i) = dEPLim_dgammaDNA(i);
            dfProtAOpt_dgammaDNA(i) = 0;

        else
            LimState = 3; %Colimitation: Weak P Limitation. value changed from 2 for debugging
            EMin = 1e-4;
            EMax = (1-gammaS)/CI(i);
			[EOpt_i, dE_dCI(i), dE_dalphaS(i), dE_dfRibE(i),pdE_kST,dE_daP,dE_dgammaDNA(i)] = feval(EColim_handle,EMin,EMax,i);
			dE_dkST(i) = pdE_kST + dE_dCI(i).*dCI_dkST(i);
			dE_dDIP(i) = dE_daP*daP_dDIP(i);

            % solve for AColim(E), must be between 0 and 1
			[fProtAOpti, dfProtAOpt_dCI(i), dfProtAOpt_dE(i), pdfProtAOpt_fRibE(i),pdfProtAOpt_DIP,pdfProtAOpt_gammaS,dfProtAOpt_dNS,dfProtAOpt_dPS] = feval(AColim_handle,EOpt_i,CI(i),i);
			dfProtAOpt_dkST(i) = dfProtAOpt_dE(i)*dE_dkST(i)+dfProtAOpt_dCI(i)*dCI_dkST(i);
			dfProtAOpt_dDIP(i) = pdfProtAOpt_DIP + dfProtAOpt_dE(i)*dE_dDIP(i);
            dfProtAOpt_dgammaDNA(i) = dfProtAOpt_dE(i)*dE_dgammaDNA(i) + pdfProtAOpt_gammaS*1 + dfProtAOpt_dNS*dNS_dgammaDNA + dfProtAOpt_dPS*dPS_dgammaDNA;
        end
    end % end cases
	EOpt(i) = EOpt_i;
    fProtAOpt(i) = fProtAOpti;
    LimType(i) = LimState;
end % end for loop

if errcount >0
    fprintf('complex_cubic solver error count: %d \n',errcount)
end

AOpt = (1-gammaS-CI.*EOpt)./2;        % optimal periplasm fraction of cell
MOpt = AOpt;                          % optimal membrane fraction of cell
rOpt = alphaS./(1-gammaS - CI.*EOpt); % optimal cell radius
muOpt = kST.*EOpt;                    % optimal growth rate

%% %%%%% Calculate Cell quotas in moles P, C, and N %%%%%%%%%%%%%%%
QPfunctional = (EOpt.*PE + gammaDNA*PDNA )/molarP;

QCfunctional = (EOpt.*CE +EOpt.*(CI-1)*CL +gammaS*CS +fProtAOpt.*AOpt.*CProt +MOpt.*CM)./molarC; 
QCreserve = QCfunctional.*(24.0*.25*kST.*EOpt.*(1+PhiS));
            %%% C reserve pool = extra carbohydrates, stored to allow cells
            %%% to continue growing during the dark hours, modeled as a set
            %%% fraction of the C quota in all functional pools. The
            %%% fraction depends on growth rate. Faster growing cells need
            %%% to store a greater fraction of their functional C quota in
            %%% the storage reserve, to be used to synthesize biomass at
            %%% night
% QC = QCfunctional + QCreserve;
QC = QCfunctional .*(1 + 24.0*.25*kST.*EOpt.*(1+PhiS));
QN = (EOpt.*NE +EOpt.*(CI-1).*NL +gammaS*NS +fProtAOpt.*AOpt.*NProt +MOpt.*NM)./molarN;

%%% ------ Calculate Phosphorus Cell Quota ---------
%%%%%%% Add P storage and phospholipids %%%%%%%%%%%%%%%%%%%
% each scaled by logistic function with 2 parameters each: PLip: PLip_scale & PLip_PCutoff; PStor: PStor_scale & PStor_rCutoff
% can turn of PLip and PStorage by setting alphaPLip = 0 and fStorage = 0;

    % PLip: one of the main ways that small cells can vary their P-quota
    % is by substituting Phospholipids for Sulfolipids at low P.
    % This only really matters for the small cells

    % PStor: model treats storage differently for large cells
	% compared to small cells. In particular, large cells, with 
    % radius>>rCutoff, maximize luxury storage, while storage approaches 
    % zero for cells much smaller than rCutoff.

    %f_plip = alphaPLip./CPhospholipid

 PLip = alphaPLip*MOpt*PPhospholipid./CPhospholipid.*molarC./molarP.*(1./(1 + exp(-PLip_scale*(P-PLip_PCutoff)))); %unit: [molP/molC]

 PStor = fStorage.*5000.*P .*(1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff)))); %unit: [molP/molC]

%%%%%% Recalculate Phosphorus Cell Quota %%%%%%%%%%%%%%%%%
% defined so that CP = 1./(1./(QC./QP)+ PStor + PLip);
PLip = PLip.*QC;		% [molP/molC*molC = molP]
PStor = PStor.*QC;	    % [molP/molC*molC = molP]
%QP =(EOpt.*PE +(gammaDNA)*PDNA)/molarP + PLip +PStor;
QP = QPfunctional + PLip +PStor;
% -------------------------------------------------------------------
%%% ---- Oxygen Hydrogen and Sulfur quotas added for respiration quotient calculation ---- 
%%%%% DEFINE CONSTANTS (Move to begining)
% Assumes phospholipids have the molecular formula: C_37.9 H_72.5 O_9.4 N_0.43 P_1 (Geider & La Roche 2002)
H2P_PLip = 72.5;
O2P_PLip = 9.4; 
N2P_PLip = 0.43;
% Assume all phosphate is stored in inorganic form (H3PO4), 
% so PStor does not affect the oxidation state of organic matter. 
% We can exclude PStor from RQ equation. (the P,H,and O contributions to 
% RQ will cancel. i.e. adds zero in numerator of RQ equation)
O2P_PStor = 4; 
H2P_PStor = 3;


% QCHNOPS_PLip = [ 37.7 72.5 0.43 9.4 1 0];
% QCHNOPS_carb = [6 12 0 6 0 0];
% O2P_PLip = QCHNOPS_PLip(4)/QCHNOPS_PLip(5);

%%% Oxygen Quota
QOxyfunctional = (EOpt.*OE +EOpt.*(CI-1)*OL +gammaS*OS + fProtAOpt.*AOpt.*OProt + MOpt.*OM)/molarO;
QOxyreserve = QCreserve; 
            %%% Oxygen reserve should be the same as the carbon reserve,
            %%% converted to oxygen, based on the molecular formula for
            %%% carbohydrates: C6 H12 O6. In carbohydrates, there are 6
            %%% moles of Oxygen for every 6 moles of Carbon, so conversion
            %%% factor is 1:1. Thus, QOxyreserve = QCreserve
%
QOxyPStor = O2P_PStor*PStor; 
            %%% Assumes all luxury phosphorus storage has the molecular 
            %%% formula of (H3PO4)
QOxyPLip = O2P_PLip*PLip; 
            %%% Assumes phospholipids have the molecular formula: 
            %%% C_37.9 H_72.5 O_9.4 N_0.43 P_1 (Geider & La Roche 2002)
QOxy = QOxyfunctional + QOxyreserve + QOxyPStor + QOxyPLip ;

%%% Hydrogen Quota
QHfunctional = (EOpt.*HE + EOpt.*(CI-1)*HL + gammaS*HS + fProtAOpt.*AOpt.*HProt + MOpt.*HM)./molarH;
QHreserve = 2*QCreserve;
QHPLip = H2P_PLip*PLip;
QHPStor = H2P_PStor*PStor;
QH = QHfunctional + QHreserve + QHPLip + QHPStor;

%%% Sulfur quota (only in proteins)
QSulf = (EOpt.*SulfE + EOpt.*(CI-1)*SulfL + gammaS*SulfS + fProtAOpt.*AOpt.*SulfProt + MOpt.*SulfM)./molarS;


%% %%%% Calculate Stoichiometric Ratios %%%%%%%%%%%%%%%%
CP = QC./QP;
NP = QN./QP;
CN = QC./QN;

%----- Calculate Total Respiration Quotient -----
% Total includes oxygen required to remineralize NH3 to HNO3
%
%           (QC + 0.25*QH - 0.5*QO +1.25*QN + 1.25*QP + 1.5*QS)
% r_O2toC = ---------------------------------------------------
%                               QC

RQtotalO2toC = (QC + 0.25*QH - 0.5*QOxy +1.25*QN + 1.25*QP + 1.5*QSulf)./QC;
% RQ using quotas including Pstorage (above) are equivalent to RQ using
% quotas excluding Pstorage (below)

% remove PStor component of QOxy QH and QP 
RQtotalO2toC = (QC + 0.25*(QHfunctional+QHreserve+QHPLip) - 0.5*(QOxyfunctional + QOxyreserve + QOxyPLip ) +1.25*QN + 1.25*(QPfunctional+PLip) + 1.5*QSulf)./QC;

%%


% ****** Calculate derivatives of C:P ****** %

		% QC(E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kST,PhiS)
		% Take partial derivatives with respect to all these variables
		pdQC_E = (CE + (CI-1)*CL).*(1+24.0*.25*kST.*EOpt.*(1+PhiS))./molarC + (EOpt.*CE +EOpt.*(CI-1)*CL +gammaS*CS +fProtAOpt.*AOpt.*CProt +MOpt.*CM).*24.0*0.25.*kST*(1+PhiS)/molarC;
        pdQC_CE = EOpt./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
		pdQC_CI = (EOpt*CL).*(1+24.0*.25*kST.*EOpt*(1+PhiS))/molarC;
		pdQC_CL = 0; % do later if needed
		pdQC_gammaS = CS./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS)); % do later
		pdQC_CS = gammaS./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
        pdQC_fProtAOpt = AOpt.*CProt./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
        pdQC_AOpt = fProtAOpt.*CProt./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
        pdQC_MOpt =  CM./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
		%pdQC_alphaS = 0;
		%pdQC_rOpt = (-fProtAOpt*alphaS./(2.0*rOpt.^2)*CProt-alphaS./(2.0*rOpt.^2)*CM).*(1+24.0*.25*kST.*EOpt.*(1+PhiS))/molarC;
		pdQC_rOpt = 0;
        pdQC_CProt = NaN; % unused
		pdQC_CM = NaN; % unused
		pdQC_kST = 24*0.25.*EOpt.*(1+PhiS).*(EOpt.*CE +EOpt.*(CI-1)*CL +gammaS*CS +fProtAOpt.*AOpt.*CProt +MOpt.*CM)./molarC ;
		pdQC_PhiS = 24*0.25.*kST.*EOpt.*(EOpt.*CE +EOpt.*(CI-1)*CL +gammaS*CS +fProtAOpt.*AOpt.*CProt +MOpt.*CM)./molarC ;

		% QP(E,fRibE,PRib,gammaDNA,PDNA)
		% Take partial derivatives with respect to all these variables
		% QP depends on QC through the PStor andPLip terms
		pdQP_PLip = 1;
		pdQP_PStor = 1 ;
		pdPLip_QC = alphaPLip*MOpt*PPhospholipid./CPhospholipid.*molarC./molarP.*(1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))));
		pdPStor_QC = fStorage.*5000.*P .*(1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))));
		pdQP_QC = pdPLip_QC + pdPStor_QC;
		pdQP_E = PE./molarP + pdQP_QC.*pdQC_E;
		pdQP_PE = EOpt./molarP + 0; %QC does not depend on PE
		pdQP_gammaDNA = PDNA/molarP;	% QC depends on gammaS, which is defined from gammaDNA; accounted for later in pdQP_gammaS; this term is never used in the CP derivatives, so may be wrong; dgammaS_dgammaDNA = 1
		pdQP_PDNA = gammaDNA/molarP + 0; %QC does not depend on PDNA
		% PLip terms
		pdQP_MOpt = pdQP_PLip.*(PLip./MOpt) + pdQP_QC.*pdQC_MOpt;
		pdQP_alphaPLip = pdQP_PLip.*(PLip./alphaPLip);	% no alphaPLip in QC
		pdQP_PLipscale = (pdQP_PLip).*QC.*alphaPLip.*MOpt.*PPhospholipid./CPhospholipid.*molarC./molarP .* (-1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))).^2) .*exp(-PLip_scale*(P-PLip_PCutoff)) .* -1.*(P-PLip_PCutoff);
		pdQP_PLipPCutoff = (pdQP_PLip).*QC.*alphaPLip.*MOpt.*PPhospholipid./CPhospholipid.*molarC./molarP .* (-1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))).^2) .*exp(-PLip_scale*(P-PLip_PCutoff)) .*PLip_scale;
		%PStor terms
		pdQP_fStorage = pdQP_PStor.*(PStor./fStorage); 	% QC does not depend on fStorage
		pdQP_rOpt = pdQP_PStor.*QC.* fStorage.*5000.*P .*(-1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))).^2) .* exp(-PStor_scale.*(rOpt-PStor_rCutoff)) .*-PStor_scale + pdQP_QC.*pdQC_rOpt;	% QC does not depend directly on rOpt
		pdQP_rCutoff = (pdQP_PStor).*QC.*fStorage*5000.*P .*(-1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))).^2) .* exp(-PStor_scale.*(rOpt-PStor_rCutoff)) .*PStor_scale;
        pdQP_PStorscale = (pdQP_PStor) * -1*PStor.^2./(QC.*fStorage*5000.*P) .*exp(-PStor_scale.*(rOpt-PStor_rCutoff)) .* -1.*(rOpt-PStor_rCutoff);
		%DIP is in both PLip and PStor; easier to follow if seperate then add together
		pdPStor_DIP = QC.*fStorage.*5000.*(1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))));
		pdPLip_DIP = QC.*alphaPLip.*MOpt*PPhospholipid/CPhospholipid*molarC/molarP.*(-1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))).^2).*(-PLip_scale*exp(-PLip_scale*(P-PLip_PCutoff))) ;
		pdQP_DIP = pdQP_PLip.*pdPLip_DIP + pdQP_PStor.*pdPStor_DIP;
		% QC only terms
		pdQP_CE = pdQP_QC.*pdQC_CE;
		pdQP_CI = pdQP_QC.*pdQC_CI;
		pdQP_CL = pdQP_QC.*pdQC_CL; % do later
		pdQP_gammaS = pdQP_QC.*pdQC_gammaS;
		pdQP_CS = pdQP_QC.*pdQC_CS;
		pdQP_fProtAOpt = pdQP_QC.*pdQC_fProtAOpt;
		pdQP_AOpt = pdQP_QC.*pdQC_AOpt;
		pdQP_CM = pdQP_QC.*pdQC_CM;		% do later
		pdQP_CProt = pdQP_QC.*pdQC_CProt; % do later
		pdQP_kST = pdQP_QC.*pdQC_kST;
		pdQP_PhiS = pdQP_QC.*pdQC_PhiS;


		% C:P = C:P(E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kST,PhiS,fRibE,PRib,gammaDNA,PDNA)
		% Use the partial derivatives of QC and QP to calculate the partial derivatives of C:P. pdC2P_E is the only complicated ones
		pdC2P_E = (QP.*pdQC_E - QC.*pdQP_E)./QP.^2;
		pdC2P_CE = (QP.*pdQC_CE - QC.*pdQP_CE)./QP.^2; % pdQC_CE./QP;
		pdC2P_CI = (QP.*pdQC_CI - QC.*pdQP_CI)./QP.^2; % pdQC_CI./QP;
		pdC2P_CL = (QP.*pdQC_CL - QC.*pdQP_CL)./QP.^2; % pdQC_CL./QP;
		pdC2P_gammaS = (QP.*pdQC_gammaS - QC.*pdQP_gammaS)./QP.^2; % pdQC_gammaS./QP;
		pdC2P_CS = (QP.*pdQC_CS - QC.*pdQP_CS)./QP.^2; % pdQC_CS./QP;
		pdC2P_fProtAOpt = (QP.*pdQC_fProtAOpt - QC.*pdQP_fProtAOpt)./QP.^2; % pdQC_fProtAOpt./QP;
        pdC2P_AOpt = (QP.*pdQC_AOpt - QC.*pdQP_AOpt)./QP.^2; % pdQC_AOpt./QP;
        pdC2P_MOpt = (QP.*pdQC_MOpt - QC.*pdQP_MOpt)./QP.^2;
        pdC2P_rOpt = (QP.*pdQC_rOpt - QC.*pdQP_rOpt)./QP.^2; % -pdQP_rOpt.*QC./QP.^2;
		pdC2P_CProt = (QP.*pdQC_CProt - QC.*pdQP_CProt)./QP.^2; % pdQC_CProt./QP;
		pdC2P_CM = (QP.*pdQC_CM - QC.*pdQP_CM)./QP.^2; % pdQC_CM./QP;
		pdC2P_kST = (QP.*pdQC_kST - QC.*pdQP_kST)./QP.^2; % pdQC_kST./QP;
		pdC2P_PhiS = (QP.*pdQC_PhiS - QC.*pdQP_PhiS)./QP.^2; % pdQC_PhiS./QP;
		pdC2P_PE = -pdQP_PE.*QC./QP.^2;
		pdC2P_gammaDNA = -pdQP_gammaDNA.*QC./QP.^2;
		pdC2P_PDNA = -pdQP_PDNA.*QC./QP.^2;
		%pdC2P_PLip = -pdQP_PLip.*QC./QP.^2;
		pdC2P_DIP = -pdQP_DIP.*QC./QP.^2;

        dC2P_dalphaPLip = -pdQP_alphaPLip.*QC./QP.^2;
        dC2P_dfStorage = -pdQP_fStorage.*QC./QP.^2;     % derivative w.r.t. fStorage
        dC2P_dPLipPCutoff = -pdQP_PLipPCutoff.*QC./QP.^2;
        dC2P_dPLipscale = -pdQP_PLipscale.*QC./QP.^2;
        dC2P_drCutoff = -pdQP_rCutoff.*QC./QP.^2;
        dC2P_dPStorscale = -pdQP_PStorscale.*QC./QP.^2;


    % Now we need to know the derivatives of E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kST,PhiS,fRibE,PRib,gammaDNA,PDNA with respect to the variables we are changing
        pdrOpt_E = rOpt.^2./alphaS.*CI;
        pdrOpt_CI =rOpt.^2./alphaS.*EOpt;
		pdrOpt_alphaS = 1./(1-gammaS - CI.*EOpt);
		pdrOpt_gammaS = rOpt.^2/alphaS;

        pdAOpt_CI = -0.5*EOpt;  pdMOpt_CI = pdAOpt_CI;
        pdAOpt_E = -0.5*CI;     pdMOpt_E = pdAOpt_E;
        pdAOpt_gammaS = -0.5;   pdMOpt_gammaS = -0.5;

    % derivate w.r.t. Q10Photo
        dE_dQ10Photo = dE_dCI.*dCI_dQ10Photo;
        dfProtAOpt_dQ10Photo = dfProtAOpt_dE.*dE_dQ10Photo + dfProtAOpt_dCI.*dCI_dQ10Photo;
        drOpt_dQ10Photo = pdrOpt_E.*dE_dQ10Photo + pdrOpt_CI.*dCI_dQ10Photo;
        dAOpt_dQ10Photo = pdAOpt_E.*dE_dQ10Photo + pdAOpt_CI.*dCI_dQ10Photo;
        dMOpt_dQ10Photo = pdMOpt_E.*dE_dQ10Photo + pdMOpt_CI.*dCI_dQ10Photo;

		dC2P_dQ10Photo = pdC2P_E.*dE_dQ10Photo + pdC2P_CI.*dCI_dQ10Photo + pdC2P_fProtAOpt.*dfProtAOpt_dQ10Photo + pdC2P_AOpt.*dAOpt_dQ10Photo + pdC2P_MOpt.*dMOpt_dQ10Photo + pdC2P_rOpt.*drOpt_dQ10Photo;
        
	% derivative w.r.t. alphaS
        dAOpt_dalphaS = pdAOpt_E.*dE_dalphaS;
		dMOpt_dalphaS = pdMOpt_E.*dE_dalphaS;
		dfProtAOpt_dalphaS = dfProtAOpt_dE.*dE_dalphaS;
		drOpt_dalphaS = pdrOpt_E.*dE_dalphaS + pdrOpt_alphaS;

        dC2P_dalphaS = pdC2P_E.*dE_dalphaS +pdC2P_AOpt.*dAOpt_dalphaS +pdC2P_MOpt.*dMOpt_dalphaS +pdC2P_fProtAOpt.*dfProtAOpt_dalphaS +pdC2P_rOpt.*drOpt_dalphaS;

	% derivative w.r.t fRibE
		dC2P_dfRibE = pdC2P_PE.*dPE_dfRibE + pdC2P_CE.*dCE_dfRibE + pdC2P_E.*dE_dfRibE + pdC2P_fProtAOpt.*pdfProtAOpt_fRibE +pdC2P_AOpt.*pdAOpt_E.*dE_dfRibE +pdrOpt_E.*dE_dfRibE;

	% derivative wrt kST0
		drOpt_dkST = pdrOpt_E.*dE_dkST + pdrOpt_CI.*dCI_dkST;
		dAOpt_dkST = pdAOpt_E.*dE_dkST + pdAOpt_CI.*dCI_dkST;
		dMOpt_dkST = pdMOpt_E.*dE_dkST + pdMOpt_CI.*dCI_dkST;
		dC2P_dkST0 = (pdC2P_kST + pdC2P_E.*dE_dkST + pdC2P_CI.*dCI_dkST + pdC2P_AOpt.*dAOpt_dkST + pdC2P_MOpt.*dMOpt_dkST + pdC2P_fProtAOpt.*dfProtAOpt_dkST + pdC2P_rOpt.*drOpt_dkST).*dkST_dkST0;

	% derivative w.r.t. DIP
		dAOpt_dDIP = pdAOpt_E.*dE_dDIP;
		dMOpt_dDIP = pdMOpt_E.*dE_dDIP;
		drOpt_dDIP = pdrOpt_E.*dE_dDIP;
		dC2P_dDIP = pdC2P_DIP + pdC2P_E.*dE_dDIP + pdC2P_fProtAOpt.*dfProtAOpt_dDIP + pdC2P_AOpt.*dAOpt_dDIP + pdC2P_MOpt.*dMOpt_dDIP + pdC2P_rOpt.*drOpt_dDIP;

   % derivative w.r.t. gammaDNA
        dAOpt_dgammaDNA = pdAOpt_E.*dE_dgammaDNA + pdAOpt_gammaS.*dgammaS_dgammaDNA;
		dMOpt_dgammaDNA = pdMOpt_E.*dE_dgammaDNA + pdMOpt_gammaS.*dgammaS_dgammaDNA;
		drOpt_dgammaDNA = pdrOpt_E.*dE_dgammaDNA + pdrOpt_gammaS.*dgammaS_dgammaDNA;

        dC2P_dgammaDNA = pdC2P_gammaDNA + pdC2P_gammaS.*dgammaS_dgammaDNA + pdC2P_E.*dE_dgammaDNA + pdC2P_AOpt.*dAOpt_dgammaDNA +pdC2P_MOpt.*dMOpt_dgammaDNA +pdC2P_fProtAOpt.*dfProtAOpt_dgammaDNA +pdC2P_rOpt.*drOpt_dgammaDNA + pdC2P_CS.*dCS_dgammaDNA;


%% *** Calculate RQ derivatives ********
%RQtotalO2toC = (QC + 0.25*QH - 0.5*QOxy +1.25*QN + 1.25*QP + 1.5*QSulf)./QC;

%----- Calculate derivatives of RQ -----
% move section to later in code 
pdRQ_QC = (QC - (QC + 0.25*QH - 0.5*QOxy +1.25*QN + 1.25*QP + 1.5*QSulf))./QC.^2;
pdRQ_QH = 0.25./QC;
pdRQ_QOxy = -0.5./QC;
pdRQ_QN = 1.25./QC;
pdRQ_QP = 1.25./QC;
pdRQ_QSulf = 1.5./QC;

% QH derivatives (from symbolic toolbox
%%% MISSING COMPONENT FROM QHPLip
pdQH_E = (HE + HL*(CI - 1))/molarH + (12*kST.*(PhiS + 1).*(CS*gammaS + AOpt*CM + CE*EOpt + AOpt.*CProt.*fProtAOpt + CL.*EOpt.*(CI - 1)))/molarC + (12*EOpt.*kST.*(PhiS + 1).*(CE + CL.*(CI - 1)))/molarC;
pdQH_CI = (12*CL*kST.*(PhiS + 1).*EOpt.^2)/molarC + (HL*EOpt)/molarH ;
pdQH_fProtAOpt = (AOpt*HProt)/molarH + (12*AOpt.*CProt.*EOpt.*kST.*(PhiS + 1))/molarC ;
pdQH_AOpt = (HProt*fProtAOpt)/molarH + (12*CProt*EOpt.*fProtAOpt.*kST*(PhiS + 1))/molarC ;
pdQH_MOpt = HM/molarH + (12*CM*EOpt.*kST*(PhiS + 1))/molarC ;
pdQH_rOpt = 0;

% pdQN 
pdQN_E = (NE + NL*(CI - 1))/molarN; 
pdQN_CI = (EOpt*NL)/molarN ;
pdQN_fProtAOpt = (AOpt*NProt)/molarN;
pdQN_AOpt = (NProt*fProtAOpt)/molarN;
pdQN_MOpt = NM/molarN;
pdQN_rOpt = 0;

% pdQOxy 
%%% MISSING COMPONENT FROM QOxyPStor an QOxyPLip
pdQOxy_E = (OE + OL*(CI - 1))/molarO + (6*kST.*(PhiS + 1).*(CS*gammaS + CE*EOpt + CM*MOpt + AOpt.*CProt.*fProtAOpt + CL*EOpt.*(CI - 1)))/molarC + (6*EOpt.*kST.*(PhiS + 1).*(CE + CL.*(CI - 1)))/molarC;

% IN PROGRESS
% pdQOxyfunctional_E = (OE + OL*(CI - 1))/molarO;
% pqQOxyreserve_E = pdQCreserve_E;
% pdQOxyreserve_E = (6*kST.*(PhiS + 1).*(CS*gammaS + CE*EOpt + CM*MOpt + AOpt.*CProt.*fProtAOpt + CL*EOpt.*(CI - 1)))/molarC + (6*EOpt.*kST.*(PhiS + 1).*(CE + CL.*(CI - 1)))/molarC;
% pdQOxyPStor_E = (5000*O2P_PStor*P*fStorage*( ...
%     (CE + CL*(CI - 1))/molarC + (6*kST*(PhiS + 1)*(CS*gammaS + CE*EOpt + CM*MOpt + AOpt*CProt*fProtAOpt + CL*EOpt*(CI - 1)))/molarC +...
%     (6*EOpt*kST*(PhiS + 1)*(CE + CL*(CI - 1)))/molarC)...
%     )/(exp(PStor_scale*(PStor_rCutoff - rOpt)) + 1);
% pdQOxyPLip_E = (MOpt*O2P_PLip*PPhospholipid*alphaPLip*molarC*((CE + CL*(CI - 1))/molarC + ...
%     (6*kST*(PhiS + 1)*(CS*gammaS + CE*EOpt + CM*MOpt + AOpt*CProt*fProtAOpt + CL*EOpt*(CI - 1)))/molarC + ...
%     (6*EOpt*kST*(PhiS + 1)*(CE + CL*(CI - 1)))/molarC)...
%     )/(CPhospholipid*molarP*(exp(-PLip_scale*(P - PLip_PCutoff)) + 1));
% pdQOxy_E = pdQOxyfunctional_E + pdQOxyreserve_E + pdQOxyPStor_E + pdQOxyPLip_E;
%

pdQOxy_CI = (6*CL.*kST.*(PhiS + 1).*EOpt.^2)/molarC + (OL*EOpt)/molarO ;
pdQOxy_fProtAOpt = (AOpt*OProt)/molarO + (6*AOpt.*CProt.*EOpt.*kST.*(PhiS + 1))./molarC ;
pdQOxy_AOpt = (OProt*fProtAOpt)/molarO + (6*CProt.*EOpt.*fProtAOpt.*kST.*(PhiS + 1))./molarC;
pdQOxy_MOpt = OM/molarO + (6*CM.*EOpt.*kST.*(PhiS + 1))./molarC ;
pdQOxy_rOpt = 0;

% pdQSulf
pdQSulf_E = (SulfE + SulfL*(CI - 1))/molarS ;
pdQSulf_CI = (EOpt*SulfL)/molarS ;
pdQSulf_fProtAOpt = (AOpt*SulfProt)/molarS ;
pdQSulf_AOpt = (SulfProt*fProtAOpt)/molarS ;
pdQSulf_MOpt = SulfM/molarS ;
pdQSulf_rOpt = 0;

%dQuotas/dQ10Photo
dQC_dQ10Photo = pdQC_E.*dE_dQ10Photo + pdQC_CI.*dCI_dQ10Photo + pdQC_fProtAOpt.*dfProtAOpt_dQ10Photo + pdQC_AOpt.*dAOpt_dQ10Photo + pdQC_MOpt.*dMOpt_dQ10Photo + pdQC_rOpt.*drOpt_dQ10Photo;
dQH_dQ10Photo = pdQH_E.*dE_dQ10Photo + pdQH_CI.*dCI_dQ10Photo + pdQH_fProtAOpt.*dfProtAOpt_dQ10Photo + pdQH_AOpt.*dAOpt_dQ10Photo + pdQH_MOpt.*dMOpt_dQ10Photo + pdQH_rOpt.*drOpt_dQ10Photo;
dQN_dQ10Photo = pdQN_E.*dE_dQ10Photo + pdQN_CI.*dCI_dQ10Photo + pdQN_fProtAOpt.*dfProtAOpt_dQ10Photo + pdQN_AOpt.*dAOpt_dQ10Photo + pdQN_MOpt.*dMOpt_dQ10Photo + pdQN_rOpt.*drOpt_dQ10Photo;
dQOxy_dQ10Photo = pdQOxy_E.*dE_dQ10Photo + pdQOxy_CI.*dCI_dQ10Photo + pdQOxy_fProtAOpt.*dfProtAOpt_dQ10Photo + pdQOxy_AOpt.*dAOpt_dQ10Photo + pdQOxy_MOpt.*dMOpt_dQ10Photo + pdQOxy_rOpt.*drOpt_dQ10Photo;
dQP_dQ10Photo = pdQP_E.*dE_dQ10Photo + pdQP_CI.*dCI_dQ10Photo + pdQP_fProtAOpt.*dfProtAOpt_dQ10Photo + pdQP_AOpt.*dAOpt_dQ10Photo + pdQP_MOpt.*dMOpt_dQ10Photo + pdQP_rOpt.*drOpt_dQ10Photo;
dQSulf_dQ10Photo = pdQSulf_E.*dE_dQ10Photo + pdQSulf_CI.*dCI_dQ10Photo + pdQSulf_fProtAOpt.*dfProtAOpt_dQ10Photo + pdQSulf_AOpt.*dAOpt_dQ10Photo + pdQSulf_MOpt.*dMOpt_dQ10Photo + pdQSulf_rOpt.*drOpt_dQ10Photo;

dRQ_dQ10Photo = pdRQ_QC.*dQC_dQ10Photo + pdRQ_QH.*dQH_dQ10Photo + pdRQ_QN.*dQN_dQ10Photo +pdRQ_QOxy.*dQOxy_dQ10Photo +pdRQ_QP.*dQP_dQ10Photo + pdRQ_QSulf.*dQSulf_dQ10Photo;

% 
% %method 2 
% % RQ = x/y
% % pdRQ_E = (y*dx_dE - x*dy_dE)/y^2
% tmp.x = (QC + 0.25*QH - 0.5*QOxy +1.25*QN + 1.25*QP + 1.5*QSulf);
% tmp.y = QC;
% tmp.dx_E = (pdQC_E + 0.25*pdQH_E -0.5*pdQOxy_E +1.25*pdQN_E +1.25*pdQP_E +1.5*pdQSulf_E);
% tmp.dy_E = pdQC_E;
% pdRQ_E = (QC*(pdQC_E + 0.25*pdQH_E -0.5*pdQOxy_E +1.25*pdQN_E +1.25*pdQP_E +1.5*pdQSulf_E) - (QC + 0.25*QH - 0.5*QOxy +1.25*QN + 1.25*QP + 1.5*QSulf)*pdQC_E )/(QC.^2)
% pdRQ_CI
% pdRQ_fProtA
% pdRQ_AOpt
% pdRQ_MOpt
% pdRQ_rOpt
% dRQ_Q10Photo = pdRQ_E*dE_Q10Photo + pdRQ_CI*dCI_dQ10Photo + pdRQ_fProtAOpt.*dfProtAOpt_dQ10Photo + pdRQ_AOpt.*dAOpt_dQ10Photo + pdRQ_MOpt.*dMOpt_dQ10Photo + pdRQ_rOpt.*drOpt_dQ10Photo;
%    

%% Check growth rates (only needed for debugging) 
L = CI.*EOpt - EOpt;
gdrymass = 4/3*pi*rOpt.^3*pDry*rho; % dry mass of cell [grams]
muL = alphaI.*L./(1+PhiS);
muN = 4*pi*(rOpt./10^6).*DN.*(N*10^3)./(gdrymass.*QN); %convert units of r to m and N to mol/m^3 => units: 1/(hr*g)
muP = 4*pi*(rOpt./10^6).*DP.*(P*10^3).*fProtAOpt./(gdrymass.*QP);

%% Output
    out.CP = CP;
    out.CPnostor = QC./QPfunctional;
    out.NP = NP;
    out.CN = CN;
    out.LimType = LimType;
    out.r = rOpt;
    out.E = EOpt;
	out.A = AOpt;
    out.mu = muOpt;
	out.QP = QP;
	out.QC = QC;
	out.PStor = PStor;
    out.PLip = PLip;
	out.CI = CI;
	out.L = CI.*EOpt - EOpt;
	out.muN = muN;
	out.muP = muP;
	out.muL = muL;
    out.RQtotalO2toC = RQtotalO2toC ;

    % name a second structure for derivatives C2Px?
	out.dC2P_dQ10Photo = dC2P_dQ10Photo;
	out.dC2P_dfStorage = dC2P_dfStorage;
	out.dC2P_dPCutoff = dC2P_dPLipPCutoff;
	out.dC2P_drCutoff = dC2P_drCutoff;
    out.dC2P_dalphaPLip = dC2P_dalphaPLip;
    out.dC2P_dPLipscale = dC2P_dPLipscale;
    out.dC2P_dPStorscale = dC2P_dPStorscale;
	out.dC2P_dalphaS = dC2P_dalphaS;
	out.dC2P_dfRibE = dC2P_dfRibE;
	out.dC2P_dkST0 = dC2P_dkST0;
	out.dC2P_dDIP = dC2P_dDIP;
    out.dC2P_dgammaDNA = dC2P_dgammaDNA;

	out.dE_dDIP = dE_dDIP;
    % RQ derivatives
    out.dRQ_dQ10Photo = dRQ_dQ10Photo;

%% define functions
    function alphaPhoto = alphaPhoto(Irr)
        % This is a fourth order rational Chebyshev polynomial approximation
        % of the function [F1,alphaPhoto] = CellPhoto(Irr)
        % The approximation works well (except perhaps at I=0)
        % Using this interpolating function, based on Chebyshev polynomials,
        % allows us to avoid repeatedly solving a nonlinear system.
        %
	    % To run the full function, replace with: 
        % [F1opt,alphaPhoto] = CellPhoto(Irr)

        % CellPhoto calculates the carbon gathering efficiency of 1 unit of
        % the photosynthetic apparatus. It is based on a paper by Geider and
        % Talmy. This model involves optimizing photosynthetic allocations
        % between carbon fixation and photo capture/electron transport.
        %
        % This function models photosynthesis by assuming that the
        % photosynthetic apparatus is divided into two compartments F1 stands
        % for the carbon fixing compartment and F2 stands for the light
        % harvesting compartment. Associated with each of these compartments are
        % rate constants k1 and k2.
        %
        % The maximum photosynthesis rate is: Pm = min([k1*F1,k2*F2]).
        % The actual photosynthesis rate depends on the irradiance and the
        % investments F1 and F2 according to f_Photo = Pm*(1-exp(-alpha_ph *
        % PhiM * F2 * Irr / PMax))
        % We assume that F1 + F2 = 1.0, and that the cell allocates between F1
        % and F2 optimally to maximize f_Photo Thus f_Photo^*(I) = max_F1
        % f_Photo(F1,1-F1) 

        I0 = Irr/200.0; % scale by roughly magnitude of irradiance (use avg irradience of box, assuming exponential decay)
        L = 2.0;   % Length scale
        ITrans = (I0-L)./(I0+L); % transforms interval to -1 to 1 . interpolation works on this variabble

        Ch0 = 0.19590903;
        Ch1 = 0.10462839;
        Ch2 = -0.05337273;
        Ch3 = 0.01467971;

        alphaPhoto = Ch0+Ch1*(ITrans)+Ch2*(2*ITrans.^2-1)+Ch3*(4*ITrans.^3-3*ITrans);
    end

	function [EColim, dEColim_dCI, dEColim_dalphaS, dEColim_dfRibE,pdEColim_kST,dEColim_daP,dEColim_dgammaDNA] = EColim(EMin,EMax,i) % rhsFuncColim(EMin,EMax)
		CIi = CI(i);
		kSTi = kST(i);
		aNi = aN(i);
		aPi = aP(i);
        % polynomial coefficients
		c3 = kSTi^2*NProt*PE + (aPi*aNi/alphaS^4)*CIi^3 - (aPi/alphaS^2)*kSTi*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		% create function to return all derivs of c3 wrt params?
		c2 = kSTi^2*NProt*gammaS*PS - 3*(aPi*aNi/alphaS^4)*CIi^2*(1-gammaS) + (aPi/alphaS^2)*kSTi*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) - (aPi/alphaS^2)*kSTi*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		c1 = 3*(aPi*aNi/alphaS^4)*CIi*(1-gammaS)^2 + (aPi/alphaS^2)*kSTi*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		c0 = -(aPi*aNi/alphaS^4)*(1-gammaS)^3;

		E0 = (EMin+EMax)/2; % starting point for solver

		if ~all(isfinite([c3,c2,c1,c0,E0]))
			fprintf('Error in EColim. Complex values. \n')
        end
		% solve for the zeros of rhsFuncColim
        % rhsFuncColim = @(E) (c3*E^3 + c2*E^2 + c1*E +c0);

        [sol,exitflag] = complex_cubic_zero(c3,c2,c1,c0,E0);
		if exitflag <=0
            errcount = errcount + 1;
        end
		EColim = sol;

		%%% dE_dDIP
		%dE_daP
		dc3_daP = (aNi/alphaS^4)*CIi^3 - (1/alphaS^2)*kSTi*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		dc2_daP = - 3*(aNi/alphaS^4)*CIi^2*(1-gammaS) + (1/alphaS^2)*kSTi*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) - (1/alphaS^2)*kSTi*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc1_daP = 3*(aNi/alphaS^4)*CIi*(1-gammaS)^2 + (1/alphaS^2)*kSTi*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc0_daP = -(aNi/alphaS^4)*(1-gammaS)^3;
		dEColim_daP = -(EColim^3*dc3_daP+EColim^2*dc2_daP+EColim*dc1_daP+dc0_daP)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		%%% dE_dCI
		dc3_dCI = 3*(aPi*aNi/alphaS^4)*CIi^2 - (aPi/alphaS^2)*kSTi*(NE+(CIi-1)*NL-0.5*CIi*NM) - (aPi/alphaS^2)*kSTi*CIi*(NL-0.5*NM);
		dc2_dCI = -6*(aPi*aNi/alphaS^4)*CIi*(1-gammaS) + (aPi/alphaS^2)*kSTi*(1-gammaS)*(NL-0.5*NM) - (aPi/alphaS^2)*kSTi*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc1_dCI = 3*(aPi*aNi/alphaS^4)*(1-gammaS)^2;
		dc0_dCI = 0;
		dEColim_dCI = -(EColim^3*dc3_dCI+EColim^2*dc2_dCI+EColim*dc1_dCI+dc0_dCI)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		% pdE_kST   %
		pdc3_kST = 2*kSTi*NProt*PE - (aPi/alphaS^2)*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		pdc2_kST = 2*kSTi*NProt*gammaS*PS + (aPi/alphaS^2)*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) -(aPi/alphaS^2)*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		pdc1_kST = (aPi/alphaS^2)*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		pdc0_kST = 0;
		pdEColim_kST = -(EColim^3*pdc3_kST+EColim^2*pdc2_kST+EColim*pdc1_kST+pdc0_kST)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		% dalphaS
		dc3_dalphaS = -4*(aPi*aNi/alphaS^5)*CIi^3 - -2*(aPi/alphaS^3)*kSTi*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		dc2_dalphaS = -4*(-3*aPi*aNi/alphaS^5)*CIi^2*(1-gammaS) + -2*(aPi/alphaS^3)*kSTi*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) - -2*(aPi/alphaS^3)*kSTi*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc1_dalphaS = -4*3*(aPi*aNi/alphaS^5)*CIi*(1-gammaS)^2 + -2*(aPi/alphaS^3)*kSTi*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc0_dalphaS = 4*(aPi*aNi/alphaS^5)*(1-gammaS)^3;
		dEColim_dalphaS = -(EColim^3*dc3_dalphaS+EColim^2*dc2_dalphaS+EColim*dc1_dalphaS+dc0_dalphaS)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		% dfRibE
		dc3_dPE = kSTi^2*NProt;
		dEColim_dPE = -(EColim^3*dc3_dPE)/(3*EColim^2*c3 +2*EColim*c2 +c1);
		dc3_dNE = -(aPi/alphaS^2)*kSTi*CIi;
		dc2_dNE = (aPi/alphaS^2)*kSTi*(1-gammaS);
		dc1_dNE = 0;
		dc0_dNE = 0;
		dEColim_dNE = -(EColim^3*dc3_dNE+EColim^2*dc2_dNE)/(3*EColim^2*c3 +2*EColim*c2 +c1);
		dEColim_dfRibE = dEColim_dPE.*dPE_dfRibE + dEColim_dNE.*dNE_dfRibE;

        % dgammaDNA
        dc3_dgammaDNA = 0;

        pdc2_gammaS = kSTi^2*NProt*PS + 3*(aPi*aNi/alphaS^4)*CIi^2 - (aPi/alphaS^2)*kSTi*(NE+(CIi-1)*NL-0.5*CIi*NM) - (aPi/alphaS^2)*kSTi*CIi*(NS-0.5*NM);
        pdc2_PS = kSTi^2*NProt*gammaS;
        pdc2_NS = - (aPi/alphaS^2)*kSTi*CIi*(gammaS);
        dc2_dgammaDNA = pdc2_gammaS*1 + pdc2_PS*dPS_dgammaDNA + pdc2_NS*dNS_dgammaDNA;

        pdc1_gammaS = -2*3*(aPi*aNi/alphaS^4)*CIi*(1-gammaS) + (aPi/alphaS^2)*kSTi*((1-gammaS)*NS-gammaS*NS-(1-gammaS)*NM);
        pdc1_NS = (aPi/alphaS^2)*kSTi*(1-gammaS)*(gammaS);
        dc1_dgammaDNA = pdc1_gammaS*1 + pdc1_NS*dNS_dgammaDNA;

        pdc0_gammaS = 3*(aPi*aNi/alphaS^4)*(1-gammaS)^2;
        dc0_dgammaDNA = pdc0_gammaS*1;

        dEColim_dgammaDNA = -(EColim^3*dc3_dgammaDNA+EColim^2*dc2_dgammaDNA+EColim*dc1_dgammaDNA+dc0_dgammaDNA)/(3*EColim^2*c3 +2*EColim*c2 +c1);
    end


    function [fProtAColim, dfProtA_dCI, dfProtA_dE, dfProtA_dfRibE,pdfProtA_DIP,pdfProtA_gammaS,dfProtA_dNS,dfProtA_dPS]  =AColim(EOpt_i,CIi,i)
		% fProtAColim: fraction periplasm (A) devoted to protein, for P&N colimited case
        % uptake rate of phosphorous = linear function of how much protein in periplasm
        
        %solving for fProtAColim
        coefA = (1/2)*(1-gammaS-CIi.*EOpt_i)*NProt;
        coefB = (1/2)*(1-gammaS-CIi.*EOpt_i)*NM + EOpt_i*NE + (CIi-1)*EOpt_i*NL + gammaS*NS;
        coefC = -DN(i)/DP(i)*molarN/molarP*N(i)/P(i)*(EOpt_i*PE+gammaS*PS);

        fProtAColim = (-coefB+sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);
        if fProtAColim <0 | fProtAColim>1
            tmp = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);
            if tmp >0 & tmp<=1
                fProtAColim = tmp;
                fprintf(' fProtAColim = tmp \n')
            end
        end

        %%% derivatives of AColim w.r.t. parameters
		%pdfProtA_dDIP
		pdcoefA_DIP = 0;
		pdcoefB_DIP = 0;
		pdcoefC_DIP = coefC/(-P(i));
		pdfProtA_DIP = ( (2*coefA)*(-pdcoefB_DIP+0.5/sqrt(coefB^2-4*coefA*coefC)*(2*coefB*pdcoefB_DIP-4*coefA*pdcoefC_DIP-4*coefC*pdcoefA_DIP)) - (-coefB+sqrt(coefB^2-4*coefA*coefC))*2*pdcoefA_DIP ) / (2*coefA)^2;

		%dfProtA_dCI
        pdcoefA_CI = -(1/2)*EOpt_i*NProt;
        pdcoefB_CI = -(1/2)*EOpt_i*NM +EOpt_i*NL;
        pdcoefC_CI = 0;
		dfProtA_dCI = ( (2*coefA)*(-pdcoefB_CI+0.5/sqrt(coefB^2-4*coefA*coefC)*(2*coefB*pdcoefB_CI-4*coefA*pdcoefC_CI-4*coefC*pdcoefA_CI)) - (-coefB+sqrt(coefB^2-4*coefA*coefC))*2*pdcoefA_CI ) / (2*coefA)^2;

		%dfProtA_dE
        pdcoefA_E = -(1/2)*CIi*NProt;
        pdcoefB_E = -(1/2)*CIi*NM +NE +(CIi-1)*NL;
        pdcoefC_E = -DN(i)/DP(i)*molarN/molarP*N(i)/P(i) *PE;
        dfProtA_dE = (2*coefA.*(-pdcoefB_E + 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*pdcoefB_E-4*(pdcoefA_E.*coefC+pdcoefC_E*coefA))) -(-coefB+sqrt(coefB.^2-4*coefA.*coefC)).*2.*pdcoefA_E)./(4*coefA.^2);
		
        %dfRibE;
		pdcoefA_NE = 0;
		pdcoefB_NE = EOpt_i;
		dfProtA_dNE = (-pdcoefB_NE+0.5/sqrt(coefB^2-4*coefA*coefC)*(2*coefB*pdcoefB_NE)) / (2*coefA);
		pdcoefC_PE = -DN(i)/DP(i)*molarN/molarP*N(i)/P(i)*(EOpt_i);
		dfProtA_dPE =(0.5/sqrt(coefB.^2-4*coefA.*coefC).*(-4*coefA*pdcoefC_PE))./(2*coefA);
		dfProtA_dfRibE = dfProtA_dNE*dNE_dfRibE + dfProtA_dPE*dPE_dfRibE;

        %dgammaS
        pdcoefA_gammaS = -(1/2)*NProt;
        pdcoefB_gammaS = -(1/2)*NM + NS;
        pdcoefC_gammaS = -DN(i)/DP(i)*molarN/molarP*N(i)/P(i)*PS;
        pdfProtA_gammaS = (2*coefA.*(-pdcoefB_gammaS + 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*pdcoefB_gammaS-4*(pdcoefA_gammaS.*coefC+pdcoefC_gammaS*coefA))) -(-coefB+sqrt(coefB.^2-4*coefA.*coefC)).*2.*pdcoefA_gammaS)./(4*coefA.^2);

        %dNS
        pdcoefB_NS = gammaS;
        dfProtA_dNS = (-pdcoefB_NS+0.5/sqrt(coefB^2-4*coefA*coefC)*(2*coefB*pdcoefB_NS)) / (2*coefA);

        %dPS
        pdcoefC_PS = -DN(i)/DP(i)*molarN/molarP*N(i)/P(i)*gammaS;
        dfProtA_dPS = ( (2*coefA)*(0.5/sqrt(coefB^2-4*coefA*coefC)*(0-4*coefA*pdcoefC_PS-0)) - 0 ) / (2*coefA)^2;
    end

end % end of CellCNP
