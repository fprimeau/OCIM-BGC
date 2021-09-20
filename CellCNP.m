function [out,M] = CellCNP(par,x,P,N,T,Irr)
  % DESCRIPTION:
  % Previously named CPDual
  % This function calculates the C:P ratio of the maximally growing plankton
  % type in the environment specified by the function arguments
	% Currently: only relies on observed PO4,NO3,Temp,&PAR fields

  % INPUTS:
  % Irr is light level (PAR, not total irradiance) in microeinsteins m^-2 s^-1 (i.e. umol photons m^-2 s^-1)
	% %%% PAR, but with units of photosynthetic photon flux density (PPFD) [μmole photons m-2 s-1]
  % T is temperature in Celsius
  % P is phosphate concentration in mol/L (OCIM-BGC in mmol/m^3. function input was converted to mol/L in eqC2Puptake)
  % N is nitrate concentration in mol/L (same not for OCIM-BGC)

  % AUTHOR: George Hagstrom. MODIFIED BY: Megan Sullivan

    on = true; off = false;
	% unpack the parameters to be optimized
	M = par.BIO;
	nbx = 0; % count number of C model tunable parameters;

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
		M.fRibE = 0.5*(1+tanh(tfRibE)); % change fRibE = 0.5*(1+tanh(tfRibE))
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

  % Read in parameter values from par
		 %%% optimize these first %%%
     Q10Photo = M.Q10Photo;
     fRibE = M.fRibE;
     fStorage = M.fStorage;
	 kST0 = M.kST0;

	 alphaS = M.alphaS;
     gammaDNA = M.gammaDNA;
     gammaLipid = M.gammaLipid;
     %lPCutoff = M.lPCutoff;
     PLip_PCutoff = M.PLip_PCutoff;
     PLip_scale = M.PLip_scale;
     %r0Cutoff = M.r0Cutoff;
     PStor_rCutoff = M.PStor_rCutoff;
	 PStor_scale = M.PStor_scale; %(replaces CStor)
		%
     DNT0 = M.DNT0;
     DPT0 = M.DPT0;
     Q10Diffusivity = M.Q10Diffusivity;
     AMin = M.AMin;
     %CStor = M.CStor;
     PhiS = M.PhiS;
		%
	 alphaPLip = M.alphaPLip;
		%
		% %%% everything below here should remain fixed %%%
     pDry = M.pDry;  % fixed
     rho = M.rho;     % fixed
     fProtM = M.fProtM;
     fProtL = M.fProtL;
     % Define the P, N, and C content of each cellular component (protein, lipid, ribosome, etc.)
     PDNA = M.PDNA; % fixed
     PRib = M.PRib; % fixed
     PPhospholipid = M.PPhospholipid; % fixed
     NProt = M.NProt;  % fixed
	 NRib = M.NRib;  %fixed
     NDNA = M.NDNA; % fixed
     CProt = M.CProt; % fixed
     CDNA = M.CDNA; % fixed
     CPhospholipid = M.CPhospholipid; % fixed
     CLipid = M.CLipid;  % fixed
	 CRib = M.CRib;
    %
	par.BIO = M;

	nanindx = find(isnan(Irr) | isnan(T) | isnan(N) |isnan(P));
	if ~isempty(nanindx)
		fprintf('Warning: [ %i ] NaNs found in Cell Model input fields. \n',length(nanindx))
    end

    %set negative phosphate values to smallest positive concentration.
    %(negative values mess up the code)
	% moved to eqC2Puptake
    % fprintf('replacing %d negative Phosphate concentrations with the minimum positive concentration \n',length(P(P<0)))
	% negDIPindx = (P<0);
    % P(P<0)= real(min(P(P>=0))); %imaginary part = 0, so these points dont affect deriv


%% Define constants
	    molarC = 12.0;              % molar mass of carbon [g/mol]
	    molarN = 14.0;              % molar mass of nitrogen [g/mol]
	    molarP = 31.0;              % molar mass of phosphorus [g/mol]
		T0 = 25.0;                  % standard temperature [degC]


%%% optimizable Parameters %%%
  % George performed a global sensitivity analysis to determine which parameters were the most important
  % Key parameters are Q10Photo, fRibE, kST0, fStorage
  % Q10Photo determines whether "translation compensation" happens,
	% fRibE and kST0 impact the growth hypothesis (it is probably best to fix one of them to avoid an odd looking posterior),
  % and fStorage impacts frugality
  % secondary parameters are alphaS, gammaDNA (and the other gamma factors), lPCutoff,


	gammaS = gammaDNA+gammaLipid;   % mass fraction of non-membrane or periplasm structural components (gammaDNA+gammaLipid+gammaCarb)
    rFullA = alphaS/(2.0*AMin);
						 % rFullA is the radius at which if the entire periplasm was filled with
						 % uptake proteins, the percent biomass would be AMin. For cells above this
						 % size, I treat the periplasmic space as being 100% filled with P-uptake
						 % proteins. When AMin = 1, it basically forces all cells to totally invest
						 % in periplasmic proteins, so it disables the part of the model where a
						 % cell can have less uptake proteins to reduce N-quota when P is not too
						 % low.

    function alphaPhoto = alphaPhoto(Irr)
       % where is this function from? (i.e. source?)
			 % This code is a clever way to avoid a nonlinear solve which allows code to run faster:

		% alphaPhoto=F1. also computed by function CellPhoto(Irr).
		% F2 = 1-alphaPhoto is proportional to chlorphyll content of cell.

        % This function calculates the carbon gathering efficiency of 1 unit of
        % the photosynthetic apparatus. It is based on a paper by Geider and
        % Talmy. This model involves optimizing photosynthetic allocations
        % between carbon fixation and photo capture/electron transport. Since I
        % don't want to solve this nonlinear system all the time, I found an
        % interpolating function based on Chebyshev polynomials that captures
        % the pattern quite well (except perhaps at I=0)

        % This function models photosynthesis by assuming that the
        % photosynthetic apparatus is divided into two compartments F1 stands
        % for the carbon fixing compartment and F2 stands for the light
        % harvesting compartment. Associated with each of these compartments are
        % rate constants k1 and k2

        % The maximum photosynthesis rate is: Pm = min([k1*F1,k2*F2]).

        % The actual photosynthesis rate depends on the irradiance and the
        % investments F1 and F2 according to f_Photo = Pm*(1-exp(-alpha_ph *
        % PhiM * F2 * Irr / PMax))

        % We assume that F1 + F2 = 1.0, and that the cell allocates between F1
        % and F2 optimally to maximize f_Photo Thus f_Photo^*(I) = max_F1
        % f_Photo(F1,1-F1) Not wanting to solve this optimization problem every
        % time I solve the model, I calculated f_Photo^*(I) for a large range of
        % values of I and approximated the resulting function using Chebyshev
        % polynomials. Below is a fourth order chebyshev polynomial
        % approximation

        % rational Chebyshev

        I0 = Irr/200.0; % scale by roughly magnitude of irradiance (use avg irradience of box, assuming exponential decay)
        L = 2.0;   % Length of
        ITrans = (I0-L)./(I0+L); % transforms interval to -1 to 1 . interpolation works on this variabble

        Ch0 = 0.19590903;
        Ch1 = 0.10462839;
        Ch2 = -0.05337273;
        Ch3 = 0.01467971;

        alphaPhoto = Ch0+Ch1*(ITrans)+Ch2*(2*ITrans.^2-1)+Ch3*(4*ITrans.^3-3*ITrans);
    end

%keyboard;
%%

% temperature dependent rates
    DN = DNT0*Q10Diffusivity.^((T-T0)./10.0);    % diffusivity of nitrate [ m^2/hr]
    DP = DPT0*Q10Diffusivity.^((T-T0)./10.0);     % diffusivity of phosphate [ m^2/hr]

	Q10Bio = 2.0;
    kST = kST0*Q10Bio.^((T-T0)./10.0);  % specific synthesis rate of biosynthetic apparatus [1/hr]
	dkST_dkST0 = Q10Bio.^((T-T0)./10.0);

    alphaI = alphaPhoto(Irr).*Q10Photo.^((T-T0)./10.0);  % ?
						% alphaI is the amount of carbon generated by 1 unit of
						% photosynthetic apparatus. This is what the alphaPhoto function
						% calculates It is useful to solve for this now to setup solving for
						% the optimum. We multiply by Q10Photo to adjust for the temperature
						% dependence


	CI  = 1+(kST.*(1+PhiS))./alphaI;
                % CI is a constant defined for convenience that is used in finding optimal solutions.
				% It comes from balancing biosynthetic and photosynthetic rates, ie if we set:
				% kST*E = alphaI*L/(1+PhiS), then L = kST*E*(1+PhiS)/alphaI
				% Thus whenever E+L occurs in the equations, we can replace it with E*CI, which makes life easier.

	%CI  = 1+(kST.*(1+PhiS))./ (alphaPhoto(Irr).*Q10Photo.^((T-T0)./10.0)) ;
	%CI = 1+(kST.*(1+PhiS))./alphaPhoto(Irr) .* Q10Photo.^(-(T-T0)./10.0) ;

	dCI_dQ10Photo = (kST.*(1+PhiS))./alphaPhoto(Irr) .* (-(T-T0)./10.0) .*Q10Photo.^((-(T-T0)./10.0)-1) ;
	dCI_dQ10Photo_dQ10Photo= (kST.*(1+PhiS))./alphaPhoto(Irr) .* (-(T-T0)./10.0).* (-(T-T0)./10.0 -1) .*Q10Photo.^((-(T-T0)./10.0)-2) ;

	dCI_dkST = (1+PhiS)./alphaI;

% | cell membrane is made up of lipid and protein
  %  fProtM = 0.25;               % protein fraction of cell membranes
    fLipidM = 1-fProtM;          % Lipid fraction of cell membranes  % not the same as alphaPLip (alphaPLip = biomass of membrane in cell)
	%	alphaPLip = 0.12; 						% fraction of membrane that is phopholipids ???

% | light harvesting apparatus is 70% protein. The rest is Lipid
  %  fProtL = .7;            % protein fraction of light harvesting apparatus
	fLipidL = 1-fProtL; 		% lipid fraction of light harvesting apparatus

% | Biosynthetic apparatus is ribosomes and proteins
	fProtE = 1-fRibE;



%% FUNCTIONAL POOLS
    %%% _L = light/photosynthesis --> (proteins)
    %%% _M = cell membrane?       --> (lipids + proteins)
    %%% _E = biosynthesis         --> (proteins,ribosomes,) %%% That is right
    %%% _S = other structure components (not membrane)? --> (DNA +lipids)
    %%% Ctotal = CL*L + CM*M + CE*E + CS*S_other
    %%% 1 = L + M + E + S
    % periplasmic space (A)??? - prokaryotes have periplasmic space, eukaryotes dont - generally filled with proteins for nutrient uptake. can be a big fraction of mass for small bateria (maybe 40% or more)
      % cells only invest in periplasmic proteins if they need too.
      % exact formula for how periplasm relates to p uptake is too complicated
      %

    NL = fProtL*NProt;      % N fraction of light harvesting apparatus (assumes all N is in proteins)
                              %   = fraction of light apparatus that is protein * fraction of protein that is N
    CL = CProt*fProtL+CLipid*(1-fProtL);            % Carbon fraction of light harvesting apparatus
                              % || why not fProtL*CProt; ? % I think this is a big, should be (CProt*fProtL+CLipid*(1-fProtL)) . This slightly underestimates the amount of C intended to be in the L pool
							  % seperate this into F1 and F2
							  % CL1 = CProt*fProtL1 + CLipid*(1-fProtL1)
							  % CL2 = CChlor
							  % CL = CL1*alphaPhoto + CL2*(1-alphaPhoto)

    NM = fProtM*NProt;        % N fraction of membrane
                              %   = fraction of cell membrane that is protein * fraction of protein that is N
    CM = fProtM*CProt+fLipidM*CPhospholipid;     % carbon fraction of membrane
%                                          % (membrane carbon from proteins + lipids)
    %CM = fProtM*CProt+fLipidM*CLipid;
    %CM = fProtM*CProt+fLipidM*(1-alphaPLip)*CLipid+fLipidM*(alphaPLip)*CPhospholipid;

    NE = NProt*fProtE+NRib*fRibE;             % Nitrogen fraction of biosynthetic components
    CE = CProt*fProtE+CRib*fRibE;             % Carbon fraction of biosynthesis pool
                                % || why not: CE = fRibE*CRibosomeEu + fProtE*CProt;
                                % what fraction of of biosynthetic apparatus is ribosome/protein?
    PE = fRibE*PRib;        % should this be: PE = PRib*fRibE; or: PE = fRibE*PRibosomeEu;

    NS = gammaDNA/gammaS*NDNA;   % N fraction of other structural components
    CS = (gammaDNA*CDNA+gammaLipid*CLipid)/gammaS; % Carbon fraction of other structural components
          % || why not: CS = (fLipidS*CPhosphoLipid + fDNAS*CDNA);  %%% Reason is that fLipid stands for other lipids that are not phospholipids. Here the idea is that the cell has some energy reserves
          % always as part of its structural pool
	PS = gammaDNA/gammaS*PDNA;

dNE_dfRibE = -NProt+NRib;
dCE_dfRibE = -CProt+CRib;
dPE_dfRibE = PRib;
% why is P fraction not defined for each of the Light, membrane, or other structure pools?

%% Solving for PLim
    aP= (3*DP.*P)*molarP/10^6 /(pDry*rho);
	% units = P:[mol/L]*[1000L/m^3]*DP:[m^2/hr]-> [mol/m/hr]*molarP:[g/mol] -> [g/m/hr]/[10^6um/m] -> [g/um/hr] % rho:[g/um^3] % alphaS:[um]
    % units of aP = [um^2/hr]
    coefA = aP.*CI.^2/alphaS^2-kST*PRib;	% why not kST*PE???
    coefB = -2*(1-gammaS)*CI.*aP/alphaS^2-gammaDNA*PDNA*kST;
    coefC = (1-gammaS)^2.*aP/alphaS^2;

    EPLim = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);

    rPLim = alphaS./(1-gammaS - CI.*EPLim); %1.0./( (1-gammaS)/alphaS-CI.*EPLim/alphaS);
    APLim = (1-gammaS -CI.*EPLim)./2; %should equal AMin;  APLim = alphaS/(2*rPLim)
    fProtAPLim = APLim.*0+1; % =1 (if purely P limited = invest as much as you can in uptake proteins)


%%% derivatives of EPlim w.r.t. parameters
	%%%% dE_dDIP
	daP_dDIP = 3*DP*molarP/10^6 /(pDry*rho);
	dcoefA_daP = CI.^2/alphaS^2;
	dcoefB_daP = -2*(1-gammaS)*CI/alphaS^2;
	dcoefC_daP = (1-gammaS)^2/alphaS^2;

	dEPLim_daP = ((2*coefA.*(-dcoefB_daP - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_daP-4*dcoefA_daP.*coefC-4*coefA.*dcoefC_daP)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_daP)./(4*coefA.^2) );
	dEPLim_dDIP = dEPLim_daP.*daP_dDIP;

	%%%% dE_dCI
    dcoefA_dCI = 2.*aP.*CI/alphaS^2;
	dcoefB_dCI = -2*(1-gammaS).*aP/alphaS^2;

	dEPLim_dCI = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );
	%dEPLim_dQ10Photo = dEPLim_dCI.*dCI_dQ10Photo;

	%%%% dE_dkST
	dcoefA_dkST = -PRib +dcoefA_dCI.*dCI_dkST;
	dcoefB_dkST = -gammaDNA*PDNA + dcoefB_dCI.*dCI_dkST;
	dEPLim_dkST = ((2*coefA.*(-dcoefB_dkST - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dkST-4*dcoefA_dkST.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dkST)./(4*coefA.^2) );

	%%%% dE_dalphaS
	dcoefA_dalphaS = -2*aP.*CI.^2/alphaS^3;
	dcoefB_dalphaS = -2*(-2*(1-gammaS)*CI.*aP)/alphaS^3;
	dcoefC_dalphaS = -2*(1-gammaS)^2.*aP/alphaS^3;

	dEPLim_dalphaS = ((2*coefA.*(-dcoefB_dalphaS - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dalphaS-4*dcoefA_dalphaS.*coefC-4*coefA.*dcoefC_dalphaS)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dalphaS)./(4*coefA.^2) );

%% Solving for NLim
    aN = 3*DN.*N*molarN/10^6 /(pDry*rho) ; %converting top to g/um/hr
    coefA = CI.^2.*aN/alphaS^2+kST.*CI*NM/2.0-kST*NL.*(CI-1)-kST*NE;
    coefB = -2*aN*(1-gammaS).*CI/alphaS^2-kST*gammaDNA*NDNA-kST*NM*(1-gammaS)/2.0-kST*AMin*NProt;
    coefC = aN.*(1-gammaS)^2/alphaS^2;

	ENLim = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);

    rNLim = alphaS./(1-gammaS - CI.*ENLim);
	ANLim = (1-gammaS-CI.*ENLim)./2; % mass fraction of periplasm in the cell. (= alphaS./(2*rNLim))
	fProtANLim = AMin./ANLim; % protein fraction of periplasm (maybe APLim./ANLim)

%%% derivatives of ENlim w.r.t. parameters
	%%%% dE_dDIP
	daN_dDIP = 0;
	dENLim_dDIP = aN.*0;

	%%%% dE_dCI
    dcoefA_dCI = 2*CI.*aN/alphaS^2 + kST*NM/2.0 -kST*NL;
    dcoefB_dCI = -2*aN*(1-gammaS)/alphaS^2;

    dENLim_dCI = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );
	%dENLim_dQ10Photo = dENLim_dCI.*dCI_dQ10Photo;

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

%% reset ANLim to AMin if the mass fraction of periplasm is too small (for large phytoplankton),
	rLG_ind = find(rNLim>rFullA);  % equivalent to find(ANLim<AMin) ;
        ANLim(rLG_ind) = AMin;
		fProtANLim(rLG_ind) = 1.0;

        coefA = CI.^2.*aN/alphaS^2+kST.*CI*(NM+NProt)/2.0-kST*NL.*(CI-1)-kST*NE;
        coefB = -2*aN*(1-gammaS).*CI/alphaS^2-kST*gammaDNA*NProt-kST*(NM+NProt)*(1-gammaS)/2.0;
        coefC = aN.*(1-gammaS)^2/alphaS^2;

        ENLim_LG = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);
        rNLim_LG = alphaS./(1-gammaS - CI.*ENLim_LG);
        ANLim_LG = (1-gammaS-CI.*ENLim_LG)./2; % mass fraction of periplasm in the cell. (= alphaS./(2*rNLim))

		ENLim(rLG_ind) =ENLim_LG(rLG_ind);
		rNLim(rLG_ind) =rNLim_LG(rLG_ind);
        ANLim(rLG_ind) = ANLim_LG(rLG_ind);

    %%% derivatives of ENlim for large phyto w.r.t. parameters
		%%%% dENLim_dDIP = aN.*0; %stays zero in both NLim cases
		%%%% dE_dCI
        dcoefA_dCI = 2*CI.*aN/alphaS^2 + kST*(NM+NProt)/2.0 -kST*NL;
        dcoefB_dCI = -2*aN*(1-gammaS)/alphaS^2;

        dENLim_dCI_LG = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );
        dENLim_dCI(rLG_ind) = dENLim_dCI_LG(rLG_ind);

		%dENLim_dQ10Photo_LG = dENLim_dCI.*dCI_dQ10Photo;
        %dENLim_dQ10Photo(rLG_ind) = dENLim_dQ10Photo_LG(rLG_ind);

		%%%% dE_dkST   % check this!
		dcoefA_dkST = dcoefA_dCI.*dCI_dkST +CI*(NM+NProt)/2.0-NL.*(CI-1)-NE;
		dcoefB_dkST = dcoefB_dCI.*dCI_dkST -gammaDNA*NProt-(NM+NProt)*(1-gammaS)/2.0;
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

%% Calculate growth rates (biosynthetic rate) under each condition
    muNLim = kST.*ENLim;  % units: [1/hr]
    muPLim = kST.*EPLim;  % units: [1/hr]

    muNLim_P= aP.*fProtANLim./(rNLim.^2.*(ENLim*PRib+gammaDNA*PDNA));
    muPLim_N= aN./(rPLim.^2.*((CI-1)*NL.*EPLim+EPLim*NProt+gammaDNA*NDNA+alphaS./(2.0.*rPLim)*NProt+alphaS./(2.0.*rPLim)*NM));
			% muNLim_P is the "P-limited" growth rate at the strategy where
			% N-limitation, Biosynthesis, and Photosynthesis balance.
			% (P-limited growth rate equation, evaluated at the N limited solution for E)
			% Likewise, muPLim_N is the "N-Limited" growth rate at the strategy where
			% P-limitation, biosynthesis, and photosynthesis balance.

			% We need to calculate these rates because we have both N and P limitation
			% in the model. The traits included in the model naturally allow for
			% trade-offs of nutrient limitation (smaller r) against synthesis or
			% carbon limitation (which require higher E and L). This allows us to
			% solve for a strategies which are colimited by either N, protein and
			% photosynthesis or P, protein and Photosynthesis

			% But it isn't guaranteed that we can find a trait combination that will
			% balance N and P limitation. So it will turn out that depending on the
			% conditions, the optimal strategy will either be N limited, P-limited, or
			% co-limited. The method for determining which starts by looking at how
			% much P the cell gets when it optimizes for N limitation, and how much N
			% it gets when it optimizes for P limitation.

			% If, for example, muNLim_P < kST*ENLim, then we conclude that the NLim
			% strategy is suboptimal (thus, not N-limited), and vice versa for the P-lim strategy.


%% initialize vector fields
all_nan = NaN(size(P));

%CP = all_nan;
%NP = all_nan;
%CN = all_nan;
%muOpt = all_nan;
%rOpt = all_nan;
AOpt = all_nan;
fProtAOpt = all_nan;
EOpt = all_nan;
LimType = all_nan;

%QP = all_nan;
%QC = all_nan;
%QN = all_nan;

%dC2P_dQ10Photo = all_nan;
%dE_dQ10Photo = all_nan;
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

%ibad=[];
%% loop through points to calculate cell quotas
EColim_handle = @EColim;
AColim_handle = @AColim;
%tic
parfor i =1:length(P)
%for i =1:length(P)
	EOpt_i = NaN;
	fProtAOpti = NaN;
	LimState = -1;
	CIi = NaN;
	kSTi = NaN;
if ~isnan(Irr(i)) & ~isnan(T(i)) & ~isnan(P(i)) & ~isnan(N(i))
	%disp(i)
    if muNLim(i)<muPLim(i)
        if muNLim(i)<muNLim_P(i) % N Limitation
            LimState = 0;
			CIi = CI(i);
			kSTi = kST(i);
            EOpt_i = ENLim(i);
            fProtAOpti = fProtANLim(i); % = 2*AMin./(1-gammaS-CI.*ENLim)
            %AOpti= ANLim(i);
            %rOpt_i = rNLim(i);

            %%% derivatives
            dE_dCI(i) = dENLim_dCI(i);
			dE_dalphaS(i) = dENLim_dalphaS(i);
			dE_dfRibE(i) = dENLim_dfRibE(i);
			dE_dkST(i) = dENLim_dkST(i);
			dE_dDIP(i) = 0; % dENLim_dDIP(i) = 0

			if ANLim(i)>AMin
				dfProtAOpt_dE(i) = fProtAOpti^2*CIi/(2*AMin); %if ANLim> AMin
				dfProtAOpt_dCI(i) = fProtAOpti^2/(2*AMin)*(EOpt_i);
                %dfProtAOpt_dCI(i) = -2*AMin./(1-gammaS -CIi.*EOpt_i).^2.*(-1.*EOpt_i); %same as before
                %dfProtAOpt_dE(i) = -2*AMin./(1-gammaS -CIi.*EOpt_i).^2.*(-1.*CIi);
                pdfProtAOpt_fRibE(i) = 0;
				dfProtAOpt_dkST(i) = dfProtAOpt_dE(i)*dE_dkST(i)+dfProtAOpt_dCI(i)*dCI_dkST(i);
				dfProtAOpt_dDIP(i) = dfProtAOpt_dE(i)*dE_dDIP(i);
			else
				dfProtAOpt_dE(i) = 0;
				dfProtAOpt_dCI(i) = 0;
                pdfProtAOpt_fRibE(i) = 0;
				dfProtAOpt_dkST(i) = 0;
				dfProtAOpt_dDIP(i) = 0;
            end

            %dE_dQ10Photo(i) = dENLim_dQ10Photo(i);

        else
            LimState = 2;
			%CIi = CI(i);
			%kSTi = kST(i);
			%aNi = aN(i);
			%aPi = aP(i);

            EMin=1e-4;
            EMax=(1-gammaS)/CI(i);
			%[EOpt_i, dE_dCI(i), dE_dalphaS(i), dE_dfRibE(i),pdE_kST,dE_daP] = EColim(EMin,EMax);
			[EOpt_i, dE_dCI(i), dE_dalphaS(i), dE_dfRibE(i),pdE_kST,dE_daP] = feval(EColim_handle,EMin,EMax,i); %if using parfor loop
			dE_dkST(i) = pdE_kST + dE_dCI(i).*dCI_dkST(i);
			dE_dDIP(i) = dE_daP*daP_dDIP(i);

            %fProtAOpti = AMin*2*rOpt_i/alphaS;  %  =2*AMin/(1-gammaS -CIi.*EOpt_i) this is only for NLim
            %case, when minimum P uptake is used

			%[fProtAOpti, dfProtAOpt_dCI(i), dfProtAOpt_dE(i), pdfProtAOpt_fRibE(i),pdfProtAOpt_DIP] = AColim(EOpt_i);
			[fProtAOpti, dfProtAOpt_dCI(i), dfProtAOpt_dE(i), pdfProtAOpt_fRibE(i),pdfProtAOpt_DIP] = feval(AColim_handle,EOpt_i,i);
			dfProtAOpt_dkST(i) = dfProtAOpt_dE(i)*dE_dkST(i)+dfProtAOpt_dCI(i)*dCI_dkST(i);
			dfProtAOpt_dDIP(i) = pdfProtAOpt_DIP + dfProtAOpt_dE(i)*dE_dDIP(i);

            %dE_dQ10Photo(i) = dE_dCI(i)*dCI_dQ10Photo(i);
            %muPColim = fProtAOpt_i.*aP/rOpt_i.^2/((EOpt*PE+(gammaDNA)*PDNA ));
            %muNColim = aNi./rOpt_i.^2/(EOpt*NE+EOpt.*(CIi-1)*NL+gammaS*NS+fProtAOpt_i*alphaS./(2.0*rOpt_i)*NProt+alphaS/(2.0*rOpt_i)*NM);
        end
    elseif muPLim(i) <= muNLim(i)
        if muPLim(i)<muPLim_N(i) % P Limitation
            LimState = 1;
            EOpt_i = EPLim(i);
			fProtAOpti = fProtAPLim(i); % = 1
            %AOpti= APLim(i);
            %rOpt_i = rPLim(i);  % = alphaS./(1-gammaS - CI.*EPLim);
			CIi = CI(i);
			kSTi = kST(i);

			%%% derivatives
            dE_dCI(i) = dEPLim_dCI(i);
			dE_dalphaS(i) = dEPLim_dalphaS(i);
            dfProtAOpt_dE(i) = 0;
			dfProtAOpt_dCI(i) = 0;
			dE_dfRibE(i) = 0;
            pdfProtAOpt_fRibE(i) = 0;
			%dE_dQ10Photo(i) = dEPLim_dQ10Photo(i);
			dE_dkST(i) = dEPLim_dkST(i);
			dfProtAOpt_dkST(i) = 0;
			dE_dDIP(i) = dEPLim_dDIP(i);
			dfProtAOpt_dDIP(i) = 0;

        else
            LimState = 3; %same as 2. changed for debugging
			%CIi = CI(i);
			%kSTi = kST(i);
			%aNi = aN(i);
			%aPi = aP(i);

            EMin = 1e-4;
            EMax = (1-gammaS)/CI(i);
			%[EOpt_i, dE_dCI(i), dE_dalphaS(i),dE_dfRibE(i),pdE_kST,dE_daP] = EColim(EMin,EMax); %make i an input and change func to use indices
			[EOpt_i, dE_dCI(i), dE_dalphaS(i), dE_dfRibE(i),pdE_kST,dE_daP] = feval(EColim_handle,EMin,EMax,i);
			dE_dkST(i) = pdE_kST + dE_dCI(i).*dCI_dkST(i);
			%dE_dfRibE(i)
			dE_dDIP(i) = dE_daP*daP_dDIP(i);

            %[fProtAOpti, dfProtAOpt_dCI(i), dfProtAOpt_dE(i), pdfProtAOpt_fRibE(i),pdfProtAOpt_DIP] = AColim(EOpt_i); %AColim(EOpt_i);
                %AColim(E) must be between 0 and 1
			[fProtAOpti, dfProtAOpt_dCI(i), dfProtAOpt_dE(i), pdfProtAOpt_fRibE(i),pdfProtAOpt_DIP] = feval(AColim_handle,EOpt_i,i);
			dfProtAOpt_dkST(i) = dfProtAOpt_dE(i)*dE_dkST(i)+dfProtAOpt_dCI(i)*dCI_dkST(i);
			dfProtAOpt_dDIP(i) = pdfProtAOpt_DIP + dfProtAOpt_dE(i)*dE_dDIP(i);


            %muPColim = fProtAOpti.*aPi/rOpt_i^2/((EOpt_i*PE+(gammaDNA)*PDNA ));
            %muNColim = aN./rOpt_i.^2/(EOpt_i*NE+EOpt_i.*(CIi-1)*NL+gammaS*NS+fProtAOpti*alphaS./(2.0*rOpt_i)*NProt+alphaS/(2.0*rOpt_i)*NM);
        end
    end % end cases

%%% save  cell quotas of P, N, and C
            %%% note: AOpti = alphaS./(2.0*rOpt_i) = MOpt

    %        % AProt = fProtAOpt * AOpt;
    %        if fProtAOpt(i)==1
    %            AProt = AOpti; % = 0.5*(1-gammaS-CIi.*EOpt_i); % = AOpti
    %        else %(fProtAOpt = AMin*2*rOpt_i/alphaS)
    %            AProt = AMin;
                % NLim always has minimum allocation to Puptake proteins
                % currently same for Colim. How to  allocate more to uptake
                % proteins? AColim function give stange values
    %        end
    %%%% replace in QN & QC eq: fProtAOpti*alphaS./(2.0*rOpt_i) = fProtAOpti*AOpti = AMin
		%QP(i) =((EOpt_i*fRibE*PRib+(gammaDNA)*PDNA ))/molarP;
		%QN(i) = (EOpt_i*NE+EOpt_i.*(CIi-1)*NL+gammaS*NS+fProtAOpti*AOpti*NProt+alphaS/(2.0*rOpt_i)*NM)/molarN;
		%QC(i) = (EOpt_i*CE+EOpt_i.*(CIi-1)*CL+gammaS*CS+fProtAOpti*AOpti*CProt+alphaS/(2.0*rOpt_i)*CM)/molarC .*(1+24.0*.25*kSTi.*EOpt_i*(1+PhiS));

		EOpt(i) = EOpt_i;
        %fprintf('EOpt(i) is: % 3.4e \n', EOpt(i));
        fProtAOpt(i) = fProtAOpti;
        %rOpt(i) = alphaS./(1-gammaS - CIi.*EOpt_i); %= rOpt_i
		%muOpt(i) = kSTi.*EOpt_i;
        LimType(i) = LimState;
        %dE_dQ10Photo(i) = dE_dCI(i)*dCI_dQ10Photo(i);
        %%dE_dQ10Photo_dQ10Photo = dE_dCI(i)*dCI_dQ10Photo_dQ10Photo(i) + dCI_dQ10Photo(i)* dEColim_dCI_dQ10Photo


% This warns that allocations are unphysical. will not work when using a parfor loop
%            if fProtAOpt(i) > 1 | fProtAOpt(i) < 0
%				ibad = [ibad, i ];
%				% fProtAOpt(i) =1;
%				if length(ibad) < 4
%                    fprintf('at i = %i , fProtAOpt  = %0.5g \n',i,fProtAOpt(i))
%
%					%keyboard
%				end
%            end

else
    fprintf('NaN value in cell model i = %d : Irr(i) = %0.5g , T(i) = %0.5g, P(i) = %0.5g, N(i) = %0.5g', i,Irr(i),T(i),P(i),N(i))

end % if ~isnan(...)
end % end for loop
%toc

AOpt = (1-gammaS-CI.*EOpt)./2;  % periplasm fraction of cell
MOpt = AOpt;                    % membrane fraction of cell
rOpt = alphaS./(1-gammaS - CI.*EOpt); % cell radius
muOpt = kST.*EOpt;

QPnostor = (EOpt.*PE + gammaDNA*PDNA )./molarP;
%C2P = (molarP/molarC)*((EOpt*CE+EOpt.*(CI-1)*CL+gammaS*CS+(AMin./AOpt).*(1-gammaS - CI.*EOpt)/2.0*CProt+ (1-gammaS - CI.*EOpt)/2.0*CM).*(1+24.0*.25*kST.*EOpt*(1+PhiS)) ) ./ (EOpt*PE+(gammaDNA)*PDNA); % only tiny random error (10^-13) from CP

%%%%%%%%% Add P storage and phospholipids %%%%%%%%%%%%%%%%%%%
 % each scaled by logistic function with 2 parameters each: PLip: PLip_scale & PLip_PCutoff; PStor: PStor_scale & PStor_rCutoff
 % can turn of PLip and PStorage by setting alphaPLip = 0 and fStorage = 0;
 PLip = alphaPLip*MOpt*PPhospholipid .*(1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))));
 PStor = fStorage.*5000.*P .*(1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff)))); % change max value? change units of P (was P*5000)

 % should PStor only exist if nitrogen limited?
%       Nindx = find(LimType ==0);
%       PStor = zeros(size(P));
%       PStor(Nindx) = fStorage .*P .*(1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))));

%%%%%%%% Calculate Cell Quotas %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% what is the extra part added on to QC? %L=kST.*E*(1+PhiS)/alphaI  % CI  = 1+(kST.*(1+PhiS))./alphaI;
	% QC = QC*(1 +L*alphaI*0.25*24)
	% alphaI = alphaPhoto*Q10Photo.^((T-T0)./10) = F1*(temperature dependence of carbon fixation)
	% L*F1  = total allocation to carbon fixation component of photosynthetic apparatus
	% L*alphaI*0.25*24 = % of cell dedicated to photosynthetic carbon fixation * temperature dependence * rate constant of 0.25 hr^-1 * 24 hr/day
	% = amount of carbon fixed by photosynthesis in a day ? is this equivalent to the carbohydrate quota of the cell?
%	% totalC = ((.3*alphaS/r+.05)*CPhospholipid+.7*alphaS*CProtein/r+.10*CLipid+.04*CCarb+(1-alphaE)*E*CProtein+alphaE*E*CRibosomeEu+(1-E-alphaS/r-gamma)*CProtein+.01*CDNA)/12.0;

% Cell quotas in moles P, C, and N
QP =((EOpt.*PE +(gammaDNA)*PDNA + PLip +PStor ))/molarP;
QC = (EOpt.*CE +EOpt.*(CI-1)*CL +gammaS*CS +fProtAOpt.*AOpt.*CProt +MOpt.*CM)./molarC .*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
QN = (EOpt.*NE +EOpt.*(CI-1).*NL +gammaS*NS +fProtAOpt.*AOpt.*NProt +MOpt.*NM)./molarN;


CP = QC./QP;
NP = QN./QP;
CN = QC./QN;

% ****** Calculate derivatives of C:P ****** %
 %if (par.optim == on)

		% QC = QC(E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kSTi,PhiS)
		% Take partial derivatives with respect to all these variables
		%pdQC_E = (CE + (CI-1)*CL).*(1+24.0*.25*kST.*EOpt.*(1+PhiS))./molarC + (EOpt.*CE+EOpt.*(CI-1).*CL+gammaS*CS+fProtAOpt.*alphaS./(2.0*rOpt).*CProt+alphaS./(2.0*rOpt)*CM).*24.0*0.25.*kST*(1+PhiS)/molarC;
		pdQC_E = (CE + (CI-1)*CL).*(1+24.0*.25*kST.*EOpt.*(1+PhiS))./molarC + (EOpt.*CE +EOpt.*(CI-1)*CL +gammaS*CS +fProtAOpt.*AOpt.*CProt +MOpt.*CM).*24.0*0.25.*kST*(1+PhiS)/molarC;
		pdQC_CE = EOpt./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
		pdQC_CI = (EOpt*CL).*(1+24.0*.25*kST.*EOpt*(1+PhiS))/molarC;
		pdQC_CL = 0; % do later
		pdQC_gammaS = CS./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS)); % do later
		pdQC_CS = 0; % do later
        pdQC_fProtAOpt = AOpt.*CProt./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
        pdQC_AOpt = fProtAOpt.*CProt./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
        pdQC_MOpt =  CM./molarC.*(1+24.0*.25*kST.*EOpt.*(1+PhiS));
		%pdQC_alphaS = 0;
		%pdQC_rOpt = (-fProtAOpt*alphaS./(2.0*rOpt.^2)*CProt-alphaS./(2.0*rOpt.^2)*CM).*(1+24.0*.25*kST.*EOpt.*(1+PhiS))/molarC;
		pdQC_rOpt = 0;
        pdQC_CProt = 0; % do later
		pdQC_CM = 0; % do later
		pdQC_kST = 24*0.25.*EOpt.*(1+PhiS).*(EOpt.*CE +EOpt.*(CI-1)*CL +gammaS*CS +fProtAOpt.*AOpt.*CProt +MOpt.*CM)./molarC ;
		pdQC_PhiS = 24*0.25.*kST.*EOpt.*(EOpt.*CE +EOpt.*(CI-1)*CL +gammaS*CS +fProtAOpt.*AOpt.*CProt +MOpt.*CM)./molarC ;
		% QP = QP(E,fRibE,PRib,gammaDNA,PDNA)
		% Take partial derivatives with respect to all these variables
		pdQP_E = PE/molarP;
		pdQP_PE = EOpt./molarP;
		pdQP_gammaDNA = PDNA/molarP;
		pdQP_PDNA = gammaDNA/molarP;
		% expand PLip to treat QP as if PLip eqn is written out fully, so marams in PLIp are directly in QP
		pdQP_PLip = 1./molarP;
        pdQP_MOpt = pdQP_PLip*(PLip./MOpt);
        pdQP_alphaPLip = pdQP_PLip*(PLip./alphaPLip);
        pdQP_PLipscale = (pdQP_PLip)*alphaPLip*MOpt*PPhospholipid .* (-1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))).^2) .*exp(-PLip_scale*(P-PLip_PCutoff)) .* -1.*(P-PLip_PCutoff);
        pdQP_PLipPCutoff = (pdQP_PLip)*alphaPLip*MOpt*PPhospholipid .* (-1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))).^2) .*exp(-PLip_scale*(P-PLip_PCutoff)) .*PLip_scale;
		% expand PStor
        pdQP_PStor = 1./molarP;
        pdQP_fStorage = pdQP_PStor*(PStor./fStorage);
        pdQP_rOpt = pdQP_PStor * fStorage*5000.*P .*(-1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))).^2) .* exp(-PStor_scale.*(rOpt-PStor_rCutoff)) .*-PStor_scale;
        pdQP_rCutoff = (pdQP_PStor)*fStorage*5000.*P .*(-1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))).^2) .* exp(-PStor_scale.*(rOpt-PStor_rCutoff)) .*PStor_scale;
        pdQP_PStorscale = (pdQP_PStor) * -1*PStor.^2./(fStorage*5000.*P) .*exp(-PStor_scale.*(rOpt-PStor_rCutoff)) .* -1.*(rOpt-PStor_rCutoff);
		%DIP is in both PLip and PStor; easier to follow if seperate then add together
		pdPLip_DIP = alphaPLip*MOpt*PPhospholipid .*(-1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))).^2).*(-PLip_scale*exp(-PLip_scale*(P-PLip_PCutoff))) ;
		pdPStor_DIP = fStorage.*5000.*(1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))));
		pdQP_DIP = pdQP_PLip*pdPLip_DIP + pdQP_PStor*pdPStor_DIP;


		% C:P = C:P(E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kSTi,PhiS,fRibE,PRib,gammaDNA,PDNA)
		% Use the partial derivatives of QC and QP to calculate the partial derivatives of C:P. pdC2P_E is the only complicated ones
		pdC2P_E = (QP.*pdQC_E - QC.*pdQP_E)./QP.^2;
		pdC2P_CE = pdQC_CE./QP;
		pdC2P_CI = pdQC_CI./QP;
		pdC2P_CL = pdQC_CL./QP;
		pdC2P_gammaS = pdQC_gammaS./QP;
		pdC2P_CS = pdQC_CS./QP;
		pdC2P_fProtAOpt = pdQC_fProtAOpt./QP;
        pdC2P_AOpt = pdQC_AOpt./QP;
        pdC2P_MOpt = (QP.*pdQC_MOpt - QC.*pdQP_MOpt)./QP.^2;
        pdC2P_rOpt = -pdQP_rOpt.*QC./QP.^2;
		pdC2P_CProt = pdQC_CProt./QP;
		pdC2P_CM = pdQC_CM./QP;
		pdC2P_kST = pdQC_kST./QP;
		pdC2P_PhiS = pdQC_PhiS./QP;
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



		% Now we need to know the derivatives of E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kSTi,PhiS,fRibE,PRib,gammaDNA,PDNA with respect to the variables we are changing
        %dE_dCI, dfProtAOpt_dE, dfProtAOpt_dCI

		pdPLip_DIP = alphaPLip*MOpt*PPhospholipid .*(-1./(1 + exp(-PLip_scale*(P-PLip_PCutoff))).^2).*(-PLip_scale*exp(-PLip_scale*(P-PLip_PCutoff))) ;
		pdPStor_DIP = fStorage.*5000.*(1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))));
		pdQP_DIP = pdQP_PLip*pdPLip_DIP + pdQP_PStor*pdPStor_DIP;

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
        %d2C2P_dQ10Photo2 =pdC2P_E.*dE_dQ10Photo_dQ10Photo + dE_dQ10Photo.*dC2P_dE_dQ10Photo + pdC2P_CI.*dCI_dQ10Photo_dQ10Photo +dC2P_dCI_dQ10Photo.*dCI_dQ10Photo;

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
		% for points where DIP was reset from a negative value, set the derivative to 0; done in eqC2Puptake now
		% dC2P_dDIP(negDIPindx) = 0;

%% add PStorage
%%% Why is Pstorage tacked on at the end? %%%
%%% dont know how to link P storage to fitness
% could be reason to have max for PStorage - places with a lot of PStorage are often places with Iron limitation (but no iron limitation in model)
% if P storage is taken into account, model will think cell is p limited.
% E.G. southern ocean: P and N repllete, but P is drawn down a little more

% alphaPLip = 0.12 (fraction of membrane that is phopholipids)
% lPCutoff is implemented because one of the main ways that small cells can
			% vary their P-quota is by substituting Phospholipids for Sulfolipids at low
			% P. This only really matters for the small cells
% CStor s
	    % model treats storage differently for large cells
	    % compared to small cells. In particular, large cells, with radius above
	    % r0Cutoff (which seems not to be in this file anymore- we can discuss
	    % whether it is needed or not), have luxury storage while cells below do
	    % not. In order to implement this, I multiplied the storage term (the part
	    % proportional to environmental P) by a function that goes to 0 for small
	    % radius and goes to 1 for large radius. This transition happens around the
	    % r=r0cutoff, and CStor is the "sharpness" of the transition region.  See
	    % the expression: PStorage =
	    % exp(lfStorage)*5000*P/(1+exp(-CStor*(rOpt-r0Cutoff)))+PLip

%
%         %PLip = alphaPLip*.5./rOpt*PPhospholipid/CPhospholipid*molarC/molarP./(1+exp((P./lPCutoff - 1)))   ; % original, doesn;t make sense
%         PLip = alphaPLip./(2*rOpt)*PPhospholipid/CPhospholipid*molarC/molarP./(1+exp(-1*(P./exp(lPCutoff) - 1)))   ; %low P ->low pLip
%         %PLip = alphaPLip*.5./rOpt*PPhospholipid/CPhospholipid*molarC/molarP./(1+exp(-1*(log(P)-lPCutoff)))   ; %at low P, very little PLip
%         % log(P)-lPCutoff = log(P/PCutoff)
%
%         % PLip = how much P from phospholipids is in cell
%         % A = alphaS/(2*rOpt) ; fProtAPLim = AMin./APLim; %A = fraction of cell dry mass that is nutrient uptake proteins
%         % alphaS/(2*rOpt)=(1-gammaS - CIi.*EOpt_i)/2
%         % PLip = alphaPLip/(2*rOpt) *... % PLip =fraction of cell dry mass that is phospholipids
%         %      =alphaPLip/alphaS * (1-gammaS - CIi.*EOpt_i)/2
%         %     alphaPLip = alphaS*fPLipM
%         %     fPLipM = alphaPLip/alphaS
%
%         PStorage = fStorage*5000*P./(1+exp(-CStor.*(rOpt-r0Cutoff)))+PLip;% units [P/C] % at small r, low PStorage
%
%         CPwithPStor = 1./(1./CP+PStorage);
%         NPwithPStor = 1./(1./NP + PStorage.*CN);
%
%         %MOpt = alphaS./(2.0.*rOpt); %membrane fraction of cell
%         fPlipM = alphaPLip./(1+exp(-1*(P./exp(lPCutoff) - 1)));
%         PStor = fStorage*5000*P./(1+exp(-CStor.*(rOpt-r0Cutoff)));
%         %fPLipCell = fPlipM * alphaS./(2.0*rOpt);  % alphaS./(2.0*rOpt)= membrane fraction of cell; fPlipM = phospholipid fraction of membrane
%         PM=fPlipM*PPhospholipid;
%         QPstor =((EOpt.*PE+(gammaDNA)*PDNA + alphaS./(2.0*rOpt).*PM +PStor ))/molarP;
%         QCstor = (EOpt.*CE+EOpt.*(CI-1)*CL+gammaS*CS+fProtAOpt.*alphaS./(2.0*rOpt)*CProt+alphaS./(2.0*rOpt)*CM)/molarC .*(1+24.0*.25*kSTi.*EOpt_i*(1+PhiS));
%
%         CPstor = QCstor./QPstor;



       % what is the max PStorage value

        %dC2P_dfStorage = QCstor./(QPstor.^2).*(1./molarP).*5000.*P .*(1./(1 + exp(-PStor_scale.*(rOpt-PStor_rCutoff))));
%keyboard;

%% Calculate C:P derivatives including PStorage %%%
		%dPLip_drOpt = -PLip./rOpt;
		%dPStorage_drOpt = -(PStorage-PLip)./(1+exp(-CStor.*(rOpt-r0Cutoff))).*(-CStor.*exp(-CStor.*(rOpt-r0Cutoff))) +dPLip_drOpt;
		% dQ10Photo
		%dPStorage_dQ10Photo = dPStorage_drOpt.* drOpt_dQ10Photo;
		%dC2Pstor_dQ10Photo = -CPwithPStor.^2 .*(-1./CP.^2 .*dC2P_dQ10Photo +dPStorage_dQ10Photo);
		% dfStorage
		%dPStorage_dfStorage = 5000*P./(1+exp(-CStor.*(rOpt-r0Cutoff)));
		%dC2Pstor_dfStorage = -CPwithPStor.^2 .*dPStorage_dfStorage;

        %dlPCutoff
        %dPLip_dlPCutoff = -1*alphaPLip./(2*rOpt)*PPhospholipid/CPhospholipid*molarC/molarP./(1+exp(-1*(P./exp(lPCutoff) - 1))).^2 .*exp(-1*(P./exp(lPCutoff) - 1)) .* (-1*-1*P./exp(lPCutoff)^2).* exp(lPCutoff);


%%% change Nans to zeros (trying to fix error in newton solver)
%CP(isnan(CP)) = 0;

%% Output
    out.CP = CP;
    out.CPnostor = QC./QPnostor;
    out.NP = NP;
    out.CN = CN;
    out.LimType = LimType;
    out.r = rOpt;
    out.E = EOpt;
	out.A = AOpt;
    out.mu = muOpt;
	out.QP = QP;
	out.PStor = PStor;
    out.PLip = PLip;
	out.CI = CI;
	out.L = CI.*EOpt - EOpt;

    % name a second structure for derivatives C2Px?
	out.dC2P_dQ10Photo = dC2P_dQ10Photo;
	%out.dC2P_dQ10Photo = dC2Pstor_dQ10Photo;
	%out.d2C2P_dQ10Photo2;
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

	out.dE_dDIP = dE_dDIP;
    %out.dE_dCI = dE_dCI;
    %out.dE_dQ10Photo = dE_dQ10Photo;
    %out.dCI_dQ10Photo =dCI_dQ10Photo;
    %out.dfProtAOpt_dE = dfProtAOpt_dE;
    %out.dfProtAOpt_dCI = dfProtAOpt_dCI;
    %out.dCI_dQ10Photo = dCI_dQ10Photo;
    %out.fProtA = fProtAOpt;
    %out.dfProtA_dQ10Photo  = dfProtAOpt_dQ10Photo;

	disp('Cell Model Done')
    % L, A, S?
    %save par values

%% define functions
	function [EColim, dEColim_dCI, dEColim_dalphaS, dEColim_dfRibE,pdEColim_kST,dEColim_daP] = EColim(EMin,EMax,i) % rhsFuncColim(EMin,EMax)
		CIi = CI(i);
		kSTi = kST(i);
		aNi = aN(i);
		aPi = aP(i);
        % polynomial coefficients
		c3 = kSTi^2*NProt*PE + (aPi*aNi/alphaS^4)*CIi^3 - (aPi/alphaS^2)*kSTi*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		% function to return all derivs of c3 wrt params
		c2 = kSTi^2*NProt*gammaS*PS - 3*(aPi*aNi/alphaS^4)*CIi^2*(1-gammaS) + (aPi/alphaS^2)*kSTi*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) - (aPi/alphaS^2)*kSTi*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		c1 = 3*(aPi*aNi/alphaS^4)*CIi*(1-gammaS)^2 + (aPi/alphaS^2)*kSTi*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		c0 = -(aPi*aNi/alphaS^4)*(1-gammaS)^3;
		%rhsFuncColim = c3*E.^3 + c2*E.^2 + c1*E +c0;

		E0 = (EMin+EMax)/2; %starting point for solver

		if ~all(isfinite([c3,c2,c1,c0,E0]))
			fprintf('c3 = %3.5f \n', c3)
			fprintf('c2 = %3.5f \n', c2)
			fprintf('c1 = %3.5f \n', c1)
			fprintf('c0 = %3.5f \n', c0)
			fprintf('E0 = %3.5f \n', E0)
		end
%{
		% solve for the zeros of rhsFuncColim
        rhsFuncColim = @(E) (c3*E^3 + c2*E^2 + c1*E +c0);

        f0 = rhsFuncColim(E0); %initial function value scale by function value to avoid tiny
        %f = @(E) rhsFuncColim(E)/f0; % scale to avoid errors in solver caused by tiny function values
		% f = @(E) rhsFuncColim(E)/1e-7; % scale to avoid errors in solver caused by tiny function values
		%print out f(E) for sensible values of E (make sure 1e-7 works for whole range)
        function [out, J] = f(E)
            out = (c3*E^3 + c2*E^2 + c1*E +c0)/1e-7;
            if nargout >1
               J = (3*c3*E^2 + 2*c2*E + c1)/1e-7;
            end
        end
		sol = fsolve_cmplx(@f,E0);
%}

        sol = complex_cubic_zero(c3,c2,c1,c0,E0);
		EColim = sol;

        %fprintf('EColim is: % 3.4e \n', EColim);
		%fprintf('f0 is: % 3.4e \n', f0);
		%fprintf('E0 is: % 3.4e \n', E0);
		%fprintf('real(EColim) is: % 3.4e \n', real(EColim));
		%fprintf('imag(EColim) is: % 3.4e \n', imag(EColim));
		%%% dE_dDIP
		%dE_daP
		dc3_daP = (aNi/alphaS^4)*CIi^3 - (1/alphaS^2)*kSTi*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		dc2_daP = - 3*(aNi/alphaS^4)*CIi^2*(1-gammaS) + (1/alphaS^2)*kSTi*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) - (1/alphaS^2)*kSTi*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc1_daP = 3*(aNi/alphaS^4)*CIi*(1-gammaS)^2 + (1/alphaS^2)*kSTi*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc0_daP = -(aNi/alphaS^4)*(1-gammaS)^3;
		dEColim_daP = -(EColim^3*dc3_daP+EColim^2*dc2_daP+EColim*dc1_daP+dc0_daP)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		% double checking everything
		%fprintf('c3 is: % 3.2e \n', dc3_daP/(imag(c3)/eps^3));
		%fprintf('c2 is: % 3.2e \n', dc2_daP/(imag(c2)/eps^3));
		%fprintf('c1 is: % 3.2e \n', dc1_daP/(imag(c1)/eps^3));
		%fprintf('c0 is: % 3.2e \n', dc0_daP/(imag(c0)/eps^3));
		%fprintf('dEColim is: % 3.2e \n', dEColim_daP/(imag(EColim)/eps^3)); % dE_daP/dE_DIP = 1/daP_dDIP

		%%% dE_dCI
		dc3_dCI = 3*(aPi*aNi/alphaS^4)*CIi^2 - (aPi/alphaS^2)*kSTi*(NE+(CIi-1)*NL-0.5*CIi*NM) - (aPi/alphaS^2)*kSTi*CIi*(NL-0.5*NM);
		dc2_dCI = -6*(aPi*aNi/alphaS^4)*CIi*(1-gammaS) + (aPi/alphaS^2)*kSTi*(1-gammaS)*(NL-0.5*NM) - (aPi/alphaS^2)*kSTi*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc1_dCI = 3*(aPi*aNi/alphaS^4)*(1-gammaS)^2;
		dc0_dCI = 0;
		dEColim_dCI = -(EColim^3*dc3_dCI+EColim^2*dc2_dCI+EColim*dc1_dCI+dc0_dCI)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		%dE_dkST   %check this
		% dc3_dkST = 2*kSTi*NProt*PE - (aPi/alphaS^2)*CIi*(NE+(CIi-1)*NL + dc3_dCI*dCI_dkST(i);
		% dc2_dkST = 2*kSTi*NProt*gammaS*PS + (aPi/alphaS^2)*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) -(aPi/alphaS^2)*CIi*(gammaS*NS+0.5*(1-gammaS)*NM) +dc2_dCI*dCI_dkST(i);
		% dc1_dkST = (aPi/alphaS^2)*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM) +dc1_dCI*dCI_dkST(i);
		% dc0_dkST = 0;
		% dEColim_dkST = -(EColim^3*dc3_dkST+EColim^2*dc2_dkST+EColim*dc1_dkST+dc0_dkST)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		% pdE_kST   %check this
		pdc3_kST = 2*kSTi*NProt*PE - (aPi/alphaS^2)*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		pdc2_kST = 2*kSTi*NProt*gammaS*PS + (aPi/alphaS^2)*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) -(aPi/alphaS^2)*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		pdc1_kST = (aPi/alphaS^2)*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		pdc0_kST = 0;
		pdEColim_kST = -(EColim^3*pdc3_kST+EColim^2*pdc2_kST+EColim*pdc1_kST+pdc0_kST)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		%dalphaS
		dc3_dalphaS = -4*(aPi*aNi/alphaS^5)*CIi^3 - -2*(aPi/alphaS^3)*kSTi*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		dc2_dalphaS = -4*(-3*aPi*aNi/alphaS^5)*CIi^2*(1-gammaS) + -2*(aPi/alphaS^3)*kSTi*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) - -2*(aPi/alphaS^3)*kSTi*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc1_dalphaS = -4*3*(aPi*aNi/alphaS^5)*CIi*(1-gammaS)^2 + -2*(aPi/alphaS^3)*kSTi*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc0_dalphaS = 4*(aPi*aNi/alphaS^5)*(1-gammaS)^3;
		dEColim_dalphaS = -(EColim^3*dc3_dalphaS+EColim^2*dc2_dalphaS+EColim*dc1_dalphaS+dc0_dalphaS)/(3*EColim^2*c3 +2*EColim*c2 +c1);

		%dfRibE
		dc3_dPE = kSTi^2*NProt;
		dEColim_dPE = -(EColim^3*dc3_dPE)/(3*EColim^2*c3 +2*EColim*c2 +c1);
		dc3_dNE = -(aPi/alphaS^2)*kSTi*CIi;
		dc2_dNE = (aPi/alphaS^2)*kSTi*(1-gammaS);
		dc1_dNE = 0;
		dc0_dNE = 0;
		dEColim_dNE = -(EColim^3*dc3_dNE+EColim^2*dc2_dNE)/(3*EColim^2*c3 +2*EColim*c2 +c1);
		dEColim_dfRibE = dEColim_dPE.*dPE_dfRibE + dEColim_dNE.*dNE_dfRibE;
    end

	function dEColim_dCI_dQ10Photo = dEColim_dCI_dQ10Photo(E)
		dEColim_dCI_dQ10Photo=dEColim_dCI_dCI*dCI_dQ10Photo +dEColim_dCI_dE*dE_dQ10Photo
	end

%dEdalphaS_colim

    function [fProtAColim, dfProtA_dCI, dfProtA_dE, dfProtA_dfRibE,pdfProtA_DIP]  =AColim(EOpt_i,i)
		CIi = CI(i);
		kSTi = kST(i);
		aNi = aN(i);
		aPi = aP(i);
        % fProtA (fraction periplasm (A) devoted to protein)
        % uptake rate of phosphorous = linear function of how much prot in periplasm
        %AColim = (-2*aN*(-1+CI*E+gammaS)^2+alphaS^2*E*kST*(2*E*(NE+(CI-1)*NL)+NM*(1-CI*E-gammaS)+2*gammaDNA*NProt))/(alphaS^2*E*(CI*E-1+gammaS)*kST*NProt);

        %fProtAPlim =kSTi*EOpt_i*(EOpt_i*PE+gammaS*PS)/(aPi/alphaS^2*(1-gammaS-CIi*EOpt_i)^2)
        %fProtAPlim = kSTi*EPLim(i)*(EPLim(i)*PE+gammaS*PS)/(aPi*2*APLim(i)^2/alphaS^2)
        %fProtANlim = 2*aNi/alphaS^2*(1-gammaS-CIi*ENLim(i))/(kSTi*ENLim(i)) -2*(ENLim(i)*NE+(CIi-1)*ENLim(i)*NL+gammaS*NS)/(NProt*(1-gammaS-CIi*ENLim(i))) - NM/NProt

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
            %else
			%	fProtAColim = 1;
			%	fprintf('Error in AColim at i=%i ,fProtAColim=%.4g \n',i,fProtAColim)
            end
        end

        %%% derivatives of AColim w.r.t. parameters
		%pdfProtA_dDIP
		pdcoefA_DIP = 0;
		pdcoefB_DIP = 0;
		%pdcoefC_DIP = -DN(i)/DP(i)*molarN/molarP*N(i)*(EOpt_i*PE+gammaS*PS)/(P(i)^2*-1);
		pdcoefC_DIP = coefC/(-P(i));
		pdfProtA_DIP = ( (2*coefA)*(-pdcoefB_DIP+0.5/sqrt(coefB^2-4*coefA*coefC)*(2*coefB*pdcoefB_DIP-4*coefA*pdcoefC_DIP-4*coefC*pdcoefA_DIP)) - (-coefB+sqrt(coefB^2-4*coefA*coefC))*2*pdcoefA_DIP ) / (2*coefA)^2;

		%dfProtA_dCI
        pdcoefA_CI = -(1/2)*EOpt_i*NProt;
        pdcoefB_CI = -(1/2)*EOpt_i*NM +EOpt_i*NL;
        pdcoefC_CI = 0;
        %dfProtA_dCI = ((2*coefA.*(-pdcoefB_CI + 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*pdcoefB_CI-4*pdcoefA_CI.*coefC)) -(-coefB +sqrt(coefB.^2-4*coefA.*coefC)).*2.*pdcoefA_CI)./(4*coefA.^2) );

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

		%dkST
		%dfProtA_dkST = dfProtA_dE*dE_dkST(i)+dfProtA_dCI*dCI_dkST(i);

        %dAColim_dQ10Photo = dAColim_dCI.*dCI_dQ10Photo(i)+dAColim_dE*dE_dCI(i)*dCI_dQ10Photo(i);
    end

%%
end
