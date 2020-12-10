function [out,M] = CellCNP(par,x,P,N,T,Irr)
  % DESCRIPTION:
  % Previously named CPDual
  % This function calculates the C:P ratio of the maximally growing plankton
  % type in the environment specified by the function arguments
	% Currently: only relies on observed PO4,NO3,Temp,&PAR fields (does NOT use model DIP and DIN)

  % INPUTS:
  % Irr is light level (PAR, not total irradiance) in microeinsteins m^-2 s^-1 (i.e. umol photons m^-2 s^-1)
	% %%% PAR, but with units of photosynthetic photon flux density (PPFD) [μmole photons m-2 s-1]
  % T is temperature in Celsius
  % P is phosphate concentration in mol/L

  % AUTHOR: George Hagstrom. MODIFIED BY: Megan Sullivan

    on = true; off = false;
	% unpack the parameters to be optimized
	M = par.BIO;
	nbx = 0; % count number of C model tunable parameters;

	%Q10Photo
	if (M.opt_Q10Photo == on)
	    nbx = nbx + 1;
	    lQ10Photo  = x(par.pindx.lQ10Photo);
		M.Q10Photo  = exp(lQ10Photo);
	else
	    M.Q10Photo = M.Q10Photo;
	end
	%fRibE
	if (par.BIO.opt_fRibE == on)
		nbx = nbx + 1;
		lfRibE = x(par.pindx.lfRibE);
		M.fRibE = exp(lfRibE);
	else
		M.fRibE = M.fRibE;
	end
	%fStorage
	if (M.opt_fStorage == on)
		nbx = nbx + 1;
		lfStorage = x(par.pindx.lfStorage);
		M.fStorage = exp(lfStorage);
	else
		M.fStorage =M.fStorage;
	end
	%kST0
	if (par.BIO.opt_kST0 == on);
		nbx = nbx + 1;
		lkST0 = x(par.pindx.lkST0);
		M.kST0 = exp(lkST0);
	else
		M.kST0 =M.kST0;
	end

	% return number count back to neglogpost
	par.nbx = nbx;

  % Read in parameter values from par
		 %%% optimize these first %%%
     Q10Photo = M.Q10Photo;
     fRibE = M.fRibE;
     fStorage = M.fStorage;
	 kST0 = M.kST0;

	 alphaS = M.alphaS;
     gammaDNA = M.gammaDNA;
     gammaLipid = M.gammaLipid;
     lPCutoff = M.lPCutoff;
     r0Cutoff = M.r0Cutoff;
		%
     DNT0 = M.DNT0;
     DPT0 = M.DPT0;
     Q10Diffusivity = M.Q10Diffusivity;
     AMin = M.AMin;
     CStor = M.CStor;
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

	nanindx = find(isnan(Irr) | isnan(T) | isnan(N) |isnan(P));
	if ~isempty(nanindx)
		fprintf('Warning: [ %i ] NaNs found in Cell Model input fields. \n',length(nanindx))
	end

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

%%

% temperature dependent rates
    DN = DNT0*Q10Diffusivity.^((T-T0)./10.0);    % diffusivity of nitrate [ m^2/hr]
    DP = DPT0*Q10Diffusivity.^((T-T0)./10.0);     % diffusivity of phosphate [ m^2/hr]

    kST = kST0*2.0.^((T-T0)./10.0);  % specific synthesis rate of biosynthetic apparatus [1/hr]

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

    NM = fProtM*NProt;        % N fraction of membrane
                              %   = fraction of cell membrane that is protein * fraction of protein that is N
    CM = fProtM*CProt+fLipidM*CPhospholipid;     % carbon fraction of membrane
                                          % (membrane carbon from proteins + lipids)
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

% why is P fraction not defined for each of the Light, membrane, or other structure pools?

%% Solving for PLim
    aP= (3*DP.*P)*molarP/10^6 /(pDry*rho);  %units of P:[mol/L]*1e3*molarP->[g/m^3]*DP-> [g/m/hr]/10^6 -> [g/um/hr] % rho:[g/um^3] % alphaS:[um]
    coefA = aP.*CI.^2/alphaS^2-kST*PRib;
    coefB = -2*(1-gammaS)*CI.*aP/alphaS^2-gammaDNA*PDNA*kST;
    coefC = (1-gammaS)^2.*aP/alphaS^2;

    EPLim = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);

    rPLim = alphaS./(1-gammaS - CI.*EPLim); %1.0./( (1-gammaS)/alphaS-CI.*EPLim/alphaS);
    APLim = (1-gammaS -CI.*EPLim)./2; %should equal AMin
    fProtAPLim = AMin./APLim; %should equal 1 (if purely P limited = invest as much as you can in uptake proteins)

%%% derivatives of EPlim w.r.t. parameters
    dcoefA_dCI = 2.*aP.*CI/alphaS^2;
	dcoefB_dCI = -2*(1-gammaS).*aP/alphaS^2;

	dEPLim_dCI = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );
	dEPLim_dQ10Photo = dEPLim_dCI.*dCI_dQ10Photo;

%% Solving for NLim
    aN = 3*DN.*N*molarN/10^6 /(pDry*rho) ; %converting top to g/um/hr
    coefA = CI.^2.*aN/alphaS^2+kST.*CI*NM/2.0-kST*NL.*(CI-1)-kST*NE;
    coefB = -2*aN*(1-gammaS).*CI/alphaS^2-kST*gammaDNA*NDNA-kST*NM*(1-gammaS)/2.0-kST*AMin*NProt;
    coefC = aN.*(1-gammaS)^2/alphaS^2;

	ENLim = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);

    rNLim = alphaS./(1-gammaS - CI.*ENLim);
	ANLim = (1-gammaS-CI.*ENLim)./2; % mass fraction of periplasm in the cell. (= alphaS./(2*rNLim))
	fProtANLim = AMin./ANLim; % protein fraction of periplasm

%%% derivatives of ENlim w.r.t. parameters
    dcoefA_dCI = 2*CI.*aN/alphaS^2 + kST*NM/2.0 -kST*NL;
    dcoefB_dCI = -2*aN*(1-gammaS)/alphaS^2;

    dENLim_dCI = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );
	dENLim_dQ10Photo = dENLim_dCI.*dCI_dQ10Photo;


%% reset ANLim to AMin if the mass fraction of periplasm is too small (for large phytoplankton),
	rLG_ind = find(ANLim<AMin);  % equivalent to find(rNLim>rFullA);
        ANLim(rLG_ind) = AMin;
		fProtANLim(rLG_ind) = 1.0;

        coefA = CI.^2.*aN/alphaS^2+kST.*CI*(NM+NProt)/2.0-kST*NL.*(CI-1)-kST*NE;
        coefB = -2*aN*(1-gammaS).*CI/alphaS^2-kST*gammaDNA*NProt-kST*(NM+NProt)*(1-gammaS)/2.0;
        coefC = aN.*(1-gammaS)^2/alphaS^2;

        ENLim_LG = (-coefB-sqrt(coefB.^2-4*coefA.*coefC))./(2*coefA);
        rNLim_LG = alphaS./(1-gammaS - CI.*ENLim_LG);

		ENLim(rLG_ind) =ENLim_LG(rLG_ind);
		rNLim(rLG_ind) =rNLim_LG(rLG_ind);

    %%% derivatives of ENlim for large phyto w.r.t. parameters
        dcoefA_dCI = 2*CI.*aN/alphaS^2 + kST*(NM+NProt)/2.0 -kST*NL;
        dcoefB_dCI = -2*aN*(1-gammaS)/alphaS^2;

        dENLim_dCI_LG = ((2*coefA.*(-dcoefB_dCI - 1./(2*sqrt(coefB.^2-4*coefA.*coefC)).*(2*coefB.*dcoefB_dCI-4*dcoefA_dCI.*coefC)) -(-coefB-sqrt(coefB.^2-4*coefA.*coefC)).*2.*dcoefA_dCI)./(4*coefA.^2) );
        dENLim_dQ10Photo_LG = dENLim_dCI.*dCI_dQ10Photo;

        dENLim_dCI(rLG_ind) = dENLim_dCI_LG(rLG_ind);
        dENLim_dQ10Photo(rLG_ind) = dENLim_dQ10Photo_LG(rLG_ind);

%% Calculate growth rates (biosynthetic rate) under each condition
    muNLim = kST.*ENLim;  % units: [1/hr]
    muPLim = kST.*EPLim;  % units: [1/hr]

    muNLim_P= aP.*fProtANLim./(rNLim.^2.*(ENLim*PRib+gammaDNA*PDNA));
    muPLim_N= aN./(rPLim.^2.*((CI-1)*NL.*EPLim+EPLim*NProt+gammaDNA*NDNA+alphaS./(2.0.*rPLim)*NProt+alphaS./(2.0.*rPLim)*NM));
			% muNLim_P is the "P-limited" growth rate at the strategy where
			% N-limitation, Biosynthesis, and Photosynthesis balance. Likewise,
			% muPLim_N is the "N-Limited" growth rate at the strategy where
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
			% strategy is suboptimal, and vice versa for the P-lim strategy.


%% initialize vector fields
all_nan = NaN(size(P));

%CP = all_nan;
%NP = all_nan;
%CN = all_nan;
muOpt = all_nan;
rOpt = all_nan;
AOpt = all_nan;
fProtAOpt = all_nan;
EOpt = all_nan;
LimType = all_nan;

QP = all_nan;
QC = all_nan;
QN = all_nan;

dC2P_dQ10Photo = all_nan;
dE_dQ10Photo = all_nan;
dfProtAOpt_dCI = all_nan;
dfProtAOpt_dE = all_nan;
dE_dCI = all_nan;

%% loop through points to calculate cell quotas
for i =1:length(P)
if ~isnan(Irr(i)) & ~isnan(T(i)) & ~isnan(P(i)) & ~isnan(N(i)) & (P(i)>0)
	disp(i)
    if muNLim(i)<muPLim(i)
        if muNLim(i)<muNLim_P(i) % N Limitation
            LimState = 0;
			CIi = CI(i);
			kSTi = kST(i);
            EOpt_i = ENLim(i);
            rOpt_i = rNLim(i);
            fProtAOpti = fProtANLim(i);

            %%% derivatives
            dE_dCI(i) = dENLim_dCI(i);
			dE_dQ10Photo(i) = dENLim_dQ10Photo(i);
			if ANLim(i)>AMin
				dfProtAOpt_dE(i) = fProtAOpti^2*CIi/(2*AMin); %if ANLim> AMin
				dfProtAOpt_dCI(i) = fProtAOpti^2/(2*AMin)*(EOpt_i);
			else
				dfProtAOpt_dE(i) = 0;
				dfProtAOpt_dCI(i) = 0;
            end

        else
            LimState = 2;
			CIi = CI(i);
			kSTi = kST(i);
			aNi = aN(i);
			aPi = aP(i);

            EMin=1e-4;
            EMax=(1-gammaS)/CIi;
			[EOpt_i, dE_dCI(i)] = EColim(EMin,EMax);
			dE_dQ10Photo(i) = dE_dCI(i)*dCI_dQ10Photo(i);

            rOpt_i = alphaS./(1-gammaS - CIi.*EOpt_i);

			[fProtAOpti, dfProtAOpt_dCI(i), dfProtAOpt_dE(i)] = AColim(EOpt_i);

            %muPColim = fProtAOpt_i.*aP/rOpt_i.^2/((EOpt*fRibE*PRib+(gammaDNA)*PDNA ));
            %muNColim = aNi./rOpt_i.^2/(EOpt*NE+EOpt.*(CIi-1)*NL+gammaS*NS+fProtAOpt_i*alphaS./(2.0*rOpt_i)*NProt+alphaS/(2.0*rOpt_i)*NM);
        end
    elseif muPLim(i) <= muNLim(i)
        if muPLim(i)<muPLim_N(i)
            LimState = 1;
            EOpt_i = EPLim(i);
            rOpt_i = rPLim(i);
            %fProtAOpti = 1;
			fProtAOpti = fProtAPLim(i); % = 2*AMin/(1-gammaS -CI.*EPLim)
			CIi = CI(i);
			kSTi = kST(i);

            dE_dCI(i) = dEPLim_dCI(i);
			dE_dQ10Photo(i) = dEPLim_dQ10Photo(i);
			if fProtAOpti<1
				dfProtAOpt_dE(i) = fProtAOpti^2*CIi/(2*AMin);
				dfProtAOpt_dCI(i) = fProtAOpti^2/(2*AMin)*(EOpt_i); %dE_dQ10Photo(i)/dCI_dQ10Photo
			else
				dfProtAOpt_dE(i) = 0;
				dfProtAOpt_dCI(i) = 0;
            end

        else
            LimState = 3; %same as 2. changed for debugging
			CIi = CI(i);
			kSTi = kST(i);
			aNi = aN(i);
			aPi = aP(i);

            EMin = 1e-4;
            EMax = (1-gammaS)/CIi;
			[EOpt_i, dE_dCI(i)] = EColim(EMin,EMax);
			dE_dQ10Photo(i) = dE_dCI(i)*dCI_dQ10Photo(i);
            %dE_dQ10Photo_dQ10Photo = dE_dCI(i)*dCI_dQ10Photo_dQ10Photo(i) + dCI_dQ10Photo(i)* dEColim_dCI_dQ10Photo

            rOpt_i = alphaS./(1-gammaS - CIi.*EOpt_i);

			[fProtAOpti, dfProtAOpt_dCI(i), dfProtAOpt_dE(i)] = AColim(EOpt_i);


            %muPColim = fProtAOpti.*aPi/rOpt_i^2/((EOpt_i*fRibE*PRib+(gammaDNA)*PDNA ));
            %muNColim = aN./rOpt_i.^2/(EOpt_i*NE+EOpt_i.*(CIi-1)*NL+gammaS*NS+fProtAOpti*alphaS./(2.0*rOpt_i)*NProt+alphaS/(2.0*rOpt_i)*NM);
        end
    end % end cases

%%% save  cell quotas of P, N, and C
		QP(i) =((EOpt_i*fRibE*PRib+(gammaDNA)*PDNA ))/molarP;
		QN(i) = (EOpt_i*NE+EOpt_i.*(CIi-1)*NL+gammaS*NS+fProtAOpti*alphaS./(2.0*rOpt_i)*NProt+alphaS/(2.0*rOpt_i)*NM)/molarN;
		QC(i) = (EOpt_i*CE+EOpt_i.*(CIi-1)*CL+gammaS*CS+fProtAOpti*alphaS./(2.0*rOpt_i)*CProt+alphaS/(2.0*rOpt_i)*CM)/molarC .*(1+24.0*.25*kSTi.*EOpt_i*(1+PhiS));

		EOpt(i) = EOpt_i;
		rOpt(i) = rOpt_i;
		AOpt(i) = AMin./fProtAOpti;
        fProtAOpt(i) = fProtAOpti;
		muOpt(i) = kSTi.*EOpt_i;
		rOpt(i) = rOpt_i;
        LimType(i) = LimState;


end % if ~isnan(...)
end % end for loop

CP = QC./QP;
NP = QN./QP;
CN = QC./QN;
%C2P = (molarP/molarC)*((EOpt*CE+EOpt.*(CI-1)*CL+gammaS*CS+(AMin./AOpt).*(1-gammaS - CI.*EOpt)/2.0*CProt+ (1-gammaS - CI.*EOpt)/2.0*CM).*(1+24.0*.25*kST.*EOpt*(1+PhiS)) ) ./ (EOpt*fRibE*PRib+(gammaDNA)*PDNA); % only tiny random error (10^-13) from CP

% ****** Calculate derivatives of C:P ****** %

        pdrOpt_E = rOpt.^2./alphaS.*CI;
        pdrOpt_CI =rOpt.^2./alphaS.*EOpt;

        dfProtAOpt_dQ10Photo = dfProtAOpt_dE.*dE_dQ10Photo + dfProtAOpt_dCI.*dCI_dQ10Photo;
        drOpt_dQ10Photo = pdrOpt_E.*dE_dQ10Photo + pdrOpt_CI.*dCI_dQ10Photo;

		% QC = QC(E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kSTi,PhiS)
		% Take partial derivatives with respect to all these variables
		pdQC_E = (CE + (CI-1)*CL).*(1+24.0*.25*kST.*EOpt.*(1+PhiS))./molarC + (EOpt.*CE+EOpt.*(CI-1).*CL+gammaS*CS+fProtAOpt.*alphaS./(2.0*rOpt).*CProt+alphaS./(2.0*rOpt)*CM).*24.0*0.25.*kST*(1+PhiS)/molarC;
		pdQC_CE = 0; % do later
		pdQC_CI = (EOpt*CL).*(1+24.0*.25*kST.*EOpt*(1+PhiS))/molarC;
		pdQC_CL = 0; % do later
		pdQC_gammaS = 0; % do later
		pdQC_CS = 0; % do later
		pdQC_fProtAOpt = (alphaS./(2.0*rOpt)*CProt).*(1+24.0*.25*kST.*EOpt.*(1+PhiS))/molarC;
		pdQC_alphaS = 0; % do later
		pdQC_rOpt = (-fProtAOpt*alphaS./(2.0*rOpt.^2)*CProt-alphaS./(2.0*rOpt.^2)*CM).*(1+24.0*.25*kST.*EOpt.*(1+PhiS))/molarC;
		pdQC_CProt = 0; % do later
		pdQC_CM = 0; % do later
		pdQC_kSTi = 0; % do later
		pdQC_PhiS = 0; %  do later
		% QP = QP(E,fRibE,PRib,gammaDNA,PDNA)
		% Take partial derivatives with respect to all these variables
		pdQP_E = fRibE*PRib/molarP;
		pdQP_fRibE = EOpt.*PRib/molarP;
		pdQP_PRib = EOpt.*fRibE/molarP;
		pdQP_gammaDNA = PDNA/molarP;
		pdQP_PDNA = gammaDNA/molarP;
		% C:P = C:P(E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kSTi,PhiS,fRibE,PRib,gammaDNA,PDNA)
		% Use the partial derivatives of QC and QP to calculate the partial derivatives of C:P. pdC2P_E is the only complicated ones
		pdC2P_E = (QP.*pdQC_E - QC.*pdQP_E)./QP.^2;
		pdC2P_CE = pdQC_CE./QP;
		pdC2P_CI = pdQC_CI./QP;
		pdC2P_CL = pdQC_CL./QP;
		pdC2P_gammaS = pdQC_gammaS./QP;
		pdC2P_CS = pdQC_CS./QP;
		pdC2P_fProtAOpt = pdQC_fProtAOpt./QP;
		pdC2P_alphaS = pdQC_alphaS./QP;
		pdC2P_rOpt = pdQC_rOpt./QP;
		pdC2P_CProt = pdQC_CProt./QP;
		pdC2P_CM = pdQC_CM./QP;
		pdC2P_kSTi = pdQC_kSTi./QP;
		pdC2P_PhiS = pdQC_PhiS./QP;
		pdC2P_fRibE = -pdQP_fRibE.*QC./QP.^2;
		pdC2P_PRib = -pdQP_PRib.*QC./QP.^2;
		pdC2P_gammaDNA = -pdQP_gammaDNA.*QC./QP.^2;
		pdC2P_PDNA = -pdQP_PDNA.*QC./QP.^2;
		% Now we need to know the derivatives of E,CE,CI,CL,gammaS,CS,fProtAOpt,alphaS,rOpt,CProt,CM,kSTi,PhiS,fRibE,PRib,gammaDNA,PDNA with respect to the variables we are changing

		dC2P_dQ10Photo = pdC2P_E.*dE_dQ10Photo + pdC2P_CI.*dCI_dQ10Photo + pdC2P_fProtAOpt.*dfProtAOpt_dQ10Photo + pdC2P_rOpt.*drOpt_dQ10Photo;
        %d2C2P_dQ10Photo2 =pdC2P_E.*dE_dQ10Photo_dQ10Photo + dE_dQ10Photo.*dC2P_dE_dQ10Photo + pdC2P_CI.*dCI_dQ10Photo_dQ10Photo +dC2P_dCI_dQ10Photo.*dCI_dQ10Photo;

%% add PStorage
%%% Why is Pstorage tacked on at the end? %%%
%%% dont know how to link P storage to fitness
% could be reason to have max for PStorage - places with a lot of PStorage are often places with Iron limitation (but no iron limitation in model)
% if P storage is taken into account, model will think cell is p limited.
% E.G. southern ocean: P and N repllete, but P is drawn down a little more

% alphaPLip = 0.12 (fraction of membrane that is phopholipids)

        PLip = alphaPLip*.5./rOpt*PPhospholipid/CPhospholipid*molarC/molarP./(1+exp((P./lPCutoff - 1)))   ;
        % PLip = how much P from phospholipis is in cell
				% magnitude of PLip seems wrong (~5*10^-7)
				% magnitude of P ~ 1
				%magnitude of rOpt ~1000
        PStorage = fStorage*5000*P./(1+exp(-CStor.*(rOpt-r0Cutoff)))+PLip;% units [P/C]

        CPwithPStor = 1./(1./CP+PStorage);
        NPwithPStor = 1./(1./NP + PStorage.*CN);

%% Calculate C:P derivatives including PStorage %%%
		dPLip_drOpt = -PLip./rOpt;
		dPStorage_drOpt = -(PStorage-PLip)./(1+exp(-CStor.*(rOpt-r0Cutoff))).*(-CStor.*exp(-CStor.*(rOpt-r0Cutoff))) +dPLip_drOpt;
		% dQ10Photo
		dPStorage_dQ10Photo = dPStorage_drOpt.* drOpt_dQ10Photo;
		dC2Pstor_dQ10Photo = -CPwithPStor.^2 .*(-1./CP.^2 .*dC2P_dQ10Photo +dPStorage_dQ10Photo);
		% dfStorage
		dPStorage_dfStorage = 5000*P./(1+exp(-CStor.*(rOpt-r0Cutoff)));
		dC2Pstor_dfStorage = -CPwithPStor.^2 .*dPStorage_dfStorage;


%%% change Nans to zeros (trying to fix error in newton solver)
%CP(isnan(CP)) = 0;

%% Output
    out.CP = CPwithPStor;
    out.CPnostor = CP;
    out.NP = NPwithPStor;
    out.CN = CN;
    out.LimType = LimType;
    out.r = rOpt;
    out.E = EOpt;
	out.A = AOpt;
    out.mu = muOpt;
	out.PStorage = PStorage;
	out.dC2Pnostor_dQ10Photo = dC2P_dQ10Photo;
	out.dC2P_dQ10Photo = dC2Pstor_dQ10Photo;
	%out.d2C2P_dQ10Photo2;
	out.dC2P_dfStorage = dC2Pstor_dfStorage;

	disp('Cell Model Done')
    % L, A, S?
    %save par values

%% define functions
	function [EColim, dEColim_dCI] = EColim(EMin,EMax) % rhsFuncColim(EMin,EMax)
        % polynomial coefficients
		c3 = kSTi^2*NProt*PE + (aPi*aNi/alphaS^4)*CIi^3 - (aPi/alphaS^2)*kSTi*CIi*(NE+(CIi-1)*NL-0.5*CIi*NM);
		% function to return all derivs of c3 wrt params
		c2 = kSTi^2*NProt*gammaS*PS - 3*(aPi*aNi/alphaS^4)*CIi^2*(1-gammaS) + (aPi/alphaS^2)*kSTi*(1-gammaS)*(NE+(CIi-1)*NL-0.5*CIi*NM) - (aPi/alphaS^2)*kSTi*CIi*(gammaS*NS+0.5*(1-gammaS)*NM);
		c1 = 3*(aPi*aNi/alphaS^4)*CIi*(1-gammaS)^2 + (aPi/alphaS^2)*kSTi*(1-gammaS)*(gammaS*NS+0.5*(1-gammaS)*NM);
		c0 = -(aPi*aNi/alphaS^4)*(1-gammaS)^3;
		%rhsFuncColim = c3*E.^3 + c2*E.^2 + c1*E +c0;

        % solve for the zeros of rhsFuncColim
        rhsFuncColim = @(E) (c3*E^3 + c2*E^2 + c1*E +c0);
        E0 = (EMin+EMax)/2; %starting point for solver

        sol = fsolve_cmplx(rhsFuncColim,E0);
		EColim = sol;

		dc3_dCI = 3*(aPi*aNi/alphaS^4)*CIi^2 - (aPi/alphaS^2)*kSTi*(NE+(CIi-1)*NL-0.5*CIi*NM) - (aPi/alphaS^2)*kSTi*CIi*(NL-0.5*NM);
		dc2_dCI = -6*(aPi*aNi/alphaS^4)*CIi*(1-gammaS) + (aPi/alphaS^2)*kSTi*(1-gammaS)*(NL-0.5*NM) - (aPi/alphaS^2)*kSTi*(gammaS*NS+0.5*(1-gammaS)*NM);
		dc1_dCI = 3*(aPi*aNi/alphaS^4)*(1-gammaS)^2;
		dc0_dCI = 0;
		dEColim_dCI = -(EColim^3*dc3_dCI+EColim^2*dc2_dCI+EColim*dc1_dCI+dc0_dCI)/(3*EColim^2*c3 +2*EColim*c2 +c1);
    end

	function dEColim_dCI_dQ10Photo = dEColim_dCI_dQ10Photo(E)
		dEColim_dCI_dQ10Photo=dEColim_dCI_dCI*dCI_dQ10Photo +dEColim_dCI_dE*dE_dQ10Photo
	end

%dEdalphaS_colim

    function [AColim, dA_dCI, dA_dE]  =AColim(E)
        % amount of A (fraction periplasm devoted to protein)
        % uptake rate of phosphorous = linear function of how much prot in periplasm
        %AColim = (-2*aN*(-1+CI*E+gammaS)^2+alphaS^2*E*kST*(2*E*(NE+(CI-1)*NL)+NM*(1-CI*E-gammaS)+2*gammaDNA*NProt))/(alphaS^2*E*(CI*E-1+gammaS)*kST*NProt);
		AColim = (-2*aNi*(gammaS-1+CIi.*E).^2+alphaS^2.*E*kSTi.*(2.*E.*(NE+(CIi-1)*NL)+NM*(1-gammaS-CIi.*E)+2*gammaS*NS))./(alphaS^2.*E.*(CIi.*E-1+gammaS)*kSTi*NProt);

		Atop = -2*aNi*(gammaS-1+CIi.*E).^2 + 2*alphaS^2*kSTi*E^2*(NE+(CIi-1)*NL) +alphaS^2*kSTi*E*(1-gammaS-CIi*E)*NM + 2*alphaS^2*kSTi*E*gammaS*NS;
		Abot = (alphaS^2.*E.*(CIi.*E-1+gammaS)*kSTi*NProt);

		dtop_dCI = -4*aNi*(gammaS-1+CIi.*E)*E + 2*alphaS^2*kSTi*E^2*NL + alphaS^2*kSTi*E*(-E)*NM;
		dbot_dCI = alphaS^2.*E.*E*kSTi*NProt;

		dA_dCI = (Abot*dtop_dCI - Atop*dbot_dCI)/(Abot^2);

		dtop_dE = -4*aNi*(gammaS-1+CIi.*E)*CIi + 4*alphaS^2*kSTi*E*(NE+(CIi-1)*NL) + alphaS^2*kSTi*(1-gammaS-CIi*E)*NM + alphaS^2*kSTi*E*(-CIi*NM) + 2*alphaS^2*kSTi*gammaS*NS;
		dbot_dE = alphaS^2*(CIi.*E-1+gammaS)*kSTi*NProt + alphaS^2.*E.*CIi*kSTi*NProt;

		dA_dE = (Abot*dtop_dE - Atop*dbot_dE)/(Abot^2);
    end

%%
end
