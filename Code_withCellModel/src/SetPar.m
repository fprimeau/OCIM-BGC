function par = SetPar(par)
    on   = true    ;
    off  = false   ;
    sph  = 60^2    ;
    spd  = 24*sph  ;
    spa  = 365*spd ;

    % fixed parameters 
    par.kappa_g  = 1/(1e6*spa)  ; % geological restoring time [1/s] ;
    par.kappa_l  = 1/(12*sph )  ; % labile DOM remi time [1/s]     ;
    par.taup     = 30*spd       ; % (s) pic dissolution time-scale ;
    par.gamma    = 0 ; % 2/5 ;
    par.kappa_p  = 1/par.taup ;
    % PIC dissolution constant 0.38 day^-1 based on first-order
    % reaction kinetics according to Sarmiento
    % and Gruber book (p.271);
    par.tauPIC = 30*spd ; 
    par.kPIC   = 1/par.tauPIC ;

%-------C13---------
    % define C13 C14 fractionation factors and ratios in ocean and atmosphere
    % set up fractionation factors fro C13 and R13
    %  A. Schmittner et al.: Distribution of carbon isotope ratios (δ13C) in the ocean
    % 기존 code에는 par.c13.R13a = 0.01116303448. 그러나 출처 모름.
    % 아래는 -6.5 permil in 1859로 계산한 것. REF: A revised 1000 year atmospheric ı13C-CO2 record from Law Dome and South Pole, Antarctica. JGR Atmosphere, 2013, Rubino et al., doi:10.1002/jgrd.50668, 2013
    par.c13.R13a = 0.011164158200000;
    par.c13.alpha_k = 0.99915; % kenetic fractionation factor 
    par.c13.alpha_g2aq = 0.998764; % gas to water fractionation factor. It negelcts the minor temperature dependency of the isotopic fractionation factor from gaseous to aqueous CO2 (5x10-6 oC). Thus, it is a constant value corresponding to a mean temperature of 15oC.

    %------C14-----------
    par.c14.fc14 = 1.9 ;  % The fractionation  for 14C is 1.9 times the fractionation for 13C. REF: Reassessment of the 13C/12C and 14C/12C isotopic fractionaion ratio and its impact on high-precision radiocarbon dating, GCA (It can also be tuned.)
    par.lambda14 = 1/spa*log(2)/5730; % radiocarbon decay rate (yr^(-1) to s^(-1))
    par.c14.R14a = 1.220805*1e-12 ; %REF 찾기.
    par.c14.alpha_k    = 1 - (1 - 0.99915)*par.c14.fc14 ; % kenetic fractionation factor
    par.c14.alpha_g2aq = 1- (1 - 0.998764)*par.c14.fc14 ; % gas to water fractionation factor

        % load optimal parameters if they exist
    if isfile(par.fxhat) & par.LoadOpt == on 
        load(par.fxhat)
    end
    
    if exist('xhat') & isfield(xhat,'sigP')
        par.sigP = xhat.sigP ;
    else 
        par.sigP = 3.500e-1 ;
    end
    if exist('xhat') & isfield(xhat,'Q10P')
        par.Q10P = xhat.Q10P ;
    else 
        par.Q10P = 2.28 ; 
    end 
    if exist('xhat') & isfield(xhat,'kdP')
        par.kdP = xhat.kdP ;
    else 
        par.kdP = 1.86e-08; 
    end 
    if exist('xhat') & isfield(xhat,'bP_T')
        par.bP_T = xhat.bP_T ;
    else 
        par.bP_T = 0.319124 ;
    end 
    if exist('xhat') & isfield(xhat,'bP')
        par.bP  = xhat.bP ;
    else 
        par.bP  = 0.7701130e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'alpha')
        par.alpha = xhat.alpha ;
    else 
        par.alpha = 2.45e-08   ;
    end 
    if exist('xhat') & isfield(xhat,'beta')
        par.beta = xhat.beta ;
    else
        par.beta = 4.50e-01  ;
    end 

    % C model parameters
    if exist('xhat') & isfield(xhat,'sigC')
        par.sigC = xhat.sigC ;
    else 
        par.sigC = 0.90e-1 ;
    end
    if exist('xhat') & isfield(xhat,'kru')
        par.kru = xhat.kru ;
    else
        par.kru = 5.74e-12 ; % corresponding to 2000 years.
    end
    if exist('xhat') & isfield(xhat,'krd')
        par.krd = xhat.krd ;
    else
        par.krd = 2.99e-12 ; % corresponding to 2000 years.
    end
    if exist('xhat') & isfield(xhat,'etau')
        par.etau = xhat.etau ;
    else
        par.etau = 0.9800 ; 
    end 
    if exist('xhat') & isfield(xhat,'etad')
        par.etad = xhat.etad ;
    else
        par.etad = 0.978911524882553 ; 
    end 
    if exist('xhat') & isfield(xhat,'bC_T')
        par.bC_T = xhat.bC_T ;
    else
        par.bC_T =  0.957372e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'bC')
        par.bC = xhat.bC ;
    else
        par.bC = 6.50e-01    ;
    end
    if exist('xhat') & isfield(xhat,'d')
        par.d = xhat.d   ;
    else
        par.d = 4.55e+03 ;
    end 
    if exist('xhat') & isfield(xhat,'Q10C')
        par.Q10C = xhat.Q10C ;
    else 
        par.Q10C = 1.05e+00 ;
    end 
    if exist('xhat') & isfield(xhat,'kdC')
        par.kdC = xhat.kdC ;
    else 
        par.kdC =  5.42e-09; % from N nature paper,same as kdN;
    end
    if exist('xhat') & isfield(xhat,'R_Si')
        par.R_Si = xhat.R_Si ;
    else
        par.R_Si = 0.10 ;
    end
    if exist('xhat') & isfield(xhat,'rR')
        par.rR = xhat.rR  ;
    else
        par.rR = 2.34e-02 ;
    end
    if exist('xhat') & isfield(xhat,'cc')
        par.cc = xhat.cc  ;
    else
        par.cc = 8.38e-4 ;
    end 
    if exist('xhat') & isfield(xhat,'dd')
        par.dd = xhat.dd  ;
    else 
        par.dd = 8.83e-03 ;
    end 

    % O model parameters
    if exist('xhat') & isfield(xhat,'O2C_T')
        par.O2C_T = xhat.O2C_T ;
    else 
        par.O2C_T = 0.00 ;
    end 
    if exist('xhat') & isfield(xhat,'rO2C')
        par.rO2C = xhat.rO2C ;
    else 
        par.rO2C = 1.77e+00 ;
    end 
    %
    % Si model parameters
    if exist('xhat') & isfield(xhat,'dsi')
        par.dsi = xhat.dsi ;
    else
        par.dsi = 3300     ;
    end 
    if exist('xhat') & isfield(xhat,'at')
        par.at = xhat.at   ;
    else
        par.at = 1.32e16/spd;
    end 
    if exist('xhat') & isfield(xhat,'bt')
        par.bt = xhat.bt   ;
    else 
        par.bt = 11481     ;
    end 
    if exist('xhat') & isfield(xhat,'aa')
        par.aa = xhat.aa   ;
    else
        par.aa = 1         ;
    end 
    if exist('xhat') & isfield(xhat,'bb')
        par.bb = xhat.bb   ;
    else 
        par.bb = 0.968     ;
    end 

    % Cell model parameters
	if exist('xhat') & isfield(xhat,'Q10Photo')
        % Q10 of photosynthesis
		par.BIO.Q10Photo = real(xhat.Q10Photo);
    else
		par.BIO.Q10Photo = 1.46; 	% Anderson et al. 2021 (median for all phytoplankton; compilation of culture studies)
	end
	if exist('xhat') & isfield(xhat,'fStorage')
        % strength of luxury P storage [L/molC]
		par.BIO.fStorage = real(xhat.fStorage);
	else
		par.BIO.fStorage = 0.1;  
	end
	if exist('xhat') & isfield(xhat,'fRibE')
        % ribosome fraction of biosynthetic apparatus
		par.BIO.fRibE = real(xhat.fRibE);
    else
		par.BIO.fRibE = .60;			% biosynthesis = 60% ribosome; 40% protein
	end
	if exist('xhat') & isfield(xhat,'kST0')
        % specific synthesis rate of synthetic apparatus at 25degC [1/hr]
		par.BIO.kST0 = real(xhat.kST0);
	else
		par.BIO.kST0 =0.185;            % empirically derived specific efficiency: 0.168 hr^-1 (shuter 1979)								
	end
	if exist('xhat') & isfield(xhat,'PLip_PCutoff')
        % [P (mol/L)] below which more PLipids are substituted with Slipids (the x value of the logistic function sigmoid midpoint)
		par.BIO.PLip_PCutoff = real(xhat.PLip_PCutoff);
    else
		par.BIO.PLip_PCutoff = 1.0e-06;		
	end
	if exist('xhat') & isfield(xhat,'PLip_scale')
        % logistic growth rate or steepness of the curve for logistic function controlling phospholipid quota)
		par.BIO.PLip_scale = real(xhat.PLip_scale);
	else
		par.BIO.PLip_scale = 3.0e6;  
	end
	if exist('xhat') & isfield(xhat,'PStor_rCutoff')
         % radius [um] above which cell stores luxury phosphorus? (the x value of the logistic function sigmoid midpoint)
		par.BIO.PStor_rCutoff = real(xhat.PStor_rCutoff);
	else
		par.BIO.PStor_rCutoff = 2.25;
	end
	if exist('xhat') & isfield(xhat,'PStor_scale')
        % logistic growth rate or steepness of the curve for logistic function controlling luxury phosphorus storage
		par.BIO.PStor_scale = real(xhat.PStor_scale);
	else
		par.BIO.PStor_scale = 2.00; 
	end
	if exist('xhat') & isfield(xhat,'alphaS')
        % radius at which cell is all periplasm and membrane [um]
		par.BIO.alphaS = real(xhat.alphaS);
	else
		par.BIO.alphaS = 0.225;          
	end
	if exist('xhat') & isfield(xhat,'gammaDNA')
        % DNA Fraction of cell dry mass
		par.BIO.gammaDNA = real(xhat.gammaDNA);
	else
		par.BIO.gammaDNA = 0.016;          
	end

% fixed cell model parameters
	if (par.Cellmodel==on)
        par.BIO.alphaPLip 	= 0.12;		% scale factor for maximum phospholipid quota
		par.BIO.gammaLipid 	= .173;		% structural Lipid (non-membrane or periplasm) fraction of cell
		par.BIO.DNT0 		= 1e-12*3.6e2*3600;    % Diffusivity of Nitrate at 25degC [m^2/hr]
		par.BIO.DPT0 		= 1e-12*3.6e2*3600;    % Diffusivity of Phosphate at 25degC [m^2/hr]
		par.BIO.Q10Diffusivity = 1.5;	% Q10 temperature dependence of diffusivity
		par.BIO.Q10Bio 		= 2.0;		% Q10 temperature dependence of biosynthesis
		par.BIO.AMin 		= .05;		% minimal fraction of cell dry mass that is nutrient uptake proteins
		par.BIO.PhiS 		= .67;		% specific carbon cost of synthesis [gC/gC] (the cost of synthesizing functional organic material (taken to be 60% protein, 40% lipid) from glucose and nitrate; from Shuter, 1979)
		par.BIO.pDry 		= .47;		% Dry mass fraction of the cell
		par.BIO.rho 		= 1e-12;	% cell density [g/um^3]
		par.BIO.fProtM 		= 0.25;		% protein fraction of cell membranes
		par.BIO.fProtL 		= .7;		% protein fraction of light harvesting apparatus
		par.BIO.PDNA 		= .095;		% phosphorus mass fraction in DNA [gP/g]
		par.BIO.PRib 		= 0.047;	% phosphorus mass fraction in ribosomes [gP/g]
		par.BIO.PPhospholipid = 0.042;	% phosphorus mass fraction in phospholipids [gP/g]
		par.BIO.NProt 		= .16;		% nitrogen mass fraction in proteins [gN/g]
		par.BIO.NDNA 		= .16;		% nitrogen mass fraction in DNA [gN/g]
		par.BIO.NRib 		= .16;		% nitrogen mass fraction in Ribosomes [gN/g]
		%par.BIO.NPhospholipid = 0.008;	% nitrogen mass fraction in phospholipids [gP/g] (0.008 from Geider and La Roche 2002)
		par.BIO.CProt 		= .53;		% carbon mass fraction in proteins [gC/g]
		par.BIO.CDNA 		= .36;		% carbon mass fraction in DNA [gC/g]
		par.BIO.CPhospholipid = .65;	% carbon mass fraction in phospholipids [gC/g]
		par.BIO.CLipid 		= .76;		% carbon mass fraction in other lipids (that are not phospholipids) [gC/g]
		par.BIO.CRib 		= .419;		% carbon mass fraction in ribosomes [gC/g] (technically, only correct for eukaryotes)		
	end
end

