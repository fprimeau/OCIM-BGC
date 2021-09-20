function par = SetPar(par)
    on   = true    ;  off  = false   ;
    spd  = 24*60^2 ;  spa  = 365*spd ;
    % fixed parameters
    par.kappa_g  = 1/(1e6*spa)  ; % geological restoring time [1/s];
    par.taup     = 30*24*60^2     ; % (s) pic dissolution time-scale
    % par.tau_TA   = 1./par.taup  ;
    par.kappa_p  = 1/(30*24*60^2) ;
	% the parameters kappa_p and b, which affect the sinking flux attenuation profile,
	% are not independently identified by the data. We therefore prescribe the value kappa_p of 1/(30 days)
    par.tauPIC = 30*spd ;
    par.kPIC   = 1/par.tauPIC ;
	% PIC dissolution constant 0.38 day^-1 based on first-order
    % reaction kinetics according to Sarmiento
    % and Gruber book (p.271);

    % load optimal parameters if they exist
    if isfile(par.fxhat) & par.LoadOpt == on
        load(par.fxhat)
    end
    % P model parameters
    if exist('xhat') & isfield(xhat,'sigma')
        par.sigma = xhat.sigma ;
    else
        par.sigma = 0.30 ; % default was 1.0e-01 , however 0.30 is used by (Wang, 2019: Nitrogen fixation)
    end
    if exist('xhat') & isfield(xhat,'kP_T')
        par.kP_T = xhat.kP_T ;
    else
        par.kP_T = 0.00 ;
    end
    if exist('xhat') & isfield(xhat,'kdP')
        par.kdP = xhat.kdP ;
    else
        par.kdP = 2.42e-08 ;
    end
    if exist('xhat') & isfield(xhat,'bP_T')
        par.bP_T = xhat.bP_T ;
    else
        par.bP_T = 0.00e+00 ;
    end
    if exist('xhat') & isfield(xhat,'bP')
        par.bP  = xhat.bP ;
    else
        par.bP  = 9.60e-01 ;
    end
    if exist('xhat') & isfield(xhat,'alpha')
        par.alpha = xhat.alpha ;
    else
        par.alpha = 3.37e-08   ;
    end
    if exist('xhat') & isfield(xhat,'beta')
        par.beta = xhat.beta ;
    else
        par.beta = 1.65e-02  ;
    end

    % C model parameters
    if exist('xhat') & isfield(xhat,'bC_T')
        par.bC_T = xhat.bC_T ;
    else
        par.bC_T =  0.00e+00 ;
    end
    if exist('xhat') & isfield(xhat,'bC')
        par.bC = xhat.bC ;
    else
        par.bC = 9.38e-01    ;
    end
    if exist('xhat') & isfield(xhat,'d')
        par.d = xhat.d   ;
    else
        par.d = 2.00e3    % default was par.d = 4.54e+03 ;
    end
    if exist('xhat') & isfield(xhat,'kC_T')
        par.kC_T = xhat.kC_T ;
    else
        par.kC_T = 0.00e+00 ;
    end
    if exist('xhat') & isfield(xhat,'kdC')
        par.kdC = xhat.kdC ;
    else
        par.kdC = 7.35e-08 ;
    end
    if exist('xhat') & isfield(xhat,'R_Si')
        par.R_Si = xhat.R_Si ;
    else
        par.R_Si = 0.00 ;
    end
    if exist('xhat') & isfield(xhat,'rR')
        par.rR = xhat.rR  ;
    else
        par.rR = 2.64e-02 ;
    end
    if exist('xhat') & isfield(xhat,'cc')
        par.cc = xhat.cc  ;
    else
        par.cc = 7.51e-4 ;
    end
    if exist('xhat') & isfield(xhat,'dd')
        par.dd = xhat.dd  ;
    else
        par.dd = 5.56e-03 ;
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
        par.rO2C = 1.10e+00 ;
    end
    if exist('xhat') & isfield(xhat,'O2P_T')
        par.O2P_T = xhat.O2P_T ;
    else
        par.O2P_T = 0.00 ;
    end
    if exist('xhat') & isfield(xhat,'rO2P')
        par.rO2P = xhat.rO2P ;
    else
        par.rO2P = 1.70e+02 ;
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
	%
	% Cell model parameters
	if exist('xhat') & isfield(xhat,'Q10Photo') % Q10 of photosynthesis
		par.BIO.Q10Photo = xhat.Q10Photo;
	else
		par.BIO.Q10Photo = 1.983;		% Q10 of photosynthesis
	end
	if exist('xhat') & isfield(xhat,'fStorage')
		par.BIO.fStorage = xhat.fStorage;
	else
		par.BIO.fStorage = exp(-.358);  % strength of luxury P storage [L/molC]
	end
	if exist('xhat') & isfield(xhat,'fRibE')
		par.BIO.fRibE = xhat.fRibE;
	else
		par.BIO.fRibE = .618;           % ribosome fraction of biosynthetic apparatus
	end
	if exist('xhat') & isfield(xhat,'kST0')
		par.BIO.kST0 = xhat.kST0;
	else
		par.BIO.kST0 =0.185;            % specific synthesis rate of synthetic apparatus at 25degC [1/hr]
	end
	if exist('xhat') & isfield(xhat,'PLip_PCutoff')
		par.BIO.PLip_PCutoff = xhat.PLip_PCutoff;
	else
		par.BIO.PLip_PCutoff = exp(-14.408);  % log of [P] below which more PLipids are substituted with Slipids
	end
	if exist('xhat') & isfield(xhat,'PLip_scale')
		par.BIO.PLip_scale = xhat.PLip_scale;
	else
		par.BIO.PLip_scale = 3.0e6;  % scale factor for logistic function controlling phospholipid quota (changed from 1.0 b.c. not using log(P). changed from 1e6 to 3e6 to make transition sharper)
	end
	if exist('xhat') & isfield(xhat,'PStor_rCutoff')
		par.BIO.PStor_rCutoff = xhat.PStor_rCutoff;
	else
		par.BIO.PStor_rCutoff = 2.25;  % radius [um] above which cell stores luxury phosphorus?
	end
	if exist('xhat') & isfield(xhat,'PStor_scale')
		par.BIO.PStor_scale = xhat.PStor_scale;
	else
		par.BIO.PStor_scale = 3.00;  % scale factor for logistic function controlling luxury phosphorus storage (changed default from 1 to 3 to give sharper transition)
	end
	if exist('xhat') & isfield(xhat,'alphaS')
		par.BIO.alphaS = xhat.alphaS;
	else
		par.BIO.alphaS = 0.225;          % radius at which cell is all periplasm and membrane [um]
	end

% cell model parameters that don't change
	if (par.Cellmodel==on)
		par.BIO.gammaDNA = .016;        % DNA fraction of cell
		par.BIO.gammaLipid = .173       % structural Lipid (non-membrane or periplasm) fraction of cell
		%par.BIO.lPCutoff = -7.252;		% log of max [P] for which Plipids will be substituted with Slipids
		%par.BIO.r0Cutoff = 2.25;		% % NEED TO REDEFINE: r0Cutoff =  rFullA; ASK GEORGE % now PStor_rCutoff
		par.BIO.DNT0 = 1e-12*3.6e2*3600;    % Diffusivity of Nitrate at 25degC [m^2/hr]
		par.BIO.DPT0 = 1e-12*3.6e2*3600;    % Diffusivity of Phosphate at 25degC [m^2/hr]
		par.BIO.Q10Diffusivity = 1.5;
		par.BIO.AMin =.05;              % minimal fraction of cell dry mass that is nutrient uptake proteins
		%par.BIO.CStor = 1.00;           %replaced by PStor_scale
		par.BIO.PhiS = .67;             % specific carbon cost of synthesis [gC/gC]
		%%% BIO parameters below should remain fixed
		par.BIO.pDry = .47;             % Dry mass fraction of the cell
		par.BIO.rho = 1e-12;            % cell density [g/um^3]
		par.BIO.fProtM = 0.25;          % protein fraction of cell membranes
		par.BIO.fProtL = .7;            % protein fraction of light harvesting apparatus
		par.BIO.PDNA = .095;            % phosphorus mass fraction in DNA [gP/g]
		par.BIO.PRib = 0.047;           % phosphorus mass fraction in ribosomes [gP/g]
		par.BIO.PPhospholipid = 0.042;  % phosphorus mass fraction in phospholipids [gP/g]
		par.BIO.NProt = .16;            % nitrogen mass fraction in proteins [gN/g]
		par.BIO.NDNA = .16;             % nitrogen mass fraction in DNA [gN/g]
		par.BIO.NRib = .16;             % nitrogen mass fraction in Ribosomes [gN/g]
		par.BIO.CProt = .53;            % carbon mass fraction in proteins [gC/g]
		par.BIO.CDNA = .36;             % carbon mass fraction in DNA [gC/g]
		par.BIO.CPhospholipid = .65;    % carbon mass fraction in phospholipids [gC/g] - why seperate form lipids?
		par.BIO.CLipid = .76;			% carbon mass fraction in other lipids (that are not phospholipids) [gC/g]
		par.BIO.CRib = .419;     		% carbon mass fraction in ribosomes [gC/g] (technically, only correct for eukaryotes)
		par.BIO.alphaPLip = 0.12;
	end
end
