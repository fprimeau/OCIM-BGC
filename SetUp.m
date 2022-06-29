global GC GO iter
iter = 0       ;
spd  = 24*60^2 ;
spa  = 365*spd ;
on   = true    ;
off  = false   ;
% addpath according to opterating system
if ismac
    addpath('~/Dropbox/myfunc'     )
    addpath('~/Documents/DATA/'    )
    addpath('~/Documents/DATA/OCIM')
elseif isunix
    addpath('/DFS-L/DATA/primeau/meganrs/DATA/')
	addpath('/DFS-L/DATA/primeau/weilewang/my_func/'  )
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/'     )
    addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
end

if GridVer == 90
    TRdivVer = 'Tv4' ;
elseif GridVer == 91
    switch(operator)
      case 'A'
        TRdivVer = 'CTL_He'             ; % choose this one %has along isopycnal diffusion
      case 'B'
        TRdivVer = 'CTL_noHe'           ;
      case 'C'
        TRdivVer = 'KiHIGH_He'          ;
      case 'D'
        TRdivVer = 'KiHIGH_noHe'        ;
      case 'E'
        TRdivVer = 'KvHIGH_KiLOW_He'    ;
      case 'F'
        TRdivVer = 'KvHIGH_KiLOW_noHe'  ;
      case 'G'
        TRdivVer = 'KiLOW_He'           ;
      case 'H'
        TRdivVer = 'KiLOW_noHe'         ;
      case 'I'
        TRdivVer = 'KvHIGH_He'          ;
      case 'J'
        TRdivVer = 'KvHIGH_noHe'        ;
      case 'K'
        TRdivVer = 'KvHIGH_KiHIGH_noHe' ;
    end
end

fprintf('Transport version: % s \n', TRdivVer)
if par.Cmodel == on
    fprintf('------- C model is on ------ \n')
    fprintf('DOP scaling factor is %2.2e \n', par.pscale)
    fprintf('DOC scaling factor is %2.2e \n', par.cscale)
	if ~isfield(par,'dynamicP')
		par.dynamicP = off;
		fprintf('--- default: P:C is a linear function of WOA observed phosphate ------ \n')
	elseif par.dynamicP == on
		fprintf('--- P:C is a linear function of modelled DIP ------ \n')
	else
		fprintf('--- P:C is a linear function of WOA observed phosphate ------ \n')
	end
end
if par.Omodel == on
    fprintf('------ O model is on ------- \n')
end
if par.Simodel == on
    fprintf('------ Si model is on ------ \n')
end
if par.Cellmodel == on
    fprintf('------ Cell stoichiometry model is on ------ \n')
	if ~isfield(par,'dynamicP')
		fprintf('--- default: Cell model depends on modelled DIP ------ \n')
		par.dynamicP = on;
	elseif par.dynamicP == on
		fprintf('--- Cell model depends on modelled DIP ------ \n')
	else
		fprintf('--- Cell model depends on observed nutrient fields ------ \n')
	end
end

fprintf('\n')

%
if GridVer == 90
    load transport_v4.mat grid M3d TR
    load M3d90x180x24v2.mat MSKS
    load Sobs_90x180x24.mat
    load tempobs_90x180x24.mat
    load po4obs_90x180x24.mat       % WOA PO4 observation [units: umol/kg]
	load no3obs_90x180x24.mat		% WOA NO3 obs [units: umol/kg]
    load Siobs_90x180x24.mat Siobs
    load Mouw_POC_90x180x24.mat  % sediment trap data MOUW
    %
    load GLODAPv2_talk.mat			% [units: umol/kg]
    load PME_TS_90x180x24.mat  pme
    load DICant_90x180x24.mat
    load GLODAPv2_90x180x24raw.mat  % GLODAP Nutrient units = [umol/kg]
	load raw_no3obs_90x180x24.mat   % GLODAP NO3
    load splco2_mod_monthly.mat     % monthly CO2 data
    load co2syspar90.mat co2syspar
    load cbpm_npp_annual_90x180.mat	% [units: mg C/m^2/day]
    load DOMobs_90x180x24.mat
    load kw660_90x180.mat
	PARobs = load('annual_PAR_90x180.mat'); %PAR data in units of [Einstein m-2 d-1] (units converted for cell model later in SetUp)
	load Kd490_MODIS_90x180.mat		% diffuse attenuation [m^-1]
    if ismac
        load MSK90/fixedPO_C.mat
        load MSK90/fixedPO_O2.mat
    elseif isunix
        load initCO_90x180x24.mat
    end
    grd = grid ;

elseif GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    %
    load M3d91x180x24.mat MSKS
    load Sobs_91x180x24.mat
    load po4obs_91x180x24.mat % WOA PO4 observation
	load no3obs_91x180x24.mat % WOA NO3 obs
    load tempobs_91x180x24.mat
    load Siobs_91x180x24.mat Siobs
    load Mouw_POC_91x180x24.mat  % sediment trap data MOUW
    %
    load GLODAPv2_talk_91x180x24.mat
    load PME_TS_91x180x24.mat  pme
    load DICant_91x180x24.mat
    load GLODAPv2_91x180x24raw.mat
	load raw_no3obs_91x180x24.mat  %GLODAP NO3
    load splco2_mod_monthly % monthly CO2 data
    load co2syspar91.mat co2syspar
    load cbpm_npp_annual_91x180.mat
    load DOMobs_91x180x24.mat
    load kw660_91x180.mat
	PARobs = load('annual_PAR_91x180.mat'); %PAR data in units of [Einstein m-2 d-1] (units converted for cell model later in SetUp)
	PARobs.par = PARobs.PAR;
	load Kd490_MODIS_91x180.mat		% diffuse attenuation [m^-1]
    if ismac
        load MSK91/CTL_He_PCO.mat
        load MSK91/CTL_He_PCO.mat
    elseif isunix
        load initCO_91x180x24.mat
    end
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
end

%---------------------- constants -------------------
par.spd = spd ;
par.spa = spa ;
[par.kw,par.P] = kw(M3d,grd);
par.Kw660 = Kw660   ;
par.p4    = p4      ;
par.c2p   = 110     ;	   % constant C:P ratio
par.rho   = 1024.5       ; % seawater density [kg/m^3];  	% consider using WOA18 density field instead of a constant?
permil    = par.rho*1e-3 ; % from umol/kg to mmol/m3;
par.permil = permil ;
% transiant CO2 concentraion;
par.year      = splco2_mod(:,1) ;
par.pco2_air  = splco2_mod(:,2) ;
par.pco2atm   = splco2_mod(1,2) ;
par.co2syspar = co2syspar       ;

%------------------- model grid info ----------------
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
dAt  = grd.DXT3d.*grd.DYT3d;
dVt  = dAt.*grd.DZT3d;
ARC  = MSKS.ARC     ;
MED  = MSKS.MED     ;
PAC  = MSKS.PAC     ;
ATL  = MSKS.ATL     ;
IND  = MSKS.IND     ;
iarc = find(ARC(:)) ;
imed = find(MED(:)) ;
ipac = find(PAC(:)) ;
iatl = find(ATL(:)) ;
iind = find(IND(:)) ;
%
par.dAt   = dAt     ;
par.dVt   = dVt     ;
par.M3d   = M3d     ;
par.iwet  = iwet    ;
par.nwet  = nwet    ;
par.TRdiv = -TR     ;
par.grd   = grd     ;
par.MSKS  = MSKS    ;
par.I = speye(nwet) ;

%------------------------  data  -----------------------
% get rid of arctic o2 observations
% DOCobs(iarc)  = nan ;
% DOCobs(imed)  = nan ;
alkraw(iarc)  = nan ;
alkraw(imed)  = nan ;
dicraw(iarc)  = nan ;
dicraw(imed)  = nan ;
% po4raw(iarc)  = nan ;
% po4raw(imed)  = nan ;
% sio4raw(iarc) = nan ;
% sio4raw(imed) = nan ;
% o2raw(iarc)   = nan ;
% o2raw(imed)   = nan ;
tempobs(tempobs(:)<-2.0) = -2.0 ;
par.Temp    = tempobs ;
par.Salt    = Sobs    ;
par.DSi     = Siobs   ;
par.po4obs  = po4obs  ;
par.no3obs  = no3obs  ;
par.o2raw   = o2raw   ;
par.po4raw  = po4raw  ;
par.no3raw  = no3raw  ; %no3raw field not in GLODAPv2_90x180x24raw.mat
par.sio4raw = sio4raw ;
par.DOCobs  = DOCobs  ;
par.alkraw  = alkraw*permil ;		% mmol/m^3
par.dicraw  = dicraw*permil ;		% mmol/m^3
par.dicant  = DICant*permil ;		% mmol/m^3
par.dopraw  = DOPobs - 0.03 ; % less refractory DOP
DOCclean   = RemoveRef(par) ;
ibad = find( DOCclean(iarc) > 50 ) ;
DOCclean(iarc(ibad)) = nan ;
par.docraw = DOCclean ;

%-------------------- prepare virtual flux -------------------
% PME part;
% [modT,modS,pme] = pme(par) ;
par.pme     = pme  ;
junk = M3d ;
junk(:,:,2:end) = 0 ;
isrf = find(junk(iwet)) ;
sdic = find(par.dicraw(iwet(isrf)) > 0);
salk = find(par.alkraw(iwet(isrf)) > 0);
% surface mean concentraions
par.sDICbar = sum(par.dicraw((iwet(isrf(sdic)))).* ...
                  dVt(iwet(isrf(sdic))))./sum(dVt(iwet(isrf(sdic))));
par.sALKbar = sum(par.alkraw((iwet(isrf(salk)))).* ...
                  dVt(iwet(isrf(salk))))./sum(dVt(iwet(isrf(salk))));

%-------------------- normalize temperature --------------------
for ji = 1:24
    t2d = par.Temp(:,:,ji);
    par.Temp(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
end
vT = par.Temp(iwet) ;
Tz = (vT - min(vT))./(max(vT) - min(vT)) ;
% Tz = zscore(vT)  ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par.Tz     = Tz*1e-8 ;
par.aveT   = nanmean(Tz3d(:,:,1:2),3) ;

%------------------- prepare for restoring ---------------------
% calculating global mean DIP, ALK, and DSi concentraions for
% restoring
idip = find(par.po4raw(iwet)>0)  ;
ialk = find(par.alkraw(iwet)>0)  ;
isil = find(par.sio4raw(iwet)>0) ;

par.DIPbar = sum(par.po4raw(iwet(idip)).*dVt(iwet(idip)))/sum(dVt(iwet(idip))) ;
par.ALKbar = sum(par.alkraw(iwet(ialk)).*dVt(iwet(ialk)))/sum(dVt(iwet(ialk))) ;
par.DSibar = sum(par.sio4raw(iwet(isil)).*dVt(iwet(isil)))/sum(dVt(iwet(isil)));

%------------------ Prepare Light field the model --------------------
par.nzo = 2;

PARobs = PARobs.par;

% fill in missing values along coastlines
PARsurf = inpaint_nans(PARobs);
SURF = M3d(:,:,1);
ilnd = find(SURF(:) == 0);
PARsurf(ilnd) = NaN;
PARsurf(PARsurf<=0) = min(PARobs(:));

% convert PAR [Einstein m^-2 d^-1] into units of [umol photon m^-2 s^-1] for cell model
PARsurf = PARsurf*10^6/spd; % PAR at surface
clear PARobs

% % PAR [Einstein m-2 d-1] is converted into units of [umol photon m^-2 s^-1] for cell model
% PARobs_PPFD = PARobs.par*10^6/spd; % PAR at surface
% clear PARobs
%
% % remove NaNs along coastlines in Light field (maybe move this into a seperate function)
% PARsurf = cleanPARobs(PARobs_PPFD,M3d); %local function at end of file % replaced with inpaint_nans
% clear PARobs_PPFD

% extrapolate light to bottom of euphotic zone
par.kI = 0.04;   % Light attenuation coefficient in seawater [m^-1]
%CESM light attenuation: average K for par wavelengths, plus absorbtion due to water; as a function of chlorophyll
%CHL = [mg/m^3]
%kI = 0.09*CHL.^0.4 ; % [m^-1]
%CbPM model uses Satellite diffuse attenuation coefficient
%	 median mixed layer light level = surface irradiance * exp (-k490 * MLD/2)
%	kI = k490
%	revised CbPM uses 9 different wavelengths
%par.kI = Kd490;

	PAR        = 0*M3d;
	for ii=1:par.nzo % only needed in euphotic zone for cell growth
		PAR(:,:,ii) = PARsurf.*exp(-par.kI*grd.zt(ii)); %PAR at mid depth of grid box [ii]
	end

par.PARobs = PAR;
		clear PAR
	% average PAR from surface to bottom of grid box 1
	%PAR(:,:,1)=PARobs_PPFD./(par.kI*grd.dzt(1)).*(1-exp(-par.kI*grd.dzt(1)));

%------------------ prepare NPP for the model --------------------
par.nzo = 2 ;
par.p2c = 0.006+0.0069*po4obs ;
inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan) = 0 ;

par.npp   = npp/(12*spd) ;		% units: mmol C/m^2/s
par.npp1  = (0.5*par.npp./grd.dzt(1)).*par.p2c(:,:,1) ; % units: mmol P/m^2/s
par.npp2  = (0.5*par.npp./grd.dzt(2)).*par.p2c(:,:,2) ;
par.Lambda = M3d*0 ;
par.Lambda(:,:,1) = 1./(1e-6+po4obs(:,:,1)) ;
par.Lambda(:,:,2) = 1./(1e-6+po4obs(:,:,2)) ;
par.Lambda(:,:,3:end) = 0 ;

%---------------------------- end ----------------------------------


% function PARsurf = cleanPARobs(PARobs_PPFD,M3d)
%     iwet1 = find(M3d(:,:,1));
%     %ibad = find(isnan(PARobs_PPFD(iwet1)));
%
%     PARsurf        = 0*M3d(:,:,1);
%     PARsurf(iwet1)=PARobs_PPFD(iwet1);
%
%     [ibady,ibadx] =find(isnan(PARsurf));
%
%     for ii=1:length(ibadx)
%         ilat = ibady(ii);
%         ilon = ibadx(ii);
%         if ilon>1 & ilon <180
%             nearby=PARsurf(ilat-1:ilat+1,ilon-1:ilon+1); % would be better to weight values on same latitude more
%
%             if length(nearby(nearby>0)) > 0
%                 PARsurf(ilat,ilon) = nanmean(nearby(nearby>0));
%             else
%                 nearby=PARsurf(ilat-2:ilat+2,ilon-2:ilon+2);
%                 if length(nearby(nearby>0)) > 0
%                     PARsurf(ilat,ilon) = nanmean(nearby(nearby>0));
%                 else
%                     keyboard
%                 end
%             end
%
%         else % if ilon==180
%             nearby=PARsurf(ilat-1:ilat+1,ilon-2:ilon);
%             if length(nearby(nearby>0)) > 0
%                 PARsurf(ilat,ilon) = nanmean(nearby(nearby>0));
%             else
%                 keyboard
%             end
%
%         end
%     end
% end
