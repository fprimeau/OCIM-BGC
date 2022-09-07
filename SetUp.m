global GC GC13 GO iter
iter = 0       ;
spd  = 24*60^2 ;
spa  = 365*spd ;
on   = true    ;
off  = false   ;
% addpath according to opterating system
if ismac 
    addpath('../myfunc/'  )
    addpath('../DATA/'    )
elseif isunix
    addpath('/DFS-L/DATA/primeau/oceandata/OCIM_BGC_DATA/myfunc/'  )
    addpath('/DFS-L/DATA/primeau/oceandata/OCIM_BGC_DATA/DATA/'     )
    % addpath('/DFS-L/DATA/primeau/oceandata/OCIM_BGC_DATA/DATA/OCIM' )
end

if GridVer == 90
    TRdivVer = 'Tv4' ;
elseif GridVer == 91 
    switch(operator)
      case 'A'
        TRdivVer = 'CTL_He'             ;
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
    fprintf('-------- C model is on -------- \n')
    fprintf('DOP scaling factor is %2.2e \n', par.pscale)
    fprintf('DOC scaling factor is %2.2e \n', par.cscale)
end
if par.Omodel == on
    fprintf('-------- O model is on -------- \n')
end 
if par.Simodel == on
    fprintf('-------- Si model is on -------- \n')
end 
fprintf('\n')

%
if GridVer == 90
    load transport_v4.mat grid M3d TR
    load M3d90x180x24v2.mat MSKS 
    load teng_region_90x180.mat R
    load tempobs_90x180x24.mat
    load Sobs_90x180x24.mat
    load po4obs_90x180x24.mat       % WOA PO4 observations
    load o2obs_90x180x24.mat        % WOA O2 observations
    load Siobs_90x180x24.mat Siobs
    load Mouw_POC_90x180x24.mat  % sediment trap data MOUW
    load initCO_90x180x24.mat
    %
    load GLODAPv2_talk.mat 
    load PME_TS_90x180x24.mat  pme  
    load DICant_90x180x24.mat
    load GLODAPv2_90x180x24raw.mat
    load splco2_mod_monthly.mat     % monthly CO2 data
    load co2syspar90.mat co2syspar
    load NPP_VGPM_Mar252022_90x180.mat
    % load cbpm_npp_annual_90x180.mat
    load DOPobs_90x180x24.mat
    load DOCobs_Feb2022_90x180x24.mat
    load kw660_90x180.mat
    alkraw = talk ;
    load initCO_90x180x24.mat
    
    grd = grid ;
    
elseif GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output') ;
    %
    % load river_flux_DOCALK.mat doc_riv_flux alk_riv_flux dic_riv_flux
    load M3d91x180x24.mat MSKS
    load teng_region_91x180.mat R
    load Sobs_91x180x24.mat
    load po4obs_91x180x24.mat   % WOA PO4 observations
    load o2obs_91x180x24.mat    % WOA O2 observations
    load tempobs_91x180x24.mat
    load Siobs_91x180x24.mat Siobs
    load Mouw_POC_91x180x24.mat  % sediment trap data MOUW
    %
    load GLODAPv2_talk_91x180x24.mat
    load PME_TS_91x180x24.mat  pme 
    load DICant_91x180x24.mat
    load GLODAPv2_91x180x24raw.mat
    load splco2_mod_monthly % monthly CO2 data
    load co2syspar91.mat co2syspar
    % load SeaWiFS_CbPM_May172022_91x180.mat
    load MODIS_CbPM_May192022_91x180.mat
    % load mkFigs/gamma_MODIS_CbPM.mat gamma 
    % load NPP_VGPM_Mar252022.mat
    load DOPobs_91x180x24.mat DOPobs 
    load DOCobs_Feb2022_91x180x24 DOCobs
    load kw660_91x180.mat
    load initCO_91x180x24.mat % initial guess for C and O fileds

    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
end

%---------------------- constants -------------------
% par.gamma = gamma ;
par.spd = spd ;
par.spa = spa ;
[par.kw, par.P] = kw(M3d, grd) ;
par.Kw660 = Kw660   ;
par.p4    = p4      ;
par.c2p   = 110     ;
par.rho   = 1024.5       ; % seawater density;
permil    = par.rho*1e-3 ; % from umol/kg to mmol/m3; 
par.permil = permil      ;
% transiant CO2 concentraion ;
par.year      = splco2_mod(:,1) ;
par.pco2_air  = splco2_mod(:,2) ;
par.pco2atm   = splco2_mod(1,2) ;
par.co2syspar = co2syspar       ;

%------------------- model grid info ----------------
ilnd = find(M3d(:) == 0);
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
par.R     = R       ;
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
% get rid of Arctic o2 observations
DOCobs(iarc)  = nan ;
DOCobs(imed)  = nan ;
alkraw(iarc)  = nan ;
alkraw(imed)  = nan ; 
dicraw(iarc)  = nan ;
dicraw(imed)  = nan ;
% po4raw(iarc)  = nan ;
% po4raw(imed)  = nan ; 
% o2raw(iarc)   = nan ;
% o2raw(imed)   = nan ;
% sio4raw(iarc) = nan ;
% sio4raw(imed) = nan ;

% remove outlayers of DOC obs
DOCobs = rmOutliers(DOCobs, M3d) ;
DOPobs = rmOutliers(DOPobs, M3d) ;
po4raw = rmOutliers(po4raw, M3d) ;
dicraw = rmOutliers(dicraw, M3d) ;
alkraw = rmOutliers(alkraw, M3d) ;
o2raw  = rmOutliers(o2raw, M3d)  ;

tempobs(tempobs(:)<-2.0) = -2.0 ;
po4obs(po4obs(:)<0.05) = 0.05   ;
po4raw(po4raw(:)<0.05) = nan    ;      
for ji = 1:24
    p2d = po4obs(:,:,ji);
    po4obs(:,:,ji) = smoothit(grd,M3d,p2d,3,1e5);
end

par.Temp    = tempobs ;
par.Salt    = Sobs    ;
par.DSi     = Siobs   ;
par.po4obs  = po4obs  ;
par.o2raw   = o2raw   ;
par.po4raw  = po4raw  ;
par.sio4raw = sio4raw ;
par.DOCobs  = DOCobs  ;
par.alkraw  = alkraw*permil ;
par.dicraw  = dicraw*permil ; 
par.dicant  = DICant*permil ;
par.dopraw  = DOPobs - 0.05 ; % less refractory DOP 
par.docraw  = DOCobs ;
% par.doc_riv_flux = doc_riv_flux ;
% par.alk_riv_flux = alk_riv_flux ;
% par.dic_riv_flux = dic_riv_flux ;
%-------------------- prepare virtual flux -------------------
% PME part;
% [modT,modS,pme] = pme(par) ; 
par.pme = pme  ;   
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
% add + 1.0 to prevent from getting infinit kP or kC 
Tz = (vT - min(vT) + 1.0)./(max(vT) - min(vT)) ;
par.Tz = Tz*1e-8 ;
par.vT = vT      ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par.aveT   = nanmean(Tz3d(:,:,1:3),3) ;

%-------------------- correct o2 concentration --------------------
o2obs_c = M3d*0;
o2obs(iwet) = o2obs(iwet).*44.661;        % convert unit form [ml/l] to [umol/l].
o2obs_c(iwet) = o2obs(iwet).*1.009-2.523; % o2 correction based on Bianchi et al.(2012) [umol/l] .
ineg = find(o2obs_c(:)<0);                % find out negative values and set them to zero.
o2obs_c(ineg) = 0.1 ;
par.o2obs = o2obs_c ;
%---------------------- prepare for restoring -----------------------
% calculating global mean DIP, ALK, and DSi concentraions for
% restoring 
idip = find(par.po4raw(iwet) > 0)  ;
ialk = find(par.alkraw(iwet) > 0)  ;
isil = find(par.sio4raw(iwet)> 0) ;

par.DIPbar = sum(par.po4raw(iwet(idip)).*dVt(iwet(idip)))/sum(dVt(iwet(idip))) ;
par.ALKbar = sum(par.alkraw(iwet(ialk)).*dVt(iwet(ialk)))/sum(dVt(iwet(ialk))) ;
par.DSibar = sum(par.sio4raw(iwet(isil)).*dVt(iwet(isil)))/sum(dVt(iwet(isil)));

%-------------------- prepare NPP for the model ----------------------
par.p2c = 0.006 + 0.0069*po4obs ; 

inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan)  = 0  ;
par.nl = 3 ;
par.Lambda = M3d*0 ;
par.nppMSK = M3d*0 ;

% create mask for upper and lower ocean.
UMSK = M3d * 0 ;
DMSK = M3d * 0 ;
par.npp = repmat(npp, [1,1,par.nl]) * 0;
[nx,ny,nz] = size(M3d) ;
for jj = 1 : nx
    for ii = 1 : ny
        tmp = squeeze(M3d(jj, ii, :)) ;
        idp = length(find(tmp) == 1) ;

        if ( idp <= 12 )
            UMSK(jj,ii,1:2) = 1 ;
            DMSK(jj,ii,3:end) = 1;
            par.nppMSK(jj,ii,1) = 2 ;
            par.nppMSK(jj,ii,2) = 2 ;

            par.npp(jj,ii,1)  = (1/2) * npp(jj,ii) / grd.dzt(1) * par.p2c(jj,ii,1);
            par.npp(jj,ii,2)  = (1/2) * npp(jj,ii) / grd.dzt(2) * par.p2c(jj,ii,2);
            par.Lambda(jj,ii,1) = 1./(1e-6+po4obs(jj,ii,1)) ;
            par.Lambda(jj,ii,2) = 1./(1e-6+po4obs(jj,ii,2)) ;
        else 
            UMSK(jj,ii,1:3) = 1 ;
            DMSK(jj,ii,4:end) = 1;
            par.nppMSK(jj,ii,1) = 3 ;
            par.nppMSK(jj,ii,2) = 3 ;
            par.nppMSK(jj,ii,3) = 3 ;

            par.npp(jj,ii,1)  = (1/3) .* npp(jj,ii) / grd.dzt(1) * par.p2c(jj,ii,1);
            par.npp(jj,ii,2)  = (1/3) .* npp(jj,ii) / grd.dzt(2) * par.p2c(jj,ii,2);
            par.npp(jj,ii,3)  = (1/3) .* npp(jj,ii) / grd.dzt(3) * par.p2c(jj,ii,3);
            par.Lambda(jj,ii,1) = 1./(1e-6+po4obs(jj,ii,1)) ; 
            par.Lambda(jj,ii,2) = 1./(1e-6+po4obs(jj,ii,2)) ;
            par.Lambda(jj,ii,3) = 1./(1e-6+po4obs(jj,ii,3)) ;
        end 
    end
end
% UMSK(iwet) = 1 ;
% DMSK(iwet) = 0 ;
UMSK(ilnd) = nan ;
DMSK(ilnd) = nan ; 
par.UM = d0(UMSK(iwet)) ;
par.DM = d0(DMSK(iwet)) ;
par.WM = d0(M3d(iwet))  ;
%---------------------------- end ---------------------------------- 
