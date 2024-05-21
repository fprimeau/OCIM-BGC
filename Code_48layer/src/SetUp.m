global GC GO iter
iter = 0       ;
spd  = 24*60^2 ;
spa  = 365*spd ;
on   = true    ;
off  = false   ;
% addpath according to opterating system
addpath('../../DATA/BGC_48layer/')

if GridVer == 91 
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
      case 'L'
        TRdivVer = 'CTL_He_48layer'     ;
    end 
end 

fprintf('Transport version: % s \n', TRdivVer)
if par.Cmodel == on
    fprintf('-------- C model is on -------- \n')
    fprintf('DOP scaling factor is %2.2e \n', par.dopscale)
    fprintf('DIP scaling factor is %2.2e \n', par.dipscale)
    fprintf('DIC scaling factor is %2.2e \n', par.dicscale)
    fprintf('DOC scaling factor is %2.2e \n', par.docscale)
    fprintf('ALK scaling factor is %2.2e \n', par.alkscale)
    fprintf('O2 scaling factor is %2.2e \n', par.o2scale)
end
if par.Omodel == on
    fprintf('-------- O model is on -------- \n')
end 
if par.Simodel == on
    fprintf('-------- Si model is on -------- \n')
end 
fprintf('\n')

%
if GridVer == 91
    OperName = sprintf('OCIM2_%s',TRdivVer);
    load(OperName,'output', 'TR') ;          % TR: convergence & yr-1
    %
    fname = 'biopump_model_output_Nowicki.nc';
    load NPP_climatology.mat NPP_48layer % mmol/m2/yr
    npp = NPP_48layer(:,:,1); % 1:CbPM_SW; 2:CAFE_SW; 3:CbPM_MODIS; 4:CAFE_MODIS
    %
    load MSKS_48layer.mat
    load TS_WOA_91x180x48.mat   % WOA temperature & salinity
    load O2_Nut_WOA_91x180x48.mat   % WOA O2 Si DIN DIP observations
    load GLODAPv2_2023_91x180x48.mat alkraw o2raw po4raw %GLODAP alk, o2, DIP(po4), Si
    load docraw_91x180x48.mat % DOC observation data
    load dopraw_91x180x48.mat % DOP observation data
    load pme_91x180_48layer.mat
    load dic_initial_91x180x48.mat      %gc12new 이용 interpolation. 추후 다시 re-update and opimization 해야함
    load splco2_mod_monthly % monthly CO2 data ---> 1765 ~ 1999까지.
    load co2syspar_48layer.mat
    load Initial_COfield_91x180x48.mat % initial guess for C and O fileds  
    
    M3d = output.M3d;
    grd = output.grid;
    TR  = TR/spa;           % yr-1   -------> s-1
end
   
%---------------------- constants -------------------
par.spd = spd ;
par.spa = spa ;
[par.kw, par.P] = kw(M3d, grd) ;
par.rho   = 1025         ; % seawater density; pme파일과 동일하게
permil    = par.rho*1e-3 ; % from umol/kg to mmol/m3; 
par.permil = permil      ;
% transiant CO2 concentraion ;
par.year      = splco2_mod(:,1) ;     
par.pco2_air  = splco2_mod(:,2) ;     
par.pco2atm   = splco2_mod(1,2) ;     % Year:1765 pCO2:277.9541 --> needed transient run
par.co2syspar = co2syspar       ;

%------------------- model grid info ----------------
ilnd = find(M3d(:) == 0)    ;
iwet = find(M3d(:))         ;
nwet = length(iwet)         ;
dAt  = output.grid.dAt      ;
dVt  = output.grid.dVt      ;
ARC  = MSKS.ARC_48layer     ;
MED  = MSKS.MED_48layer     ;
PAC  = MSKS.PAC_48layer     ;
ATL  = MSKS.ATL_48layer     ;
IND  = MSKS.IND_48layer     ;
iarc = find(ARC(:))         ;
imed = find(MED(:))         ;
ipac = find(PAC(:))         ;
iatl = find(ATL(:))         ;
iind = find(IND(:))         ;
%
par.dAt   = dAt     ;
par.dVt   = dVt     ;
par.M3d   = M3d     ;
par.iwet  = iwet    ;
par.nwet  = nwet    ;
par.TRdiv = -TR     ;             % Convergence ----> Divergence
par.grd   = grd     ;
par.MSKS  = MSKS    ;
par.I = speye(nwet) ;

%------------------------  data  -----------------------
dicraw = dic_initial;           % ---> Transient run 통해 다시 re-optimized 필요
% get rid of arctice o2 observations
% docraw(iarc)  = nan ;            % 기존 weilei 모델에서는 arc med자료 뺐음. 이번에는 다 포함시켜 보자.
% docraw(imed)  = nan ;
% alkraw(iarc)  = nan ;
% alkraw(imed)  = nan ; 
% dicraw(iarc)  = nan ;
% dicraw(imed)  = nan ;
% po4raw(iarc)  = nan ;
% po4raw(imed)  = nan ; 
% o2raw(iarc)   = nan ;
% o2raw(imed)   = nan ;


% remove outlayers of DOC obs    ----> 여기서는 모두 포함해보자.
% docraw = rmOutliers(docraw, M3d) ;
% po4raw = rmOutliers(po4raw, M3d) ;
% dicraw = rmOutliers(dicraw, M3d) ;
% alkraw = rmOutliers(alkraw, M3d) ;
% o2raw  = rmOutliers(o2raw, M3d)  ;

tempobs(tempobs(:)<-2.0) = -2.0   ;
DIP_obs(DIP_obs(:)<0.05) = 0.05   ;
po4raw(po4raw(:)<0.05)   = nan    ; 

%for ji = 1:24
%    p2d = po4obs(:,:,ji);
%    po4obs(:,:,ji) = smoothit(grd,M3d,p2d,3,1e5);   % 이게 필요한 이유...?
%end                                       ----> NPP field에 영향을 주려나?
 
par.Temp     = tempobs       ;
par.Salt     = salobs        ;
par.DSi      = Si_obs        ;
par.po4obs   = DIP_obs       ;
par.o2raw    = o2raw         ;
par.o2obs    = O2_obs        ;     %WOA O2가 어디에 쓰이는지는 모르겠음.
par.po4raw   = po4raw        ;
par.DOCobs   = docraw        ;
par.dicraw   = dicraw        ; 
par.alkraw   = alkraw*permil ;
par.dopraw   = dopraw - 0.05 ;    %less refractory DOP  ---> 0.05를 왜 빼주지..?
par.docraw   = docraw ;
%-------------------- prepare virtual flux -------------------
% PME part;
% [modT,modS,pme] = pme(par) ; 
par.pme = pme_new(iwet) ;   
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
%for ji = 1:24
%    t2d = par.Temp(:,:,ji); 
%    par.Temp(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);        % smoothit ---> inpaint_nan하면 안돼?
%end                                                        % 일단 빼고 해보기.

vT = par.Temp(iwet) ;                                     
% add + 1.0 to prevent from getting infinit kP or kC 
Tz = (vT - min(vT) + 1.0)./(max(vT) - min(vT)) ;
par.Tz = Tz      ;                                         % eqOcycle에서 수정해야함.
par.vT = vT      ;
Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par.aveT   = nanmean(Tz3d(:,:,1:3),3) ;                    % tsnanmean하고 큰 차이?
                                                             % gp에서 run할때는
                                                             % nanmean으로 수정
                                                             % 1-3번째 layer의 T평균? 

% aveT: PFD에 쓰임. vertical하게는 동일.
% vT: eqPcylce, eqCcycle에서 Q10의 tf에 쓰임.
% Tz: aveT에 쓰일때말고는 모르겠음. 1e-8을 왜하는지? ---> eqOcycle에서는 다시 10^8 곱해서 사용.
                                                         

%-------------------- correct o2 concentration --------------------
%----------> 내가 받은 O2_obs WOA는 이미 correction 되어 있음. 여긴 일단 넘어가기.
% o2obs_c = M3d*0;
% o2obs(iwet) = o2obs(iwet).*44.661;        % convert unit form [ml/l] to [umol/l].
% o2obs_c(iwet) = o2obs(iwet).*1.009-2.523; % o2 correction based on Bianchi et al.(2012) [umol/l] .
% ineg = find(o2obs_c(:)<0);                % find out negative values and set them to zero.
% o2obs_c(ineg) = 0.1 ;
% par.o2obs = o2obs_c ;

%---------------------- prepare for restoring -----------------------
% calculating global mean DIP, ALK, and DSi concentraions for
% restoring 
idip = find(par.po4raw(iwet) > 0) ;
ialk = find(par.alkraw(iwet) > 0) ;

par.DIPbar = sum(par.po4raw(iwet(idip)).*dVt(iwet(idip)))/sum(dVt(iwet(idip))) ;
par.ALKbar = sum(par.alkraw(iwet(ialk)).*dVt(iwet(ialk)))/sum(dVt(iwet(ialk))) ;


%-------------------- prepare NPP for the model ----------------------
par.p2c = 0.006 + 0.0069*DIP_obs ;         
%par.p2c = (1/117) * M3d ;                % 이건 아마 redfield ratio인듯?
inan = find(isnan(npp(:)) | npp(:) < 0) ;
npp(inan)  = 0  ;
par.nl = 8 ;                              % nl이 뭐지? NPP 발생하는 layer를 의미하나?
par.Lambda = M3d*0 ;
par.nppMSK = M3d*0 ;

% create mask for upper and lower ocean.
UMSK = M3d * 0 ;
DMSK = M3d * 0 ;
par.npp = repmat(npp, [1,1,par.nl]) * 0; %npp 91x180 ----> 91x180xeuphotic zone(~100m)
[x, a, K, C] = V_NPP(par) ;
[nx,ny,nz] = size(M3d) ;
for jj = 1 : nx
    for ii = 1 : ny
        tmp = squeeze(M3d(jj, ii, :)) ;  %----> tmp는 각 91 by 180 field에서의 깊이 (24*1) 벡터
        idp = length(find(tmp) == 1)  ;

        if ( idp <= 49 )                 %----> 필요 없는것 같은데....
            UMSK(jj,ii,1:par.nl) = 1 ;
            DMSK(jj,ii,par.nl+1:end) = 1;
            %
            for kk = 1:par.nl
              par.nppMSK(jj,ii,kk) = par.nl ;
              par.npp(jj,ii,kk)  = x(kk) * npp(jj,ii) / grd.dzt(kk) * par.p2c(jj,ii,kk);
              Cnpp(jj,ii,kk) = x(kk)* npp(jj,ii) / grd.dzt(kk) /spa ;
            end
            %
            for kk = 1:par.nl
              Pnpp(jj,ii,kk)  = par.npp(jj,ii,kk) / spa ;  
              par.Lamda(jj,ii,kk) = 1./(1e-6+DIP_obs(jj,ii,kk)) ;
            end
        end 
    end
end

[par.Pnpp, par.Cnpp] = deal(M3d*0) ;
par.Pnpp(:,:,1:par.nl) = Pnpp;
par.Cnpp(:,:,1:par.nl) = Cnpp;

UMSK(ilnd) = nan ;
DMSK(ilnd) = nan ; 
par.UM = d0(UMSK(iwet)) ;
par.DM = d0(DMSK(iwet)) ;
par.WM = d0(M3d(iwet))  ;
%---------------------------- end ---------------------------------- 


