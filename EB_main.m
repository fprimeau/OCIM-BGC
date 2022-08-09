clc; clear alul; close all
global GC
on   = true  ;
off  = false ;
format long
% 
GridVer  = 91  ;
operator = 'A' ;

Gtest = off ;
Htest = off ;
par.optim   = off ; 
par.Cmodel  = on ; 
par.Omodel  = on ; 
par.Simodel = off ;
par.LoadOpt = on ; % if load optimial par. 
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigma = on ; 
par.opt_kP_T  = off ;
par.opt_kdP   = on ;
par.opt_bP_T  = off ; 
par.opt_bP    = on ;
par.opt_alpha = on ;
par.opt_beta  = on ;
% C model parameters
par.opt_bC_T  = off ;
par.opt_bC    = on ; 
par.opt_d     = on ;
par.opt_kC_T  = off ;
par.opt_kdC   = on ; 
par.opt_R_Si  = off ; 
par.opt_rR    = on ; 
par.opt_cc    = on ;
par.opt_dd    = on ;
% O model parameters
par.opt_O2C_T = off ;
par.opt_rO2C  = on ;
par.opt_O2P_T = off ; 
par.opt_rO2P  = on ; 
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;
%
%-------------load data and set up parameters---------------------
SetUp ;

% save results 
% ATTENTION: Change this direcrtory to where you wanna save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer); 
elseif isunix
    output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/Cexp/']);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
                        % 'MSK%2d/PME4DICALK/'],GridVer);
                        % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
    % 'MSK%2d/'],GridVer);
end
VER = strcat(output_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_P');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
        base_name = strcat(VER,'_PCv1');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
        base_name = strcat(VER,'_PCOv1');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
        base_name = strcat(VER,'_PCSi');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
        base_name = strcat(VER,'_PCOSi');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    end
end
par.fname = strcat(fname,'.mat') ; 
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ; 
load(par.fname) ;
load(par.fxhat) ;
ferror = strcat(output_dir,'/MT_ERR.mat') ;

%---------------- inital guesses on C and O ---------------
if par.Cmodel == on 
    DIC = data.DIC - par.dicant ;
    
    GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); ...
           data.PIC(iwet); data.ALK(iwet)];
end 
if par.Omodel == on 
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end 

%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values 
par = SetPar(par) ;
% pack parameters into an array, assign them corresponding indices.
par = PackPar(par) ;

idip = find(par.po4raw(iwet)>0)  ;
isil = find(par.sio4raw(iwet)>0) ;
idic = find(par.dicraw(iwet)>0)  ;
idoc = find(par.docraw(iwet)>0)  ;
io2  = find(par.o2raw(iwet)>0)   ;

ndip = length(idip) ;
nsil = length(isil) ;
ndic = length(idic) ;
ndoc = length(idoc) ;
no2  = length(io2)  ;

if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
    sig = 2*xhat.f/ndip ;
    HH  = xhat.fxx/sig  ;
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
    if fscale ~= 0
        sig = (2*xhat.f)/(ndip+ndic+ndoc);
        HH  = xhat.fxx/sig  ;
    else
        sig = (2*xhat.f)/(ndip+ndic);
        HH  = xhat.fxx/sig  ;
    end 
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
    sig = (2*xhat.f)/(ndip+ndic+no2);
    HH  = xhat.fxx/sig  ;
elseif (par.Cmodel == off & par.Omodel == off & par.Simodel == on)
    sig = (2*xhat.f)/(ndip+nsil);
    HH  = xhat.fxx/sig  ; 
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
    sig = (2*xhat.f)/(ndip+ndic+no2+nsil) ;
    HH  = xhat.fxx/sig  ;
end
Hessian = inv(HH);
Hessian = (Hessian+Hessian.')/2;

iter = 1000;
TOC_exp = zeros(91,180,iter);
POC_exp = zeros(91,180,iter);
DOC_exp = zeros(91,180,iter);
for ji = 3*iter+1 : 4*iter 
    if mod(ji,100) == 0
        fprintf('current iteration is %d ....:-) \n',ji);
    end
    x0 = mvnrnd(par.p0, Hessian);
    if par.opt_sigma
        sigma = exp(x0(par.pindx.lsigma));
        if sigma >= 1
            x0(par.pindx.lsigma) = par.p0(par.pindx.lsigma)  ;
        end 
    end 
    if par.opt_cc
        cc = exp(x0(par.pindx.lcc));
        if cc > 1e-3
            x0(par.pindx.lcc) = par.p0(par.pindx.lcc);
        end 
    end 
    
    [par, P] = eqPcycle(x0, par) ;
    DIP = M3d+nan  ;  DIP(iwet) = P(1+0*nwet:1*nwet) ;
    POP = M3d+nan  ;  POP(iwet) = P(1+1*nwet:2*nwet) ;
    DOP = M3d+nan  ;  DOP(iwet) = P(1+2*nwet:3*nwet) ;
    par.DIP = DIP(iwet) ;

    [par,C] = eqCcycle(x0, par) ;
    DIC = M3d+nan ;  DIC(iwet) = C(0*nwet+1:1*nwet) ;
    POC = M3d+nan ;  POC(iwet) = C(1*nwet+1:2*nwet) ;
    DOC = M3d+nan ;  DOC(iwet) = C(2*nwet+1:3*nwet) ;
    PIC = M3d+nan ;  PIC(iwet) = C(3*nwet+1:4*nwet) ;
    ALK = M3d+nan ;  ALK(iwet) = C(4*nwet+1:5*nwet) ;
    par.POC = POC ;
    %------------------ extract parameters ---------------------------
    % POP disolution constant [s^-1];
    sigma = par.sigma ; 
    % linear parameter of npp to DIP assimilation function. 
    alpha = par.alpha ; 
    % exponential parameter of npp to DIN assimilation function.
    beta  = par.beta ; 

    %------------------ prepare NPP for the model --------------------
    % DIP assimilation
    par.G   = alpha*par.L*DIP(iwet) ; % primary production in P unit.
    
    % -------------------- call export function ----------------------
    EXP  = cexp(par);
    
    % -------------------- put results in a structure ---------------
    TOC_exp(:,:,ji) = EXP.TOCexp ;
    POC_exp(:,:,ji) = EXP.POCexp ;
    DOC_exp(:,:,ji) = EXP.TOCexp-EXP.POCexp ;

    ERR.TC_HOTS = EXP.TC_HOTS ;
    ERR.TC_BATS = EXP.TC_BATS ;
    ERR.TC_OSP  = EXP.TC_OSP  ; 

    ERR.RTOC_tropical(ji) = EXP.mean_TOC_tropical;
    ERR.RTOC_subtro(ji) = EXP.mean_TOC_subtro;
    ERR.RTOC_subtro_subpo(ji) = EXP.mean_TOC_subtro_subpo;
    ERR.RTOC_subpolar(ji) = EXP.mean_TOC_subpolar;
    
    ERR.fD2T_tropical(ji) = EXP.mean_D2T_tropical;
    ERR.fD2T_subtro(ji) = EXP.mean_D2T_subtro;
    ERR.fD2T_subtro_subpo(ji) = EXP.mean_D2T_subtro_subpo;
    ERR.fD2T_subpolar(ji) = EXP.mean_D2T_subpolar;
end
save(ferror, 'ERR', 'TOC_exp', 'POC_exp', 'DOC_exp')
fprintf('---------------END-------------------\n')

function EXP = cexp(par) ;
    nn   = 2 ; 
    sigma   =  par.sigma ;
    kdP     =  par.kdP   ;
    bP      =  par.bP    ;
    alpha   =  par.alpha ;
    beta    =  par.beta  ;
    bC      =  par.bC    ;
    d       =  par.d     ;
    kdC     =  par.kdC   ;
    rR      =  par.rR    ;
    cc      =  par.cc    ;
    dd      =  par.dd    ;
    rO2C    =  par.rO2C  ;
    rO2P    =  par.rO2P  ;

    M3d  = par.M3d  ;
    grd  = par.grd  ;
    iwet = par.iwet ;
    nwet = par.nwet ;
    dAt  = par.dAt  ;
    dVt  = par.dVt  ;
    I    = par.I    ;
    
    PFD = buildPFD(par, 'POC') ; 
    TRD = par.TRdiv ;
    W   = d0(par.dVt(iwet))  ;
    F_diag_p = inv(W)*PFD'*W ;
    T_diag   = inv(W)*TRD'*W ;

    % -------------- C:P uptake ratio --------------------------------
    PO4 = par.po4obs(iwet) ;
    C2P3D = M3d + nan ;
    C2P3D(iwet) = 1./(cc*PO4 + dd) ;
    
    junk = M3d ;
    junk(:,:,1:nn) = 0 ;
    Omega = junk(iwet) ;
    Prod  = par.G.*C2P3D(iwet) ;
    % adjoint method.
    kC    = d0(par.kC_T * par.Tz + kdC) ;
    Jex_C = d0(kC*Prod)*(sigma*I+par.kappa_p*(1-sigma) * ...
                         inv(F_diag_p+par.kappa_p*I))*((T_diag+kC)\Omega); 
    
    C3d = M3d+nan ;
    C3d(iwet) = Jex_C ;
    
    % convert export from mmol C/m^3/s to mg C/m^2/day;
    TOCexp = C3d(:,:,2).*grd.dzt(2)*12*par.spd;
    tem_Cexp = TOCexp.*dAt(:,:,2);
    sTOCexp  = nansum(tem_Cexp(:))*365*1e-18;
    
    % POC export
    [~,Gout] = buildPFD(par,'POC') ;
    w = -Gout.w ;
    POCexp     = par.POC(:,:,2).*w(:,:,2)*12*par.spd ;
    tmp_POCexp = POCexp.*dAt(:,:,2);
    sPOCexp    = nansum(tmp_POCexp(:))*365*1e-18;

    EXP.sTOC   = sTOCexp ;
    EXP.sPOC   = sPOCexp ; 
    EXP.TOCexp = TOCexp  ;
    EXP.POCexp = POCexp  ;
    %--------------- compare to ANCP -----------------------------
    DOCexp = TOCexp - POCexp; DOCexp(DOCexp(:)<0) = 0 ;
    % TOCexp = smoothit(grd,M3d,TOCexp,3,1e5);
    % POCexp = smoothit(grd,M3d,POCexp,3,1e5);
    % DOCexp = smoothit(grd,M3d,DOCexp,3,1e5);
    
    Lat_HOTS = 22+45/60; Lon_HOTS = mod(-158,360);
    Lat_BATS = 31+40/60; Lon_BATS = mod((-64-10/60),360);
    Lat_OSP  = 50+1/60;  Lon_OSP  = mod((-144-9/60),360);
    
    indx_hots = length(find(grd.xt<Lon_HOTS));
    indy_hots = length(find(grd.yt<Lat_HOTS));
    
    indx_bats = length(find(grd.xt<Lon_BATS));
    indy_bats = length(find(grd.yt<Lat_BATS));
    
    indx_osp = length(find(grd.xt<Lon_OSP));
    indy_osp = length(find(grd.yt<Lat_OSP));
    
    % find ANCP at specific location and convert unit from
    % mg/m2/day to mol/m2/year;
    EXP.TC_HOTS = TOCexp(indy_hots,indx_hots)/12/1000*365;
    EXP.TC_BATS = TOCexp(indy_bats,indx_bats)/12/1000*365;
    EXP.TC_OSP  = TOCexp(indy_osp,indx_osp)/12/1000*365;

    % calculate mean DOC export and DOC./TOC export ratio.
    msk_tropical = M3d(:,:,1)*0;
    msk_tropical(length(find(grd.yt<-15)):length(find(grd.yt<15)),:) = 1;
    
    junk1 = M3d(:,:,1)*0;
    junk1(length(find(grd.yt<-30)):length(find(grd.yt<30)),:) = 1;
    msk_subtro = junk1-msk_tropical;
    
    junk2  = M3d(:,:,1)*0;
    junk2(length(find(grd.yt<-45)):length(find(grd.yt<45)),:) = 1;
    msk_subtro_subpo = junk2-junk1;
    
    junk3 =  M3d(:,:,1)*0;
    junk3(length(find(grd.yt<-60)):length(find(grd.yt<60)),:) = 1;
    msk_subpolar = junk3-junk2;
    
    % units mg/m2/day;
    TOCexp_tropical = msk_tropical.*TOCexp;
    TOCexp_subtro = msk_subtro.*TOCexp;
    TOCexp_subtro_subpo = msk_subtro_subpo.*TOCexp;
    TOCexp_subpolar = msk_subpolar.*TOCexp;
    
    TOCexp_tropical(TOCexp_tropical(:)==0) = nan;
    TOCexp_subtro(TOCexp_subtro(:)==0) = nan;
    TOCexp_subtro_subpo(TOCexp_subtro_subpo(:)==0) = nan;
    TOCexp_subpolar(TOCexp_subpolar(:)==0) = nan;
    
    EXP.mean_TOC_tropical = nanmean(TOCexp_tropical(:))/12/1000*365;
    EXP.mean_TOC_subtro = nanmean(TOCexp_subtro(:))/12/1000*365;
    EXP.mean_TOC_subtro_subpo = nanmean(TOCexp_subtro_subpo(:))/12/1000*365;
    EXP.mean_TOC_subpolar = nanmean(TOCexp_subpolar(:))/12/1000*365;
    
    % calculate DOC to TOC export ratio for the four biomes.
    D2T = DOCexp./TOCexp;
    D2T_tropical = (msk_tropical.*D2T);
    D2T_subtro = (msk_subtro.*D2T);
    D2T_subtro_subpo = (msk_subtro_subpo.*D2T);
    D2T_subpolar = (msk_subpolar.*D2T);
    
    D2T_tropical(D2T_tropical(:)==0) = nan;
    D2T_subtro(D2T_subtro(:)==0) = nan;
    D2T_subtro_subpo(D2T_subtro_subpo(:)==0) = nan;
    D2T_subpolar(D2T_subpolar(:)==0) = nan;
    
    EXP.mean_D2T_tropical = nanmean(D2T_tropical(:));
    EXP.mean_D2T_subtro = nanmean(D2T_subtro(:));
    EXP.mean_D2T_subtro_subpo = nanmean(D2T_subtro_subpo(:));
    EXP.mean_D2T_subpolar = nanmean(D2T_subpolar(:));
end 