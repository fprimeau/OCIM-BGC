clc; clear all; close all
% addpath according to opterating system
spd  = 24*60^2 ;
spa  = 365*spd ;
GridVer  = 91  ;
operator = 'A' ;
par.pscale  = 0.0 ;
par.cscale  = 1.0 ; % factor to weigh DOC in the objective function
TRdivVer = 'CTL_He' ;
if ismac
    addpath('~/Dropbox/myfunc'     )
    addpath('~/OneDrive - xmu.edu.cn/DATA/'    )
    addpath('~/OneDrive - xmu.edu.cn/DATA/OCIM')
    output_dir = sprintf('~/OneDrive/rDOC/MSK%2d/',GridVer); 
elseif isunix
    addpath('~/myfunc/'  )
    addpath('/data2/WangWL_file/DATA/'     )
    addpath('/data2/WangWL_file/DATA/OCIM' )
    output_dir = sprintf('~/rDOC-OP/MSK%2d/',GridVer) ;    
end
VER = strcat(output_dir,TRdivVer);
%
load MLD_CESM_91x180.mat MLD 
% load MLD_MIMOC_91x180.mat MLDmax
load tempobs_91x180x24.mat
load SeaWiFS_Cexp.mat
MLD = MLDmax ;
% Creat output file names based on which model(s) is(are) optimized
base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_CbPM_aveTeu_diffSig_O2C_uniEta');
% base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_MODIS_CbPM_aveTeu_diffSig_O2C_uniEta');
% base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta');
% base_name = strcat(VER,'_PCO_Gamma1to3_POC2DIC_GM15_VGPM_aveTeu_diffSig_O2C_uniEta_noArcMed');
catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
fname = strcat(base_name,catDOC);

par.fname = strcat(fname,'.mat') ; 
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
load(fxhat) ;

OperName = sprintf('OCIM2_%s',TRdivVer);
load(OperName,'output') ;
M3d  = output.M3d ;
grd  = output.grid;
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
%-------------------- normalize temperature --------------------
par.Temp    = tempobs ;
for ji = 1:24
    t2d = par.Temp(:,:,ji);
    par.Temp(:,:,ji) = smoothit(grd,M3d,t2d,3,1e5);
end
vT = par.Temp(iwet) ;
% add + 1.0 to prevent from getting infinit kP or kC
Tz = (vT - min(vT) + 1.0)./(max(vT) - min(vT)) ;
par.Tz = Tz*1e-8 ;

Tz3d = M3d + nan ;
Tz3d(iwet) = Tz  ;
par.aveT   = nanmean(Tz3d(:,:,1:3),3) ;
% ---------------------------------------------------------------
bC = xhat.bC_T*par.aveT + xhat.bC ;
[x,y,z] = size(POCexp) ;
for kk = 1 : x
    for tt = 1 : y
        tmp = squeeze(M3d(kk,tt,:)) ;
        itmp = find(tmp == 1) ;
        if length(itmp) > 12
            POC_adj(kk,tt) = POCexp(kk,tt).*(100/grd.zw(4)).^(-bC(kk,tt)) ;
            DOC_adj(kk,tt) = DOCexp(kk,tt).*(100/grd.zw(4)).^(-bC(kk,tt)) ;
            POC_mld(kk,tt) = POCexp(kk,tt).*(MLD(kk,tt)./grd.zw(4)).(-bC(kk,tt)) ;
            PIC_adj(kk,tt) = PICexp(kk,tt).*exp(-(100 - grd.zw(4))/xhat.d);
        else
            POC_adj(kk,tt) = POCexp(kk,tt).*(100/grd.zw(3)).^(-bC(kk,tt)) ;
            DOC_adj(kk,tt) = DOCexp(kk,tt).*(100/grd.zw(3)).^(-bC(kk,tt)) ;
            POC_mld(kk,tt) = POCexp(kk,tt).*(MLD(kk,tt)./grd.zw(3)).(-bC(kk,tt)) ;
            PIC_adj(kk,tt) = PICexp(kk,tt).*exp(-(100 - grd.zw(3))/xhat.d);

        end 
    end
end
TOC_adj = POC_adj + DOC_adj ;

Lat_HOT  = 22+45/60; Lon_HOT = mod(-158,360);
Lat_BATS = 31+40/60; Lon_BATS = mod((-64-10/60),360);
Lat_OSP  = 50+1/60;  Lon_OSP  = mod((-144-9/60),360);

loc.indx_hot = length(find(grd.xt<Lon_HOT));
loc.indy_hot = length(find(grd.yt<Lat_HOT));

loc.indx_bats = length(find(grd.xt<Lon_BATS));
loc.indy_bats = length(find(grd.yt<Lat_BATS));

loc.indx_osp = length(find(grd.xt<Lon_OSP));
loc.indy_osp = length(find(grd.yt<Lat_OSP));

msk_tropical = M3d(:,:,3)*0;
msk_tropical(length(find(grd.yt<-15)):length(find(grd.yt<15)),:) = 1;
loc.iTropical = find(msk_tropical(:)) ;

junk1 = M3d(:,:,3) * 0 ;
junk1(length(find(grd.yt<-30)):length(find(grd.yt<30)),:) = 1;
msk_subtro = junk1 - msk_tropical ;
loc.iSubtro = find(msk_subtro(:)) ;

junk2  = M3d(:,:,3) * 0 ;
junk2(length(find(grd.yt<-45)):length(find(grd.yt<45)),:) = 1;
msk_subtro_subpo = junk2 - junk1 ;
loc.iSubtro_subpo = find(msk_subtro_subpo(:)) ;

junk3 =  M3d(:,:,3) * 0 ;
junk3(length(find(grd.yt<-60)):length(find(grd.yt<60)),:) = 1;
msk_subpolar = junk3 - junk2 ;
loc.iSubpolar = find(msk_subpolar(:)) ;

% Area of the third layer
loc.Area3 = grd.DXT3d(:,:,3) .* grd.DYT3d(:,:,3) ;
dAt = grd.DXT3d.*grd.DYT3d ;

[nx, ny, nz] = size( TOC_adj ) ;
for kk = 1 : nz
    TOCexp = TOC_adj(:, :, kk) ;
    POCexp = POC_adj(:, :, kk) ;
    DOCexp = DOC_adj(:, :, kk) ;
    PICexp = PIC_adj(:, :, kk) ;
    % sum of TOC export
    tem_Cexp = TOCexp.*dAt(:,:,3) ;
    ERR.sTOC(kk) = nansum(tem_Cexp(:))*12*1e-18 ;

    % sum of POC export
    tmp_POCexp = POCexp.*dAt(:,:,3) ;
    ERR.sPOC(kk)   = nansum(tmp_POCexp(:))*12*1e-18 ;

    % sum of PIC export
    tmp_PICexp = PICexp.*dAt(:,:,3) ;
    ERR.sPIC(kk)   = nansum(tmp_PICexp(:))*12*1e-18 ;

    % sum of DOC export
    tmp_DOCexp = DOCexp.*dAt(:,:,3) ;
    ERR.sDOC(kk) = nansum(tmp_DOCexp(:))*12*1e-18 ;

    DOCexp(DOCexp(:)<0) = 0 ;
    EXP = getReg( par, loc, TOCexp, DOCexp ) ;
    ERR.TC_HOT(kk)  = EXP.TC_HOT ;
    ERR.TC_BATS(kk) = EXP.TC_BATS ;
    ERR.TC_OSP(kk)  = EXP.TC_OSP ;
    ERR.mean_TOC_tropical(kk) = EXP.mean_TOC_tropical ;
    ERR.mean_TOC_subtro(kk) = EXP.mean_TOC_subtro ;
    ERR.mean_TOC_subtro_subpo(kk) = EXP.mean_TOC_subtro_subpo ;
    ERR.mean_TOC_subpolar(kk) = EXP.mean_TOC_subpolar ;
    ERR.mean_D2T_tropical(kk) = EXP.mean_D2T_tropical ;
    ERR.mean_D2T_subtro(kk) = EXP.mean_D2T_subtro ;
    ERR.mean_D2T_subtro_subpo(kk) = EXP.mean_D2T_subtro_subpo ;
    ERR.mean_D2T_subpolar(kk) = EXP.mean_D2T_subpolar ;

end
% save Cexp_adj POC_adj TOC_adj DOC_adj PIC_adj ERR
function EXP = getReg( par, loc, TOCexp,  DOCexp)
% find ANCP at specific location and convert unit from
% mg/m2/day to mol/m2/year;
    indx_hot  = loc.indx_hot  ;
    indy_hot  = loc.indy_hot  ;
    indx_osp  = loc.indx_osp  ;
    indy_osp  = loc.indy_osp  ;
    indx_bats = loc.indx_bats ;
    indy_bats = loc.indy_bats ;
    iTropical = loc.iTropical ;
    iSubtro   = loc.iSubtro   ;
    iSubtro_subpo = loc.iSubtro_subpo ;
    iSubpolar = loc.iSubpolar ;
    Area3 = loc.Area3 ; 
    EXP.TC_HOT  = TOCexp(indy_hot,indx_hot)/1000;
    EXP.TC_BATS = TOCexp(indy_bats,indx_bats)/1000;
    EXP.TC_OSP  = TOCexp(indy_osp,indx_osp)/1000;
    
    % calculate mean DOC export and DOC./TOC export ratio.
    % units mol/m2/year;
    TOCexp_tropical = nansum(TOCexp(iTropical).*Area3(iTropical)) / nansum(Area3(iTropical)) ;
    TOCexp_subtro = nansum(TOCexp(iSubtro) .* Area3(iSubtro)) / nansum(Area3(iSubtro)) ;
    TOCexp_subtro_subpo = nansum(TOCexp(iSubtro_subpo).*Area3(iSubtro_subpo)) / nansum(Area3(iSubtro_subpo)) ;
    TOCexp_subpolar = nansum(TOCexp(iSubpolar).*Area3(iSubpolar)) / nansum(Area3(iSubpolar)) ;
    
    EXP.mean_TOC_tropical = TOCexp_tropical/1000;
    EXP.mean_TOC_subtro = TOCexp_subtro/1000;
    EXP.mean_TOC_subtro_subpo = TOCexp_subtro_subpo/1000;
    EXP.mean_TOC_subpolar = TOCexp_subpolar/1000;
    
    % calculate DOC to TOC export ratio for the four biomes.
    D2T = DOCexp./TOCexp;
    D2T(isinf(D2T)) = nan ;
    D2T_tropical = D2T(iTropical) ;
    D2T_subtro = D2T(iSubtro) ;
    D2T_subtro_subpo = D2T(iSubtro_subpo) ;
    D2T_subpolar = D2T(iSubpolar) ;
    
    EXP.mean_D2T_tropical = nansum(D2T_tropical.*Area3(iTropical)) / nansum(Area3(iTropical)) ;
    EXP.mean_D2T_subtro   = nansum(D2T_subtro.*Area3(iSubtro)) / nansum(Area3(iSubtro)) ;
    EXP.mean_D2T_subtro_subpo = nansum(D2T_subtro_subpo.*Area3(iSubtro_subpo)) / nansum(Area3(iSubtro_subpo));
    EXP.mean_D2T_subpolar = nansum(D2T_subpolar.*Area3(iSubpolar)) / nansum(Area3(iSubpolar));

end

