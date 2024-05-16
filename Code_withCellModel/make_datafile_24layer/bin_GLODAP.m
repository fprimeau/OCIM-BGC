% Editing the bin_glodap.m in Myfunction directory
% Using observed DIP, DIC, O2, ALK from GLODAPv2.2023 to OCIM2 48layer grid

clc; close all; clear all;

addpath('../src')
addpath('../../Observation_rawdata/GLODAP')
addpath('../../DATA/BGC_48layer/')
load GLODAPv2_2023_Merged_Master_File.mat 
%https://glodap.info/index.php/merged-and-adjusted-data-product-v2-2023/
%Download GLODAPv2.2023 data from above link


load OCIM2_CTL_He_48layer.mat
%Holzer, M., DeVries T., and de Laverne, C (2021). 
%Diffusion controls the ventilation of a Pacific Shadow Zone above abyssal overturning
%Make new mfile including Transport operator and Grid intormation
%>>> in /Transport_Operator/48layer/ directory 

M3d = output.M3d  ;
grd = output.grid ; % grid >> XT3d, YT3d, ZT3d. 3d array (91 by 180 by 48)


x = wrapTo360(G2longitude) ;
y = G2latitude ;
z = G2depth    ;

%These data are believed to be accurate to 0.005 in salinity, 
%1% in oxygen, 2% in nitrate, 2% in silicate, 2% in phosphate, 
%4 µmol kg-1 in TCO2, 4 µmol kg-1 in TAlk, and for the halogenated transient tracers and SF6: 5%.
%excluding Nana data from GLODAP and use bin3d to bin to the OCIM grid.

io2data = find(~isnan(G2oxygen));
o2lon   = x(io2data);
o2lat   = y(io2data);
o2dep   = z(io2data);
o2raw   = G2oxygen(io2data);
o2err   = o2raw*0.01;

ipo4data = find(~isnan(G2phosphate));
po4lon   = x(ipo4data);
po4lat   = y(ipo4data);
po4dep   = z(ipo4data);
po4raw   = G2phosphate(ipo4data);
po4err   = po4raw*0.02;

idicdata = find(~isnan(G2tco2));
diclon   = x(idicdata);
diclat   = y(idicdata);
dicdep   = z(idicdata);
dicraw   = G2tco2(idicdata);      %It(Observational result) includes anthropogenic C.
dicerr   = dicraw*0.05;

ialkdata = find(~isnan(G2talk));
alklon   = x(ialkdata);
alklat   = y(ialkdata);
alkdep   = z(ialkdata);
alkraw   = G2talk(ialkdata);     
alkerr   = alkraw*0.05;
%TALKraw(TALKraw<1000) = nan;  % #19 of total 1402829 data = 0.001%.

isidata = find(~isnan(G2silicate));
silon   = x(isidata);
silat   = y(isidata);
sidep   = z(isidata);
siraw   = G2silicate(isidata);     
sierr   = siraw*0.05;


[o2raw,o2err]     = bin3d(o2lon,o2lat,o2dep,o2raw,o2err,grd.XT3d,grd.YT3d,grd.ZT3d);
[po4raw,po4err]   = bin3d(po4lon,po4lat,po4dep,po4raw,po4err,grd.XT3d,grd.YT3d,grd.ZT3d);
[dicraw,DICerr]   = bin3d(diclon,diclat,dicdep,dicraw,dicerr,grd.XT3d,grd.YT3d,grd.ZT3d);
[alkraw,TALKerr]  = bin3d(alklon,alklat,alkdep,alkraw,alkerr,grd.XT3d,grd.YT3d,grd.ZT3d);
[siraw,sierr]     = bin3d(silon,silat,sidep,siraw,sierr,grd.XT3d,grd.YT3d,grd.ZT3d);

o2raw (o2raw(:) ==-9.999)=NaN;
po4raw(po4raw(:)==-9.999)=NaN;
dicraw(dicraw(:)==-9.999)=NaN;
alkraw(alkraw(:)==-9.999)=NaN;
siraw (siraw(:) ==-9.999)=NaN;

o2raw (M3d(:)==0)=NaN;
po4raw(M3d(:)==0)=NaN;
dicraw(M3d(:)==0)=NaN;
alkraw(M3d(:)==0)=NaN;
siraw (M3d(:)==0)=NaN;


%DATA file deirectory에 저장하기
fileName = 'GLODAPv2_2023_91x180x48.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'o2raw', 'po4raw', 'dicraw', 'alkraw', 'siraw');


