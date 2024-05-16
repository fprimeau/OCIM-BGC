% Editing the bin_glodap.m in Myfunction directory
% Using observed DON from Liang et al. 2022 to OCIM2 48layer grid
% Ref: A global ocean dissolved organic phosphorus concentration database
% (DOPv2021)
% Data: https://www.bco-dmo.org/dataset/855139
clc; close all; clear all;

addpath('../src')
addpath('../../Observation_rawdata')
addpath('../../DATA/BGC_48layer/')

dataTable = readtable('dopv2021v4-1_raw.csv');
dataTable(1,:) = [];
longtitude = table2array(dataTable(:,7)); %o-180ow to 180oE. Thus, we need wrapTo360
latitude   = table2array(dataTable(:,6)); %-90 to 90
depth      = table2array(dataTable(:,8)); 
dop        = table2array(dataTable(:,13)); %umol kg-1

load OCIM2_CTL_He_48layer.mat
M3d = output.M3d  ;
grd = output.grid ; % grid >> XT3d, YT3d, ZT3d. 3d array (91 by 180 by 48)

x = wrapTo360(longtitude) ;
y = latitude ;
z = depth    ;

dop(dop <= 0) = NaN;
idopdata = find(~isnan(dop));
doplon   = x(idopdata);
doplat   = y(idopdata);
dopdep   = z(idopdata);
dopraw   = dop(idopdata);
doperr   = dopraw*0.05;      % No reports on analytical uncertainty in REF. Shoule be checked.


[dopraw,doperr]     = bin3d(doplon,doplat,dopdep,dopraw,doperr,grd.XT3d,grd.YT3d,grd.ZT3d);

dopraw(dopraw(:)==-9.999)=NaN;
dopraw(M3d(:)==0)=NaN;

%DATA file deirectory에 저장하기
fileName = 'dopraw_91x180x48.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'dopraw');