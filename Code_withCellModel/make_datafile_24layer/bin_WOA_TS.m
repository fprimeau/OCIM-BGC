% bin_WOA_TS.m
% Script to bin World Ocean Atlas 2018 observations onto the OCIM 24
% layer model grid: Temperature and Salinity climatology
%
clc; close all; clear all;

%---- Add path to Raw Data Files -------
addpath('../utils/')
% addpath('../../Observation_rawdata/WOA')
% addpath('../../DATA/BGC_48layer/')
addpath('../../../../DATASETS/WOA18/')

%--- Add path to OCIM grid file -------
addpath('../../../../DATASETS/OCIM/')

%--- Set Directory to save the resulting regridded .mat file -----
SaveToDir = '../../../../DATASETS/BGC_24layer/';

%https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/bin/woa18.pl
%Download from above link
% an:Objectively analyzed climatology mn:Statistcal mean 
% dd:Number of Observation ma:Seasonal or monthly - annual climatology
% sd:standard deviation from statistcal mean se:standard error from statistcal mean
% oa:Sstatistical mean - objectively nalyzed climatology gp:Number of mean vales within radius of influence
% from WOA2018 Product Documentation

tempWOA      = ncread('woa18_decav_t00_01.nc', 't_an');
salWOA       = ncread('woa18_decav_s00_01.nc', 's_an');
%
tlatitude    = ncread('woa18_decav_t00_01.nc', 'lat'); 
tlongitude   = wrapTo360(ncread('woa18_decav_t00_01.nc', 'lon'));
tdepth       = ncread('woa18_decav_t00_01.nc', 'depth');
%
slatitude    = ncread('woa18_decav_s00_01.nc', 'lat'); 
slongitude   = wrapTo360(ncread('woa18_decav_s00_01.nc', 'lon'));
sdepth       = ncread('woa18_decav_s00_01.nc', 'depth');

temprearrange = permute(tempWOA, [2, 1, 3]); % change the position between lat. and long.
salrearrange  = permute(salWOA,  [2, 1, 3]);

[tlongitude, tlatitude, tdepth] = meshgrid(tlongitude, tlatitude, tdepth);
[slongitude, slatitude, sdepth] = meshgrid(slongitude, slatitude, sdepth);

%--------------------------------------
load OCIM2_Grid.mat grd M3d

[ny,nx,nz] = size(M3d);
[woany,woanx,woanz]= size(temprearrange)

%--------------------------------------------
[sorted_temp_lon , sorted_temp_index]  = sort(tlongitude , 2);
[sorted_sal_lon  , sorted_sal_index ]  = sort(slongitude , 2);

for i = 1:woany;
    for k = 1:woanz;
        temprearrange(i,:,k)  = temprearrange (i,sorted_temp_index(i,:,k) ,k) ;
        salrearrange (i,:,k)  = salrearrange  (i,sorted_sal_index (i,:,k) ,k) ;
    end
end

%
fprintf('Inpaint the nans...')
for k = 1:woanz
   tempobs(:,:,k)  = inpaint_nans(temprearrange(:,:,k));
   salobs (:,:,k)  = inpaint_nans(salrearrange (:,:,k));
end
fprintf('done.\n');
%

%
fprintf('Interpolate in the vertical to the model grid...');
temptmp = zeros(woany,woanx,nz);
saltmp  = zeros(woany,woanx,nz);
for i = 1:woany
  for j = 1:woanx
    temptmp(i,j,:)  = interp1(squeeze(tdepth(1,1,:)) ,squeeze(tempobs(i,j,:)) ,squeeze(grd.ZT3d(1,1,:)), 'linear', 'extrap');
    saltmp (i,j,:)  = interp1(squeeze(sdepth(1,1,:)) ,squeeze(salobs (i,j,:)) ,squeeze(grd.ZT3d(1,1,:)), 'linear', 'extrap');
  end
end
fprintf('done.\n');
%

%
%
fprintf('Interpolate in the horizontal to the model grid...');
tempobs =0*M3d;
salobs =0*M3d;
for k = 1:nz
  tempobs(:,:,k)  = interp2(sorted_temp_lon(:,:,1) ,tlatitude(:,:,1) ,temptmp(:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
  salobs (:,:,k)  = interp2(sorted_sal_lon (:,:,1) ,slatitude(:,:,1) ,saltmp (:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
end
fprintf('done.\n');
%
inan = find(M3d(:)==0);
tempobs(inan) = NaN;
salobs(inan)  = NaN;
%
%-------------Save the files---------
fileName  = 'TS_WOA_91x180x24.mat'
filePath  = fullfile(SaveToDir, fileName);
save(filePath, 'tempobs', 'salobs');
