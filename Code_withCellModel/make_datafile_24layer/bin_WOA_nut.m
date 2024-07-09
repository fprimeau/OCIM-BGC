%% bin_WOA_nut.m
% Script to bin World Ocean Atlas 2018 observations onto the OCIM2 24
% layer model grid: Oxygen, Phosphate, Nitrate, Silicate climatology
%
%
% Using observed O2, Silicate, Nitrate, Phosphate from WOA2018 data
clc; close all; clear all;

%addpath('../src/')
addpath('../utils/')
%---- Add path to Raw Data Files -------
addpath('../../../../DATASETS/WOA18/')

%--- Add path to OCIM grid file -------
addpath('../../../../DATASETS/OCIM/')

%--- Set Directory to save the resulting regridded .mat file -----
SaveToDir = '../../../../DATASETS/BGC_24layer/';

%https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/
%Download from above link
% an:Objectively analyzed climatology mn:Statistcal mean 
% dd:Number of Observation ma:Seasonal or monthly - annual climatology
% sd:standard deviation from statistcal mean se:standard error from statistcal mean
% oa:statistical mean - objectively analyzed climatology gp:Number of mean values within radius of influence
% from WOA2018 Product Documentation

O2_WOA   = ncread('woa18_all_o00_01.nc', 'o_an');   % unit: [umol/kg]
Si_WOA   = ncread('woa18_all_i00_01.nc', 'i_an');   % unit: [umol/kg]
DIN_WOA  = ncread('woa18_all_n00_01.nc', 'n_an');   % unit: [umol/kg]
DIP_WOA  = ncread('woa18_all_p00_01.nc', 'p_an');   % unit: [umol/kg]
sigma_WOA = ncread('woa18_decav_I00_01.nc','I_an'); % unit: ['kilograms_per_cubic_meter']
% convert sea_water_sigma to in situ density kg/m^3 
%     sea_water_sigma = [(density(t,s,z) - 1000] kg/m3; 
%     range(sea_water_sigma): [0.9, 53.1]
rho_WOA   = sigma_WOA + 1000;   

% Convert nutrient units from umol/kg to mmol/m3; 
DIP_WOA = DIP_WOA.*rho_WOA.*1e-3 ;
DIN_WOA = DIN_WOA.*rho_WOA.*1e-3 ;
O2_WOA  = O2_WOA.*rho_WOA.*1e-3  ;
Si_WOA  = Si_WOA.*rho_WOA.*1e-3  ;

O2_latitude    = ncread('woa18_all_o00_01.nc', 'lat');
Si_latitude    = ncread('woa18_all_i00_01.nc', 'lat');
DIN_latitude   = ncread('woa18_all_n00_01.nc', 'lat');
DIP_latitude   = ncread('woa18_all_p00_01.nc', 'lat');
rho_latitude   = ncread('woa18_decav_I00_01.nc', 'lat');

O2_longitude   = wrapTo360(ncread('woa18_all_o00_01.nc', 'lon'));
Si_longitude   = wrapTo360(ncread('woa18_all_i00_01.nc', 'lon'));
DIN_longitude  = wrapTo360(ncread('woa18_all_n00_01.nc', 'lon'));
DIP_longitude  = wrapTo360(ncread('woa18_all_p00_01.nc', 'lon'));
rho_longitude  = wrapTo360(ncread('woa18_decav_I00_01.nc', 'lon'));

O2_depth    = ncread('woa18_all_o00_01.nc', 'depth');
Si_depth    = ncread('woa18_all_i00_01.nc', 'depth');
DIN_depth   = ncread('woa18_all_n00_01.nc', 'depth');
DIP_depth   = ncread('woa18_all_p00_01.nc', 'depth');
rho_depth   = ncread('woa18_decav_I00_01.nc', 'depth');


% change the position between lat. and long.
O2_rearrange   =  permute(O2_WOA,   [2, 1, 3]); 
Si_rearrange   =  permute(Si_WOA,   [2, 1, 3]);
DIN_rearrange  =  permute(DIN_WOA,  [2, 1, 3]); 
DIP_rearrange  =  permute(DIP_WOA,  [2, 1, 3]);
rho_rearrange  =  permute(rho_WOA,   [2, 1, 3]); 

[O2_longitude, O2_latitude, O2_depth]    = meshgrid(O2_longitude,  O2_latitude,  O2_depth) ;
[Si_longitude, Si_latitude, Si_depth]    = meshgrid(Si_longitude,  Si_latitude,  Si_depth) ;
[DIN_longitude, DIN_latitude, DIN_depth] = meshgrid(DIN_longitude, DIN_latitude, DIN_depth);
[DIP_longitude, DIP_latitude, DIP_depth] = meshgrid(DIP_longitude, DIP_latitude, DIP_depth);
[rho_longitude, rho_latitude, rho_depth] = meshgrid(rho_longitude,  rho_latitude,  rho_depth) ;

%--------------------------------------
% load OCIM2_CTL_He.mat output
% M3d = output.M3d  ;
% grd = output.grid ;
load OCIM2_Grid.mat grd M3d

[ny,nx,nz] = size(M3d);
[woany,woanx,woanz]= size(DIP_rearrange)

%--------------------------------------------
[sorted_O2_lon , sorted_O2_index]  = sort(O2_longitude , 2);
[sorted_Si_lon , sorted_Si_index]  = sort(Si_longitude , 2);
[sorted_DIN_lon, sorted_DIN_index] = sort(DIN_longitude, 2);
[sorted_DIP_lon, sorted_DIP_index] = sort(DIP_longitude, 2);
[sorted_rho_lon, sorted_rho_index] = sort(rho_longitude, 2);

for i = 1:woany;
    for k = 1:woanz;
        O2_rearrange (i,:,k)  = O2_rearrange (i,sorted_O2_index (i,:,k), k) ;
        Si_rearrange (i,:,k)  = Si_rearrange (i,sorted_Si_index (i,:,k), k) ;
        DIN_rearrange(i,:,k)  = DIN_rearrange(i,sorted_DIN_index(i,:,k), k) ;
        DIP_rearrange(i,:,k)  = DIP_rearrange(i,sorted_DIP_index(i,:,k), k) ;
        rho_rearrange(i,:,k)  = rho_rearrange(i,sorted_rho_index(i,:,k), k) ;
    end
end

%
fprintf('Inpaint the nans...')
for k = 1:woanz
   O2_obs (:,:,k)  = inpaint_nans(O2_rearrange (:,:,k));
   Si_obs (:,:,k)  = inpaint_nans(Si_rearrange (:,:,k));
   DIN_obs(:,:,k)  = inpaint_nans(DIN_rearrange(:,:,k));
   DIP_obs(:,:,k)  = inpaint_nans(DIP_rearrange(:,:,k));
   rho_obs(:,:,k)  = inpaint_nans(rho_rearrange(:,:,k));
end
fprintf('done.\n');
%

%
fprintf('Interpolate in the vertical to the model grid...');
DIPtmp = zeros(woany,woanx,nz);
for i = 1:woany
  for j = 1:woanx
    O2tmp(i,j,:)  = interp1(squeeze(O2_depth(1,1,:)) ,squeeze(O2_obs(i,j,:)) ,squeeze(grd.ZT3d(1,1,:)), 'linear', 'extrap');
    Sitmp(i,j,:)  = interp1(squeeze(Si_depth(1,1,:)) ,squeeze(Si_obs(i,j,:)) ,squeeze(grd.ZT3d(1,1,:)), 'linear', 'extrap');
    DINtmp(i,j,:) = interp1(squeeze(DIN_depth(1,1,:)),squeeze(DIN_obs(i,j,:)),squeeze(grd.ZT3d(1,1,:)), 'linear', 'extrap');
    DIPtmp(i,j,:) = interp1(squeeze(DIP_depth(1,1,:)),squeeze(DIP_obs(i,j,:)),squeeze(grd.ZT3d(1,1,:)), 'linear', 'extrap');
    rhotmp(i,j,:) = interp1(squeeze(rho_depth(1,1,:)),squeeze(rho_obs(i,j,:)),squeeze(grd.ZT3d(1,1,:)), 'linear', 'extrap');
  end
end
fprintf('done.\n');
%

%
fprintf('Interpolate in the horizontal to the model grid...');
O2_obs =0*M3d;
Si_obs =0*M3d;
DIN_obs =0*M3d;
DIP_obs =0*M3d;
rho_obs =0*M3d;
for k = 1:nz
  O2_obs (:,:,k)  = interp2(sorted_O2_lon (:,:,1) ,O2_latitude (:,:,1) ,O2tmp (:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
  Si_obs (:,:,k)  = interp2(sorted_Si_lon (:,:,1) ,Si_latitude (:,:,1) ,Sitmp (:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
  DIN_obs(:,:,k)  = interp2(sorted_DIN_lon(:,:,1) ,DIN_latitude(:,:,1) ,DINtmp(:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
  DIP_obs(:,:,k)  = interp2(sorted_DIP_lon(:,:,1) ,DIP_latitude(:,:,1) ,DIPtmp(:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
  rho_obs(:,:,k)  = interp2(sorted_rho_lon(:,:,1) ,rho_latitude(:,:,1) ,rhotmp(:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
end
fprintf('done.\n');
%
inan = find(M3d(:)==0);
O2_obs(inan) = NaN;
Si_obs(inan) = NaN;
DIN_obs(inan) = NaN;
DIP_obs(inan) = NaN;
rho_obs(inan) = NaN;
%
%-------------Save the files---------
%fileName = 'O2_Nut_WOA_91x180x24.mat'
fileName = 'O2_Nut_WOA_mmolperm3_91x180x24.mat'
filePath = fullfile(SaveToDir, fileName);
save(filePath, 'O2_obs', 'Si_obs', 'DIN_obs', 'DIP_obs','rho_obs');

