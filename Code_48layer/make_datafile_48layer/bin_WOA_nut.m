%
%
% Using observed O2, Silicate, Nitrate, Phosphate from WOA2018 data
clc; close all; clear all;

addpath('../src/')
addpath('../../Observation_rawdata/WOA/')
addpath('../../DATA/BGC_48layer/')

%https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/
%Download from above link
% an:Objectively analyzed climatology mn:Statistcal mean 
% dd:Number of Observation ma:Seasonal or monthly - annual climatology
% sd:standard deviation from statistcal mean se:standard error from statistcal mean
% oa:Sstatistical mean - objectively nalyzed climatology gp:Number of mean vales within radius of influence
% from WOA2018 Product Documentation

O2_WOA   = ncread('woa18_all_o00_01.nc', 'o_an');
Si_WOA   = ncread('woa18_all_i00_01.nc', 'i_an');
DIN_WOA  = ncread('woa18_all_n00_01.nc', 'n_an');
DIP_WOA  = ncread('woa18_all_p00_01.nc', 'p_an');

O2_latitude    = ncread('woa18_all_o00_01.nc', 'lat');
Si_latitude    = ncread('woa18_all_i00_01.nc', 'lat');
DIN_latitude   = ncread('woa18_all_n00_01.nc', 'lat');
DIP_latitude   = ncread('woa18_all_p00_01.nc', 'lat');

O2_longitude   = wrapTo360(ncread('woa18_all_o00_01.nc', 'lon'));
Si_longitude   = wrapTo360(ncread('woa18_all_i00_01.nc', 'lon'));
DIN_longitude  = wrapTo360(ncread('woa18_all_n00_01.nc', 'lon'));
DIP_longitude  = wrapTo360(ncread('woa18_all_p00_01.nc', 'lon'));

O2_depth    = ncread('woa18_all_o00_01.nc', 'depth');
Si_depth    = ncread('woa18_all_i00_01.nc', 'depth');
DIN_depth   = ncread('woa18_all_n00_01.nc', 'depth');
DIP_depth   = ncread('woa18_all_p00_01.nc', 'depth');


% change the position between lat. and long.
O2_rearrange   =  permute(O2_WOA,   [2, 1, 3]); 
Si_rearrange   =  permute(Si_WOA,   [2, 1, 3]);
DIN_rearrange  =  permute(DIN_WOA,  [2, 1, 3]); 
DIP_rearrange  =  permute(DIP_WOA,  [2, 1, 3]);

[O2_longitude, O2_latitude, O2_depth]    = meshgrid(O2_longitude,  O2_latitude,  O2_depth) ;
[Si_longitude, Si_latitude, Si_depth]    = meshgrid(Si_longitude,  Si_latitude,  Si_depth) ;
[DIN_longitude, DIN_latitude, DIN_depth] = meshgrid(DIN_longitude, DIN_latitude, DIN_depth);
[DIP_longitude, DIP_latitude, DIP_depth] = meshgrid(DIP_longitude, DIP_latitude, DIP_depth);


%--------------------------------------
load OCIM2_CTL_He_48layer.mat output
M3d = output.M3d  ;
grd = output.grid ;

[ny,nx,nz] = size(M3d);
[woany,woanx,woanz]= size(DIP_rearrange)

%--------------------------------------------
[sorted_O2_lon , sorted_O2_index]  = sort(O2_longitude , 2);
[sorted_Si_lon , sorted_Si_index]  = sort(Si_longitude , 2);
[sorted_DIN_lon, sorted_DIN_index] = sort(DIN_longitude, 2);
[sorted_DIP_lon, sorted_DIP_index] = sort(DIP_longitude, 2);

for i = 1:woany;
    for k = 1:woanz;
        O2_rearrange (i,:,k)  = O2_rearrange (i,sorted_O2_index (i,:,k), k) ;
        Si_rearrange (i,:,k)  = Si_rearrange (i,sorted_Si_index (i,:,k), k) ;
        DIN_rearrange(i,:,k)  = DIN_rearrange(i,sorted_DIN_index(i,:,k), k) ;
        DIP_rearrange(i,:,k)  = DIP_rearrange(i,sorted_DIP_index(i,:,k), k) ;
    end
end

%
fprintf('Inpaint the nans...')
for k = 1:woanz
   O2_obs (:,:,k)  = inpaint_nans(O2_rearrange (:,:,k));
   Si_obs (:,:,k)  = inpaint_nans(Si_rearrange (:,:,k));
   DIN_obs(:,:,k)  = inpaint_nans(DIN_rearrange(:,:,k));
   DIP_obs(:,:,k)  = inpaint_nans(DIP_rearrange(:,:,k));
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
for k = 1:nz
  O2_obs (:,:,k)  = interp2(sorted_O2_lon (:,:,1) ,O2_latitude (:,:,1) ,O2tmp (:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
  Si_obs (:,:,k)  = interp2(sorted_Si_lon (:,:,1) ,Si_latitude (:,:,1) ,Sitmp (:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
  DIN_obs(:,:,k)  = interp2(sorted_DIN_lon(:,:,1) ,DIN_latitude(:,:,1) ,DINtmp(:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
  DIP_obs(:,:,k)  = interp2(sorted_DIP_lon(:,:,1) ,DIP_latitude(:,:,1) ,DIPtmp(:,:,k) ,squeeze(grd.XT3d(:,:,1)),squeeze(grd.YT3d(:,:,1)));
end
fprintf('done.\n');
%
inan = find(M3d(:)==0);
O2_obs(inan) = NaN;
Si_obs(inan) = NaN;
DIN_obs(inan) = NaN;
DIP_obs(inan) = NaN;
%
%-------------Save the files---------
fileName = 'O2_Nut_WOA_91x180x48.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'O2_obs', 'Si_obs', 'DIN_obs', 'DIP_obs');

