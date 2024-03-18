clc; clear all; close all

% 
% Downloaded the data from Holzer, M., DeVries T., and de Laverne, C (2021). 
% Diffusion controls the ventilation of a Pacific Shadow Zone above abyssal overturning, 
% Nature Communications. See opensource in here.
% Dimensions: 91 by 180 by 48
% CFC11 time = 70 years, CFC12 time = 70 years


load OCIM2_48L_base_transport.mat %TR: Convergence & yr-1

% All variables
output.M3d = ncread('OCIM2_48L_base_Data.nc', 'ocnmask'); % only 'wet'box == 1
output.grid.YT3d = ncread('OCIM2_48L_base_Data.nc', 'tlat'); % latitude in each grid (from -89.x to 89.0110)
output.grid.XT3d = ncread('OCIM2_48L_base_Data.nc', 'tlon'); % longtitude in each grid (from 1 to 359)
output.grid.ZT3d = ncread('OCIM2_48L_base_Data.nc', 'tz'); % Depth layer (1st layer: 4.9345 m, 48th layer: 5581.9 m)
output.ulat = ncread('OCIM2_48L_base_Data.nc', 'ulat'); %
output.ulon = ncread('OCIM2_48L_base_Data.nc', 'ulon'); %
output.grid.ZW3d = ncread('OCIM2_48L_base_Data.nc', 'wz'); %
output.grid.dVt = ncread('OCIM2_48L_base_Data.nc', 'vol'); %
output.grid.dAt = ncread('OCIM2_48L_base_Data.nc', 'area'); %
output.kappa_para = ncread('OCIM2_48L_base_Data.nc', 'kappa_para'); %
output.kappa_perp = ncread('OCIM2_48L_base_Data.nc', 'kappa_perp'); %
output.ptemp = ncread('OCIM2_48L_base_Data.nc', 'ptemp'); %
output.salt = ncread('OCIM2_48L_base_Data.nc', 'salt'); %
output.Delta14C = ncread('OCIM2_48L_base_Data.nc', 'Delta14C'); %
output.del3He = ncread('OCIM2_48L_base_Data.nc', 'del3He'); %
output.cfc11 = ncread('OCIM2_48L_base_Data.nc', 'cfc11'); %
output.heatflux = ncread('OCIM2_48L_base_Data.nc', 'heatflux'); %
output.saltflux = ncread('OCIM2_48L_base_Data.nc', 'saltflux'); %
output.mld = ncread('OCIM2_48L_base_Data.nc', 'mld'); %
output.uvel = ncread('OCIM2_48L_base_Data.nc', 'uvel'); %
output.vvel = ncread('OCIM2_48L_base_Data.nc', 'vvel'); %
output.wvel = ncread('OCIM2_48L_base_Data.nc', 'wvel'); %
output.ssh = ncread('OCIM2_48L_base_Data.nc', 'ssh'); %
output.mantle_3he_flux = ncread('OCIM2_48L_base_Data.nc', 'mantle_3he_flux'); %

% make model grid for dzt and dzt for 48layer model
%ZW3d = output.grid.ZW3d;
%zw   = squeeze(ZW3d(1,1,:))';        % Depth coordinate of tracer grid points
%dzt  = nan(1, length(zw));           % Width of tracer grid cell in vertical direction

% make model grid for dzt and dzt for 48layer model
dVt = output.grid.dVt;
dAT = output.grid.dAt;
DZT3d = dVt./dAt;
dzt   = squeeze(DZT3d(1,1,:))';

output.grid.dzt  =  dzt;

%------------------------------------------------
fileName  = 'OCIM2_CTL_He_48layer.mat'
directory = '../DATA/BGC_48layer'
filePath  = fullfile(directory, fileName);
save(filePath, 'output', 'TR');
