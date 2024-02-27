%
%
% build co2syspar for CO2SYS.
%
%

clc; clear all; close all
addpath('../../DATA/BGC_48layer')

%
load OCIM2_CTL_He_48layer.mat output
load TS_WOA_91x180x48.mat
load O2_Nut_WOA_91x180x48.mat Si_obs DIP_obs

grd = output.grid;
M3d = output.M3d;
tmp = M3d;
tmp(:,:,2:end) = 0;
isrf = find(tmp(:));

pres =  squeeze(grd.ZT3d(:,:,1));

co2syspar.pres  = pres(isrf);
co2syspar.si    = Si_obs(isrf);
co2syspar.po4   = DIP_obs(isrf);
co2syspar.temp  = tempobs(isrf);
co2syspar.salt  = salobs(isrf);

%Save pme file in DATA directory
fileName  = 'co2syspar91_48layer.mat'
directory = '../../DATA/BGC_48layer'
filePath  = fullfile(directory, fileName);
save(filePath, 'co2syspar');
