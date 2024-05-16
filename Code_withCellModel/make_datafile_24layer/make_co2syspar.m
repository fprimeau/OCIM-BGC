%
%
% build co2syspar for CO2SYS in the OCIM2 24layer model. 
% Input pres: grd.zt(1) Output pres in CO2SYS (in the eqCcycle): 0 m
% 
clc; clear all; close all
addpath('../../DATA/BGC_24layer')
%
load OCIM2_CTL_He.mat output
load TS_WOA_91x180x24.mat
load O2_Nut_WOA_91x180x24.mat Si_obs DIP_obs
grd = output.grid;
M3d = output.M3d;
presin = grd.ZT3d(:,:,1);  %should this be pressure instead of depth?
presout = zeros(size(M3d(:,:,1)));
%
tmp = M3d;
tmp(:,:,2:end) = 0;
isrf = find(tmp(:));
%
co2syspar.presin  = presin(isrf);
co2syspar.presout = presout(isrf);
co2syspar.si      = Si_obs(isrf);
co2syspar.po4     = DIP_obs(isrf);
co2syspar.temp    = tempobs(isrf);
co2syspar.salt    = salobs(isrf);


%Save pme file in DATA directory
fileName  = 'co2syspar_91x180x24.mat'
directory = '../../DATA/BGC_24layer/'
filePath  = fullfile(directory, fileName);
save(filePath, 'co2syspar');
