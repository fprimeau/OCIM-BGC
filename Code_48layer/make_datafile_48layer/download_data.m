%
%
clc; clear all; close all;
%
% download MODIS data from ocean productivity based on CbPM and CAFE model. 

for yr = 2002:2023
    % wget 
    url  = sprintf('http://orca.science.oregonstate.edu/data/1x2/monthly/cbpm2.modis.r2022/hdf/cbpm.m.%i.tar',yr);
    tdir = '../../Observation_raw data/Satelite/MODIS/CbPM'; 
    filename = sprintf('%s/cbpm.m.%i.tar', tdir, yr);
    websave(filename, url);  
end

% untar
files = dir(fullfile(tdir, '*.tar'));
for i = 1:length(files)
    tarfiles = fullfile(tdir, files(i).name);
    untar(tarfiles, tdir);
end


% unzip
files = dir(fullfile(tdir, '*.gz'));
for i = 1:length(files)
    gzfiles = fullfile(tdir, files(i).name);
    gunzip(gzfiles);
end