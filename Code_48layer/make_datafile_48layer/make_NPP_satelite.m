%
% make NPP for OCIM2-48layer horizontal grid using satelite NPP data 
% NPP satelite data are from the following website
% Ocean Productivity website
% http://sites.science.oregonstate.edu/ocean.productivity/custom/php
% downloaded in 02/12/2024
% 1080 by 2160 Monthly HDF files from SeaWiFS and MODIS (CbPM & CAFE model)
%
%
clc; clear all; close all;
addpath(genpath('../../Observation_rawdata/Satelite/'));
addpath('../src');
addpath('../../DATA/BGC_48layer/');
addpath('../../DATA/BGC_2023Nature/');   % for comparing with NOwicki et al. 2002

% load grd for 48layer
load OCIM2_CTL_He_48layer.mat output
grd = output.grid;
msk = output.M3d;
clear output
long_48layer = squeeze(grd.XT3d(:,:,1));
lat_48layer  = squeeze(grd.YT3d(:,:,1));

% For SeaWiFS
% For CbPM,
%
fprintf('%s\n', 'load NPP satelite data...');
tic
yCbPM_SW   = 1997:2010;
NPP_CbPM_SW = NaN(1080,2160,12, length(yCbPM_SW));
%
tic
for i = 1:length(yCbPM_SW);
    flocation_CbPM_SW = sprintf('../../Observation_rawdata/Satelite/SeaWiFS/CbPM/cbpm.s.%d/', yCbPM_SW(i));
    fdir_CbPM_SW      = dir(flocation_CbPM_SW);
    %
    for j = 1:length(fdir_CbPM_SW);
       fnlist_CbPM_SW{j} = fdir_CbPM_SW(j).name;
    end
    fnlist_CbPM_SW(1:2) = [];
    %
    for k = 1 : length(fnlist_CbPM_SW);
        fname   = cell2mat(fnlist_CbPM_SW(k));
        tmpinfo = hdfinfo(fname);
        %
        dsplit  = strsplit(tmpinfo.Attributes(3).Value, '/');
        month   = str2double(dsplit{1});
        %                     
        Data_CbPM_SW             = flipud(hdfread(tmpinfo.SDS));
        NPP_CbPM_SW(:,:,month,i) = Data_CbPM_SW; 
    end
end


%
% For CAFE,
%
yCAFE_SW    = 1997:2009;
NPP_CAFE_SW = NaN(1080,2160,12, length(yCAFE_SW));
%
for i = 1:length(yCAFE_SW);
    flocation_CAFE_SW = sprintf('../../Observation_rawdata/Satelite/SeaWiFS/CAFE/cafe.s.%d/', yCAFE_SW(i));
    fdir_CAFE_SW     = dir(flocation_CAFE_SW);
    
    for j = 1:length(fdir_CAFE_SW);
       fnlist_CAFE_SW{j} = fdir_CAFE_SW(j).name;
    end
    fnlist_CAFE_SW(1:2) = [];
    %
    for k = 1 : length(fnlist_CAFE_SW);
        fname   = cell2mat(fnlist_CAFE_SW(k));
        tmpinfo = hdfinfo(fname);
        %
        dsplit  = strsplit(tmpinfo.Attributes(3).Value, '/');
        month   = str2double(dsplit{1});
        %                     
        Data_CAFE_SW             = flipud(hdfread(tmpinfo.SDS));
        NPP_CAFE_SW(:,:,month,i) =  Data_CAFE_SW; 
    end
end


% For MODIS
% For CbPM,
%
yCbPM_MODIS   = 2002:2023;
NPP_CbPM_MODIS = NaN(1080,2160,12, length(yCbPM_MODIS));
%
for i = 1:length(yCbPM_MODIS);
    flocation_CbPM_MODIS = sprintf('../../Observation_rawdata/Satelite/MODIS/CbPM/cbpm.m.%d/', yCbPM_MODIS(i));
    fdir_CbPM_MODIS      = dir(flocation_CbPM_MODIS);
    %
    for j = 1:length(fdir_CbPM_MODIS);
       fnlist_CbPM_MODIS{j} = fdir_CbPM_MODIS(j).name;
    end
    fnlist_CbPM_MODIS(1:2) = [];
    %
    for k = 1 : length(fnlist_CbPM_MODIS);
        fname   = cell2mat(fnlist_CbPM_MODIS(k));
        tmpinfo = hdfinfo(fname);
        %
        dsplit  = strsplit(tmpinfo.Attributes(3).Value, '/');
        month   = str2double(dsplit{1});
        %                     
        Data_CbPM_MODIS             = flipud(hdfread(tmpinfo.SDS));
        NPP_CbPM_MODIS(:,:,month,i) = Data_CbPM_MODIS; 
    end
end


%
% For CAFE,
%
yCAFE_MODIS    = 2002:2023;
NPP_CAFE_MODIS = NaN(1080,2160,12, length(yCAFE_MODIS));
%
for i = 1:length(yCAFE_MODIS);
    flocation_CAFE_MODIS = sprintf('../../Observation_rawdata/Satelite/MODIS/CAFE/cafe.m.%d/', yCAFE_MODIS(i));
    fdir_CAFE_MODIS      = dir(flocation_CAFE_MODIS);
    
    for j = 1:length(fdir_CAFE_MODIS);
       fnlist_CAFE_MODIS{j} = fdir_CAFE_MODIS(j).name;
    end
    fnlist_CAFE_MODIS(1:2) = [];
    %
    for k = 1 : length(fnlist_CAFE_MODIS);
        fname   = cell2mat(fnlist_CAFE_MODIS(k));
        tmpinfo = hdfinfo(fname);
        %
        dsplit  = strsplit(tmpinfo.Attributes(3).Value, '/');
        month   = str2double(dsplit{1});
        %                     
        Data_CAFE_MODIS             = flipud(hdfread(tmpinfo.SDS));
        NPP_CAFE_MODIS(:,:,month,i) = Data_CAFE_MODIS; 
    end
end
toc

clear Data_CbPM_SW;
clear Data_CAFE_SW;
clear Data_CbPM_MODIS;
clear Data_CAFE_MODIS;

% set latitude and longitude
dx    =   1/6;  dy   =   1/6;
lat   = (-90 +dy/2:dy:90 -dy/2)';
long  = (-180+dx/2:dx:180-dx/2)';
long  = wrapTo360(long);
[long lat] = meshgrid(long, lat);


%---------NO observation:NaN, observation but no value due to ice or angle of satelite:-9999-----%
% filling missing value according to the method from Nowicki et al. (2002)
% the missing values were filled with either 10% of the maximum monthly NPP
% during the year, or the lowest monthly NPP
%
% for CbPM_SW
% ..._fill matrix indicates filling the missing value from satelite NPP data
% NPP_CbPM_SW and NPP_CbPM_SW_fill ---> 1080x2160x12x14 
%
fprintf('%s\n', 'filling missing values...');
tic
NPP_CbPM_SW_fill    = NPP_CbPM_SW;
NPP_CAFE_SW_fill    = NPP_CAFE_SW;
NPP_CbPM_MODIS_fill = NPP_CbPM_MODIS;
NPP_CAFE_MODIS_fill = NPP_CAFE_MODIS;
%
% for CbPM & SeaWiFis
for i = 1:1080;
    for j = 1:2160;
        for l = 1: length(yCbPM_SW);
            %
            if any(~isnan(NPP_CbPM_SW(i,j,:,l)) & NPP_CbPM_SW(i,j,:,l) ~=-9999);
                for k = 1:12
                    monthly_NPP = NPP_CbPM_SW(i,j,k,l);
                    %
                    if monthly_NPP == -9999;
                       yearly_NPP = NPP_CbPM_SW(i,j,:,l);
                       valid_NPP  = yearly_NPP(~isnan(yearly_NPP) & yearly_NPP ~= -9999);
                       max_NPP    = 0.1 * max(valid_NPP);
                       min_NPP    = min(valid_NPP);
                       NPP_CbPM_SW_fill(i,j,k,l) = min(max_NPP, min_NPP);
                    end
                end
            end
        end
    end
end
toc

% for CAFE & SeaWiFis
tic
for i = 1:1080;
    for j = 1:2160;
        for l = 1: length(yCAFE_SW);
            %
            if any(~isnan(NPP_CAFE_SW(i,j,:,l)) & NPP_CAFE_SW(i,j,:,l) ~=-9999);
                for k = 1:12
                    monthly_NPP = NPP_CAFE_SW(i,j,k,l);
                    %
                    if monthly_NPP == -9999;
                       yearly_NPP = NPP_CAFE_SW(i,j,:,l);
                       valid_NPP  = yearly_NPP(~isnan(yearly_NPP) & yearly_NPP ~= -9999);
                       max_NPP    = 0.1 * max(valid_NPP);
                       min_NPP    = min(valid_NPP);
                       NPP_CAFE_SW_fill(i,j,k,l) = min(max_NPP, min_NPP);
                    end
                end
            end
        end
    end
end
toc

% for CbPM & MODIS
tic
for i = 1:1080;
    for j = 1:2160;
        for l = 1: length(yCbPM_MODIS);
            %
            if any(~isnan(NPP_CbPM_MODIS(i,j,:,l)) & NPP_CbPM_MODIS(i,j,:,l) ~=-9999);
                for k = 1:12
                    monthly_NPP = NPP_CbPM_MODIS(i,j,k,l);
                    %
                    if monthly_NPP == -9999;
                       yearly_NPP = NPP_CbPM_MODIS(i,j,:,l);
                       valid_NPP  = yearly_NPP(~isnan(yearly_NPP) & yearly_NPP ~= -9999);
                       max_NPP    = 0.1 * max(valid_NPP);
                       min_NPP    = min(valid_NPP);
                       NPP_CbPM_MODIS_fill(i,j,k,l) = min(max_NPP, min_NPP);
                    end
                end
            end
        end
    end
end
toc
% for CAFE & MODIS
tic
for i = 1:1080;
    for j = 1:2160;
        for l = 1: length(yCAFE_MODIS);
            %
            if any(~isnan(NPP_CAFE_MODIS(i,j,:,l)) & NPP_CAFE_MODIS(i,j,:,l) ~=-9999);
                for k = 1:12
                    monthly_NPP = NPP_CAFE_MODIS(i,j,k,l);
                    %
                    if monthly_NPP == -9999;
                       yearly_NPP = NPP_CAFE_MODIS(i,j,:,l);
                       valid_NPP  = yearly_NPP(~isnan(yearly_NPP) & yearly_NPP ~= -9999);
                       max_NPP    = 0.1 * max(valid_NPP);
                       min_NPP    = min(valid_NPP);
                       NPP_CAFE_MODIS_fill(i,j,k,l) = min(max_NPP, min_NPP);
                    end
                end
            end
        end
    end
end
toc
%
% calculating monthly climatology of NPP using NPP_fill matrix
%
fprintf('%s\n', 'calculating monthly climatology...');
%
NPP_CbPM_SW_fill_monthly      = NaN(1080, 2160, 12);
NPP_CAFE_SW_fill_monthly      = NaN(1080, 2160, 12);
NPP_CbPM_MODIS_fill_monthly   = NaN(1080, 2160, 12);
NPP_CAFE_MODIS_fill_monthly   = NaN(1080, 2160, 12);
%
ndata_CbPM_SW_fill            = NaN(1080, 2160, 12);
ndata_CAFE_SW_fill            = NaN(1080, 2160, 12);
ndata_CbPM_MODIS_fill         = NaN(1080, 2160, 12);
ndata_CAFE_MODIS_fill         = NaN(1080, 2160, 12);
%
mdata_CbPM_SW_fill            = NaN(1080, 2160, 12);
mdata_CAFE_SW_fill            = NaN(1080, 2160, 12);
mdata_CbPM_MODIS_fill         = NaN(1080, 2160, 12);
mdata_CAFE_MODIS_fill         = NaN(1080, 2160, 12);
%
% for CbPM & SeaWiFis
tic
for i = 1:1080;
    for j = 1:2160;
        for k = 1:12;
            tmpdata       = NPP_CbPM_SW_fill(i,j,k,:);
            valid_data    = tmpdata(~isnan(tmpdata) & (tmpdata ~= -9999));
            missing_data  = tmpdata( isnan(tmpdata) | (tmpdata == -9999));
            %
            NPP_CbPM_SW_fill_monthly(i,j,k) = mean(valid_data);
            ndata_CbPM_SW_fill(i,j,k)       = numel(valid_data);
            mdata_CbPM_SW_fill(i,j,k)       = numel(missing_data);
        end
    end
end
toc
% for CAFE & SeaWiFis
tic
for i = 1:1080;
    for j = 1:2160;
        for k = 1:12;
            tmpdata       = NPP_CAFE_SW_fill(i,j,k,:);
            valid_data    = tmpdata(~isnan(tmpdata) & (tmpdata ~= -9999));
            missing_data  = tmpdata( isnan(tmpdata) | (tmpdata == -9999));
            %
            NPP_CAFE_SW_fill_monthly(i,j,k) = mean(valid_data);
            ndata_CAFE_SW_fill(i,j,k)       = numel(valid_data);
            mdata_CAFE_SW_fill(i,j,k)       = numel(missing_data);
        end
    end
end
toc
% for CbPM & MODIS
tic
for i = 1:1080;
    for j = 1:2160;
        for k = 1:12;
            tmpdata       = NPP_CbPM_MODIS_fill(i,j,k,:);
            valid_data    = tmpdata(~isnan(tmpdata) & (tmpdata ~= -9999));
            missing_data  = tmpdata( isnan(tmpdata) | (tmpdata == -9999));
            %
            NPP_CbPM_MODIS_fill_monthly(i,j,k) = mean(valid_data);
            ndata_CbPM_MODIS_fill(i,j,k)       = numel(valid_data);
            mdata_CbPM_MODIS_fill(i,j,k)       = numel(missing_data);
        end
    end
end
toc
% for CAFE & MODIS
tic
for i = 1:1080;
    for j = 1:2160;
        for k = 1:12;
            tmpdata       = NPP_CAFE_MODIS_fill(i,j,k,:);
            valid_data    = tmpdata(~isnan(tmpdata) & (tmpdata ~= -9999));
            missing_data  = tmpdata( isnan(tmpdata) | (tmpdata == -9999));
            %
            NPP_CAFE_MODIS_fill_monthly(i,j,k) = mean(valid_data);
            ndata_CAFE_MODIS_fill(i,j,k)       = numel(valid_data);
            mdata_CAFE_MODIS_fill(i,j,k)       = numel(missing_data);
        end
    end
end
toc

%
NPP_monthly          = zeros(1080, 2160, 12, 4);
NPP_monthly(:,:,:,1) = NPP_CbPM_SW_fill_monthly;
NPP_monthly(:,:,:,2) = NPP_CAFE_SW_fill_monthly;
NPP_monthly(:,:,:,3) = NPP_CbPM_MODIS_fill_monthly;
NPP_monthly(:,:,:,4) = NPP_CAFE_MODIS_fill_monthly;
%
ndata_monthly          = zeros(1080, 2160, 12, 4);
ndata_monthly(:,:,:,1) = ndata_CbPM_SW_fill;
ndata_monthly(:,:,:,2) = ndata_CAFE_SW_fill;
ndata_monthly(:,:,:,3) = ndata_CbPM_MODIS_fill;
ndata_monthly(:,:,:,4) = ndata_CAFE_MODIS_fill;
%
mdata_monthly          = zeros(1080, 2160, 12, 4);
mdata_monthly(:,:,:,1) = ndata_CbPM_SW_fill;
mdata_monthly(:,:,:,2) = ndata_CAFE_SW_fill;
mdata_monthly(:,:,:,3) = ndata_CbPM_MODIS_fill;
mdata_monthly(:,:,:,4) = ndata_CAFE_MODIS_fill;
%
monthly.NPP = NPP_monthly;
monthly.ndata = ndata_monthly;
monthly.mdata = mdata_monthly;

%---------calculating Annual mean------------%
fprintf('%s\n', 'calculating Annual climatology...');
tic
NPP_CbPM_SW_annually    = zeros(1080, 2160);
NPP_CAFE_SW_annually    = zeros(1080, 2160);
NPP_CbPM_MODIS_annually = zeros(1080, 2160);
NPP_CAFE_MODIS_annually = zeros(1080, 2160);

for i = 1:1080
    for j = 1:2160
        %
        tmp1 = squeeze(NPP_CbPM_SW_fill_monthly   (i,j,:));
        tmp2 = squeeze(NPP_CAFE_SW_fill_monthly   (i,j,:));
        tmp3 = squeeze(NPP_CbPM_MODIS_fill_monthly(i,j,:));
        tmp4 = squeeze(NPP_CAFE_MODIS_fill_monthly(i,j,:));
        %
        NPP_CbPM_SW_annually   (i,j) = (365/12)*nanmean(tmp1, 'all'); %mg C m-2 d-1 ----> mmol C m-2 yr-1
        NPP_CAFE_SW_annually   (i,j) = (365/12)*nanmean(tmp2, 'all');
        NPP_CbPM_MODIS_annually(i,j) = (365/12)*nanmean(tmp3, 'all');
        NPP_CAFE_MODIS_annually(i,j) = (365/12)*nanmean(tmp4, 'all');
    end
end
toc

%
NPP_annual          = zeros(1080, 2160, 4);
NPP_annual(:,:,1)   = NPP_CbPM_SW_annually;
NPP_annual(:,:,2)   = NPP_CAFE_SW_annually;
NPP_annual(:,:,3)   = NPP_CbPM_MODIS_annually;
NPP_annual(:,:,4)   = NPP_CAFE_MODIS_annually;


%--------bin to OCIM2_48layer grid----------%
% find the latitude band that has no values in all longitude due to sea ice in Arctic
fprintf('%s\n', 'bin to OCIM2_48layer grid...');
tic
vdata_lat_CbPM_SW    = 0;
vdata_lat_CAFE_SW    = 0;
vdata_lat_CbPM_MODIS = 0;
vdata_lat_CAFE_MODIS = 0;
%
for i = 1080:-1:1;
    if ~all(isnan(NPP_CbPM_SW_annually(i,:)));
        vdata_lat_CbPM_SW = i;
        break;
    end
end

for i = 1080:-1:1;
    if ~all(isnan(NPP_CAFE_SW_annually(i,:)));
        vdata_lat_CAFE_SW = i;
        break;
    end
end

for i = 1080:-1:1;
    if ~all(isnan(NPP_CbPM_MODIS_annually(i,:)));
        vdata_lat_CbPM_MODIS = i;
        break;
    end
end

for i = 1080:-1:1;
    if ~all(isnan(NPP_CAFE_MODIS_annually(i,:)));
        vdata_lat_CAFE_MODIS = i;
        break;
    end
end
%
% inpaint_nans
NPP_CbPM_SW_nans    = [NPP_CbPM_SW_annually   , NPP_CbPM_SW_annually(:,1)];
NPP_CAFE_SW_nans    = [NPP_CAFE_SW_annually   , NPP_CAFE_SW_annually(:,1)];
NPP_CbPM_MODIS_nans = [NPP_CbPM_MODIS_annually, NPP_CbPM_MODIS_annually(:,1)];
NPP_CAFE_MODIS_nans = [NPP_CAFE_MODIS_annually, NPP_CAFE_MODIS_annually(:,1)];
%
NPP_CbPM_SW_nans   (vdata_lat_CbPM_SW   :1080, :) = 0;
NPP_CAFE_SW_nans   (vdata_lat_CAFE_SW   :1080, :) = 0;
NPP_CbPM_MODIS_nans(vdata_lat_CbPM_MODIS:1080, :) = 0;
NPP_CAFE_MODIS_nans(vdata_lat_CAFE_MODIS:1080, :) = 0;
%
NPP_CbPM_SW_nans    = inpaint_nans(NPP_CbPM_SW_nans);
NPP_CAFE_SW_nans    = inpaint_nans(NPP_CAFE_SW_nans);
NPP_CbPM_MODIS_nans = inpaint_nans(NPP_CbPM_MODIS_nans);
NPP_CAFE_MODIS_nans = inpaint_nans(NPP_CAFE_MODIS_nans);
%
NPP_CbPM_SW_nans    = NPP_CbPM_SW_nans    (:,1:2160);
NPP_CAFE_SW_nans    = NPP_CAFE_SW_nans    (:,1:2160);
NPP_CbPM_MODIS_nans = NPP_CbPM_MODIS_nans (:,1:2160);
NPP_CAFE_MODIS_nans = NPP_CAFE_MODIS_nans (:,1:2160);
%
[sorted_long, sorted_index] = sort(long, 2);
for i = 1:1080
    NPP_CbPM_SW_nans   (i,:) = NPP_CbPM_SW_nans(i,sorted_index(i,:));   
    NPP_CAFE_SW_nans   (i,:) = NPP_CAFE_SW_nans(i,sorted_index(i,:));
    NPP_CbPM_MODIS_nans(i,:) = NPP_CbPM_MODIS_nans(i,sorted_index(i,:));
    NPP_CAFE_MODIS_nans(i,:) = NPP_CAFE_MODIS_nans(i,sorted_index(i,:));
end
%
NPP_CbPM_SW_48L    = interp2(sorted_long, lat, NPP_CbPM_SW_nans   , long_48layer, lat_48layer);
NPP_CAFE_SW_48L    = interp2(sorted_long, lat, NPP_CAFE_SW_nans   , long_48layer, lat_48layer);
NPP_CbPM_MODIS_48L = interp2(sorted_long, lat, NPP_CbPM_MODIS_nans, long_48layer, lat_48layer);
NPP_CAFE_MODIS_48L = interp2(sorted_long, lat, NPP_CAFE_MODIS_nans, long_48layer, lat_48layer);
%
NPP_CbPM_SW_48L(msk(:,:,1)==0) = 0;
NPP_CAFE_SW_48L(msk(:,:,1)==0) = 0;
NPP_CbPM_MODIS_48L(msk(:,:,1)==0) = 0;
NPP_CAFE_MODIS_48L(msk(:,:,1)==0) = 0;
%
NPP_CbPM_SW_48L(NPP_CbPM_SW_48L < 0) = 0;
NPP_CAFE_SW_48L(NPP_CAFE_SW_48L < 0) = 0;
NPP_CbPM_MODIS_48L(NPP_CbPM_MODIS_48L < 0) = 0;
NPP_CAFE_MODIS_48L(NPP_CAFE_MODIS_48L < 0) = 0;
%
NPP_48layer          = zeros(91, 180, 4);
NPP_48layer(:,:,1)   = NPP_CbPM_SW_48L;
NPP_48layer(:,:,2)   = NPP_CAFE_SW_48L;
NPP_48layer(:,:,3)   = NPP_CbPM_MODIS_48L;
NPP_48layer(:,:,4)   = NPP_CAFE_MODIS_48L;


%-------------Save the file----------------------------%
fileName = 'NPP_climatology.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'monthly', 'NPP_annual', 'NPP_48layer', '-v7.3');

