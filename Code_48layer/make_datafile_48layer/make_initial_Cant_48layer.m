% estimate initial Cant field in 48layer from Weilei's Cant filed in 24layer
% Interp DIC result of 48layer will be used to re-optimzied for the anthropogenic carbon field

clc; clear all; close all
addpath('../src/')
addpath('../../DATA/BGC_2023Nature/');
addpath('../../DATA/BGC_48layer/');
load GLODAPv2_DIC_remove_cant.mat       %observed dic field in 24layer from weilei's paper
                                        %DIC field with anthropogenic DIC removed according to 14C method
load OCIM2_CTL_He.mat
output_24layer = output;
clear output
msk_24layer    = output_24layer.M3d;
XT3d_24layer   = output_24layer.grid.XT3d;
YT3d_24layer   = output_24layer.grid.YT3d;
ZT3d_24layer   = output_24layer.grid.ZT3d;
                                        
                                      
load OCIM2_CTL_He_48layer.mat
output_48layer = output;
clear output
msk_48layer      = output_48layer.M3d;
XT3d_48layer     = output_48layer.grid.XT3d;
YT3d_48layer     = output_48layer.grid.YT3d;
ZT3d_48layer     = output_48layer.grid.ZT3d;
dzt_48layer      = output_48layer.grid.dzt;

dicraw_24layer = nanmean(gc12new, 4);

%inpaint_nans for each horizontal field
for i = 1:size(msk_24layer, 3);
    dicraw_24layer(:,:,i) = inpaint_nans(dicraw_24layer(:,:,i));
end

%-------> Nan value 없어짐. land까지 전부 interpolate 됨.
% 91x180 각 gird 별로 interp함수 이용해서 interpolation.

dic_initial = zeros(91, 180, 48);

for i = 1:91
    for j = 1:180
        vertical = squeeze(dicraw_24layer(i, j, :));
        vertical_interp = interp1(1:24, vertical, linspace(1, 24, 48), 'linear');
        dic_initial(i, j, :) = vertical_interp;
    end
end

dic_initial(msk_48layer(:)==0) = NaN;

%-------------Save the file----------------------------%
fileName = 'dic_initial_91x180x48.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'dic_initial');
