%
%
% make initial C & O field for iteration
% using initial field of 24 layer and bin to OCIM2_48layer grid
clc; clear all; close all;
addpath ('../../DATA/BGC_2023Nature/');
addpath ('../../DATA/BGC_48layer/');
addpath ('../src');
load initCO_91x180x24.mat

DIC_24 = data.DIC;
DOC_24 = data.DOC;
PIC_24 = data.PIC;
POC_24 = data.POC;
ALK_24 = data.ALK;
O2_24  = data.O2 ;

load OCIM2_CTL_He_48layer.mat output
msk = output.M3d;

%inpaint_nans for each horizontal field
for i = 1:24;
    DIC_24(:,:,i) = inpaint_nans(DIC_24(:,:,i));
    DOC_24(:,:,i) = inpaint_nans(DOC_24(:,:,i));
    PIC_24(:,:,i) = inpaint_nans(PIC_24(:,:,i));
    POC_24(:,:,i) = inpaint_nans(POC_24(:,:,i));
    ALK_24(:,:,i) = inpaint_nans(ALK_24(:,:,i));
    O2_24 (:,:,i) = inpaint_nans(O2_24 (:,:,i));
end

DIC_int = zeros(91, 180, 48);
DOC_int = zeros(91, 180, 48);
PIC_int = zeros(91, 180, 48);
POC_int = zeros(91, 180, 48);
ALK_int = zeros(91, 180, 48);
O2_int  = zeros(91, 180, 48);

for i = 1:91
    for j = 1:180
        DIC_vertical = squeeze(DIC_24(i, j, :));
        DOC_vertical = squeeze(DOC_24(i, j, :));
        PIC_vertical = squeeze(PIC_24(i, j, :));
        POC_vertical = squeeze(POC_24(i, j, :));
        ALK_vertical = squeeze(ALK_24(i, j, :));
        O2_vertical  = squeeze(O2_24 (i, j, :));
        
        DIC_interp = interp1(1:24, DIC_vertical, linspace(1, 24, 48), 'linear');
        DOC_interp = interp1(1:24, DOC_vertical, linspace(1, 24, 48), 'linear');
        PIC_interp = interp1(1:24, PIC_vertical, linspace(1, 24, 48), 'linear');
        POC_interp = interp1(1:24, POC_vertical, linspace(1, 24, 48), 'linear');
        ALK_interp = interp1(1:24, ALK_vertical, linspace(1, 24, 48), 'linear');
        O2_interp  = interp1(1:24, O2_vertical , linspace(1, 24, 48), 'linear');
                
        DIC_int(i, j, :) = DIC_interp;
        DOC_int(i, j, :) = DOC_interp;
        PIC_int(i, j, :) = PIC_interp;
        POC_int(i, j, :) = POC_interp;
        ALK_int(i, j, :) = ALK_interp;
        O2_int (i, j, :) = O2_interp ;
    end
end

DIC_int(msk(:)==0) = NaN;
DOC_int(msk(:)==0) = NaN;
PIC_int(msk(:)==0) = NaN;
POC_int(msk(:)==0) = NaN;
ALK_int(msk(:)==0) = NaN;
O2_int (msk(:)==0) = NaN;

data.DIC = DIC_int;
data.DOC = DOC_int;
data.PIC = PIC_int;
data.POC = POC_int;
data.ALK = ALK_int;
data.O2  = O2_int ;

data = rmfield(data, "DIP");
data = rmfield(data, "POP");
data = rmfield(data, "DOP");

fileName = 'Initial_COfield_91x180x48.mat'
directory = '../../DATA/BGC_48layer'
filePath = fullfile(directory, fileName);
save(filePath, 'data');
