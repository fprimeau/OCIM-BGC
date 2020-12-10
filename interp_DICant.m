clc; clear all; close all
addpath('/DFS-L/DATA/primeau/weilewang/my_func/')
addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
addpath('/DFS-L/DATA/primeau/weilewang/DATA/')

load DICant.mat
load transport_v4.mat M3d iocn grid
load OCIM2_CTL_He.mat output

grd90 = grid;
grd91 = output.grid;

M3d90 = M3d;
M3d91 = output.M3d;

iwet90 = iocn;
iwet91 = find(M3d91(:));

DIC3d_90 = M3d90+nan;
DIC3d_90(iwet90) = DICant;

DIC3d_91 = interp3(grd90.XT3d, grd90.YT3d, grd90.ZT3d, ...
                   DIC3d_90, ...
                   grd91.XT3d, grd91.YT3d, grd91.ZT3d);




DICant = DIC3d_91;
save DICant_91x180x24 DICant
