clc; clear all; close all
addpath('/DFS-L/DATA/primeau/weilewang/my_func/')
addpath('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2')
addpath('/DFS-L/DATA/primeau/weilewang/DATA')
load transport_v4.mat grid M3d
grd90 = grid;
M3d90 = M3d;
%
load OCIM2_CTL_He.mat
load Sobs_91x180x24.mat
load tempobs_91x180x24.mat
load po4obs_91x180x24.mat % WOA PO4 observation
load Siobs_91x180x24.mat  Siobs
load GLODAPv2_talk_91x180x24.mat
load co2syspar90.mat co2syspar

grd91 = output.grid;
M3d = output.M3d;
tmp = M3d;
tmp(:,:,2:end) = 0;
isrf = find(tmp(:));
Pres = M3d+nan;

tmp = M3d90;
tmp(:,:,2:end) = 0;
isrf2 = find(tmp(:));
pres = M3d90+nan;
pres(isrf2) = co2syspar.pres;
FT = interp2(grd90.XT,grd90.YT, pres(:,:,1), grd91.XT, grd91.YT);
Pres(:,:,1) = inpaint_nans(FT);

co2syspar.pres = Pres(isrf);
co2syspar.si = Siobs(isrf);
co2syspar.po4 = po4obs(isrf);
co2syspar.temp = tempobs(isrf);
co2syspar.salt = Sobs(isrf);
co2syspar.alk = talk(isrf);

save co2syspar91 co2syspar
