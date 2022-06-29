clc; clear all; close all
%global iter
%iter = 0 ;
on   = true  ;
off  = false ;
format long

%ver = datestr(now,'mmmdd');
RunVer = 'testPobs_CTL_He_PCCella1b1_DOC0.25_DOP0';

GridVer  = 91  ;
operator = 'A' ;

%model output directory
outputDir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/', GridVer);
figDir = strcat(outputDir,'FIGS_testPobs_PCCell/a1b1_');
%outPath = figDir;

% load model output fields
par.fname = strcat(outputDir, RunVer, '.mat');
load(par.fname);
model = data;

% load optimal parameter values
par.fxhat = strcat(outputDir, RunVer,'_xhat.mat');
load(par.fxhat);

%perturbed parameter fields
%pfname    = strcat(outputDir, RunVer,'_pert.mat');

par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

crosssect_lon = 110; % Atlantic = 170; Pacific = 110

%-------------load data and set up parameters---------------------
par.LoadOpt = on;
SetUp ;
xhat
iwet = par.iwet;

% parameter indices
pindx = xhat.pindx;

% will need to run PackPar. or find other way to create p0
% par = SetPar(par);

%-------------set up figure properties------------------
set(groot,'defaultAxesFontName','Times',...
    'defaultAxesFontSize',14,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','on',...
    'defaultAxesYMinorTick','on');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Times',...
    'defaultTextInterpreter','latex');

%Colors
aqua = [0.2 0.8 0.8];
teal = [0 128/255 128/255];
darkgreen = [0 100/255 0];

lblue = [0 191/255 255/255];
navy = [ 0 0 128/255];
%%% ------------ OLD CODE ----------------------------
%{
GridVer  = 91  ;
operator = 'A' ;
% GridVer: choose from 90 and 91; Ver 90 is for a Transport
% operator without diapycnal mixing but optimized using DIP ;
% Ver 91 include a bunch of operators that include diapycnal
% mixing. These operators represent sensiviity tests on He
% constraint and on mixing parameterizations (DeVries et al, 2018).
% A -> CTL_He; B -> CTL_noHe; C -> KiHIGH_He; D -> KiHIGH_noHe;
% E -> KvHIGH_KiLOW_He; F -> KvHIGH_KiLOW_noHe; G -> KiLOW_He;
% H -> KiLOW_noHe; I -> KvHIGH_He; J -> KvHIGH_noHe; K -> KvHIGH_KiHIGH_noHe
Gtest = off ;
Htest = off ;
par.optim   = on ;
par.Cmodel  = on ;
par.Omodel  = on ;
par.Simodel = off ;
par.LoadOpt = off ; % if load optimial par.
par.pscale  = 0.0 ;
par.cscale  = 0.0 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigma = on ;
par.opt_kP_T  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ;
par.opt_bP    = on ;
par.opt_beta  = on ;
par.opt_alpha = on ;
% C model parameters
par.opt_bC_T  = on ;
par.opt_bC    = on ;
par.opt_d     = on ;
par.opt_kC_T  = on ;
par.opt_kdC   = on ;
par.opt_R_Si  = on ;
par.opt_rR    = on ;
par.opt_cc    = on ;
par.opt_dd    = on ;
% O model parameters
par.opt_O2C_T = on ;
par.opt_rO2C  = on ;
par.opt_O2P_T = on ;
par.opt_rO2P  = on ;
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;

% --------------load data and set up parameters -----------
SetUp ;

% save results
% ATTENTION: Change this direcrtory to where you wanna
% save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer);
elseif isunix
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
    % 'MSK%2d/'],GridVer);
    output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
                        'MSK%2d/PME4DICALK/'],GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
    % 'MSK%2d/'],GridVer);
end
VER = strcat(output_dir,TRdivVer);
% Creat output file names based on which model(s) is(are) optimized
if Gtest == on
    fname = strcat(VER,'_GHtest');
elseif Gtest == off
    if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
        fname = strcat(VER,'_P');
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
        base_name = strcat(VER,'_PCv2');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
        base_name = strcat(VER,'_PCOv1');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
        base_name = strcat(VER,'_PCSi');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
        base_name = strcat(VER,'_PCOSi');
        catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
        fname = strcat(base_name,catDOC);
    end
end
pfname    = strcat(fname,'_pert.mat');
par.fname = fname ;
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ;
load(fname)
load(fxhat)

%---------------- inital guesses on C and O ---------------
DIC = data.DIC - par.dicant ;

GC  = [DIC(iwet); data.POC(iwet); data.DOC(iwet); ...
       data.PIC(iwet); data.ALK(iwet)];
% GC  = GC + 1e-6*randn(5*nwet,1) ;
if par.Omodel == on
    GO  = real(data.O2(iwet)) + 1e-9*randn(par.nwet,1);
end
%--------------------- prepare parameters ------------------
if par.optim == on
    % load optimal parameters from a file or set them to default values
    par = SetPar(par) ;
    % pack parameters into an array, assign them corresponding indices.
    par = PackPar(par) ;
end
p0 = par.p0 ;

% ---------------------perturb temperature ---------------
vT0  = par.Temp(iwet) ;
vT1  = vT0 + 2        ;
Tz1  = (vT1-min(vT0))/(max(vT0)-min(vT0));
Tz0  = (vT0-min(vT0))/(max(vT0)-min(vT0));
Tz3d = M3d + nan      ;
Tz3d(iwet) = Tz1      ;
par.aveT   = nanmean(Tz3d(:,:,1:3),3) ;
par.Tz     = Tz1*1e-8 ;
%}

% if isfile(pfname)
%     load(pfname)
%     pDIP = pdata.DIP ;
%     pDOP = pdata.DOP ;
%     pPOP = pdata.POP ;
% else
% 	%iter = 5;
% 	%[f, fx, fxx, pdata] = neglogpost(x, par)
%     %[par, P ] = eqPcycle(p0, par) ;
% 	[par, P, Px, Pxx] = eqPcycle(x0, par) ;
%     %
%     pDIP = M3d+nan ;  pDIP(iwet) = P(1+0*nwet:1*nwet) ;
%     pPOP = M3d+nan ;  pPOP(iwet) = P(1+1*nwet:2*nwet) ;
%     pDOP = M3d+nan ;  pDOP(iwet) = P(1+2*nwet:3*nwet) ;
% 	pdata.DIP = pDIP; pdata.DOP = pDOP; pdata.POP = pPOP;
% 	pdata.Px   = Px  ;
%     pdata.Pxx  = Pxx ;
% 	fprintf('saving perturbed model solution to file: %s \n',pfname)
%     save(pfname, 'pdata')
% end
% par.DIP  = pDIP(iwet) ;

%--- Load Phosphorus Model Output ----
DIP = model.DIP ;
POP = model.POP ;
DOP = model.DOP ;

%--- Split Px into DIP, DOP, and POP
DIPx = model.Px(1+0*nwet:1*nwet,:);
POPx = model.Px(1+1*nwet:2*nwet,:) ;
DOPx = model.Px(1+2*nwet:3*nwet,:) ;

%--------define S patterns for DIP ----------

%derivative of DIP field wrt log(sigma)
if isfield(xhat,'sigma')
	S_sigma = M3d+nan;
	S_sigma(iwet) = DIPx(:,pindx.lsigma);
end

if isfield(xhat,'kP_T')
	S_kP_T = M3d+nan;
	S_kP_T(iwet) = DIPx(:,pindx.kP_T);
end

if isfield(xhat,'kdP')
	S_kdP = M3d+nan;
	S_kdP(iwet) = DIPx(:,pindx.lkdP);
end

if isfield(xhat,'bP_T')
	S_bP_T = M3d+nan;
	S_bP_T(iwet) = DIPx(:,pindx.bP_T);
end

if isfield(xhat,'bP')
	S_bP = M3d+nan;
	S_bP(iwet) = DIPx(:,pindx.lbP);
end

if isfield(xhat,'alpha')
	S_alpha = M3d+nan;
	S_alpha(iwet) = DIPx(:,pindx.lalpha);
end

if isfield(xhat,'beta')
	S_beta = M3d+nan;
	S_beta(iwet) = DIPx(:,pindx.lbeta);
end

nfig = 0 ;
% ---- make zonal cross section layout for DIP sensitivity to P params ---------
%crosssect_lon = 110;	%set longitude in pacific
aspectR = [1 0.6 1];  %define aspect ratio for plots
nfig = nfig + 1 ;
figure(nfig);
t_crosssect = tiledlayout('flow','TileSpacing','compact','Padding','compact'); 	% create tiled layout

if isfield(xhat,'sigma')
	nexttile;
	S_sigma_min = -0.1;
	S_sigma_max = 0.2;

	contourf(grd.yt,-grd.zt,squeeze(S_sigma(:,crosssect_lon,:))',[S_sigma_min:0.01:S_sigma_max])
	set(gca,'color','black')
	colorbar;
	colormap(gca,darkb2r(S_sigma_min,S_sigma_max));
	xlabel('latitude (deg)');
	ylabel('depth (m)')
	pbaspect(aspectR);
	%t = sprintf('DIP sensitivity to sigma (dDIP/dlsigma): lon = %4.1f deg', crosssect_lon);
	title('dDIP/dlsigma');
end

if isfield(xhat,'kP_T')
	nexttile;
	%S_kP_T_min = min(S_kP_T,[],'all');
	%S_kP_T_max = max(S_kP_T,[],'all');
	[S_kP_T_min, S_kP_T_max] = bounds(S_kP_T,'all');
	nlevs = 25;
	levels = linspace(S_kP_T_min,S_kP_T_max,nlevs);

	contourf(grd.yt,-grd.zt,squeeze(S_kP_T(:,crosssect_lon,:))',nlevs)
	set(gca,'color','black')
	colorbar
	colormap(gca,darkb2r(S_kP_T_min, S_kP_T_max));
	xlabel('latitude (deg)');
	ylabel('depth (m)')
	pbaspect(aspectR);
	%t = sprintf('DIP sensitivity to kP_T (dDIP/dkP_T): lon = %4.1f deg', crosssect_lon);
	title('dDIP/dkPT');
end

if isfield(xhat,'kdP')
	nexttile;
	[S_kdP_min, S_kdP_max] = bounds(S_kdP,'all');
	nlevs = 25;
	contourf(grd.yt,-grd.zt,squeeze(S_kdP(:,crosssect_lon,:))',nlevs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_kdP_min, S_kdP_max));
	colorbar;
	xlabel('latitude (deg)');
	ylabel('depth (m)');
	pbaspect(aspectR);
	title('dDIP/dlkdP');
end

if isfield(xhat,'bP_T')
	nexttile;
	[S_bP_T_min, S_bP_T_max] = bounds(S_bP_T,'all');
	nlevs = 25;
	contourf(grd.yt,-grd.zt,squeeze(S_bP_T(:,crosssect_lon,:))',nlevs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_bP_T_min, S_bP_T_max));
	colorbar;
	xlabel('latitude (deg)');
	ylabel('depth (m)');
	pbaspect(aspectR);
	title('dDIP/dbPT');
end

if isfield(xhat,'bP')
	nexttile;
	[S_bP_min, S_bP_max] = bounds(S_bP,'all');
	nlevs = 25;
	contourf(grd.yt,-grd.zt,squeeze(S_bP(:,crosssect_lon,:))',nlevs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_bP_min, S_bP_max));
	colorbar;
	xlabel('latitude (deg)');
	ylabel('depth (m)');
	pbaspect(aspectR);
	title('dDIP/dlbP');
end

if isfield(xhat,'alpha')
	nexttile;
	[S_alpha_min, S_alpha_max] = bounds(S_alpha,'all');
	nlevs = 25;
	contourf(grd.yt,-grd.zt,squeeze(S_alpha(:,crosssect_lon,:))',nlevs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_alpha_min, S_alpha_max));
	colorbar;
	xlabel('latitude (deg)');
	ylabel('depth (m)');
	pbaspect(aspectR);
	title('dDIP/dlalpha');
end

if isfield(xhat,'beta')
	nexttile;
	[S_beta_min, S_beta_max] = bounds(S_beta,'all');
	nlevs = 25;
	contourf(grd.yt,-grd.zt,squeeze(S_beta(:,crosssect_lon,:))',nlevs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_beta_min, S_beta_max));
	colorbar;
	xlabel('latitude (deg)');
	ylabel('depth (m)');
	pbaspect(aspectR);
	title('dDIP/dlbeta');
end

t = sprintf('DIP sensitivity Zonal Cross Section: lon = %4.1f deg', grd.xt(crosssect_lon));
title(t_crosssect,t);	% Add layout title
set(gcf, 'InvertHardcopy', 'off')
figTitle = ['DIPcrosssect_all' num2str(grd.xt(crosssect_lon))];
print(gcf,[figDir 'SensiTest_' figTitle '.png'],'-dpng')


% % automation ideas
% % how to give only the left plots depth labels and only the bottom plots in layout latitude labels
% nax = 0
% nax = nax+1;
% allaxes{nax} =gca;
% set(gca, 'XTickLabel', []);
%
% t = sprintf('Surface DIP sensitivity to parameters');
% title(t_surface,t);	% Add layout title
% [nrows ncols] = t_surface.GridSize;
%
% leftplot_ind = [];
% for i = 1:nrows
% 	leftplot_ind = [leftplot_ind; i +(i-1)*ncols];
% end

% ---- make surface map layout for DIP sensitivity to P params ---------
nfig = nfig + 1 ;
figure(nfig);
aspectR = [1 0.58 1];  %define aspect ratio for plots
t_surface = tiledlayout('flow','TileSpacing','compact','Padding','compact'); 	% create tiled layout

if isfield(xhat,'sigma')
    nexttile;
    [S_sigma_min, S_sigma_max] = bounds(S_sigma(:,:,1:2),'all');
	nlevs = 11;
    levs = [S_sigma_min:0.1:S_sigma_max];
	contourf(grd.xt,grd.yt,mean(S_sigma(:,:,1:2),3,'omitnan'),levs)
    set(gca,'color','black')
    colormap(darkb2r(S_sigma_min, S_sigma_max)), colorbar
    pbaspect(aspectR);
    xlabel('longitude (deg E)');
	ylabel('latitude (deg)');
	title('dDIP/dlsigma');
    %figTitle = 'DIP_sigma_surf';
    %print(gcf,[figDir 'SensiTest_' figTitle '.png'],'-dpng')
end

if isfield(xhat,'kP_T')
	nexttile;
    [S_kP_T_min, S_kP_T_max] = bounds(S_kP_T(:,:,1:2),'all');
	nlevs = 11;
    levs = [S_kP_T_min:0.1:S_kP_T_max];
	contourf(grd.xt,grd.yt,mean(S_kP_T(:,:,1:2),3,'omitnan'),levs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_kP_T_min, S_kP_T_max));
	colorbar;
    pbaspect(aspectR);
    xlabel('longitude (deg E)');
	ylabel('latitude (deg)');
	title('dDIP/dkPT');
end

if isfield(xhat,'kdP')
	nexttile;
    [S_kdP_min, S_kdP_max] = bounds(S_kdP(:,:,1:2),'all');
	nlevs = 11;
    levs = [S_kdP_min:0.1:S_kdP_max];
	contourf(grd.xt,grd.yt,mean(S_kdP(:,:,1:2),3,'omitnan'),levs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_kdP_min, S_kdP_max));
	colorbar;
    pbaspect(aspectR);
    xlabel('longitude (deg E)');
	ylabel('latitude (deg)');
	title('dDIP/dlkdP');
end

if isfield(xhat,'bP_T')
	nexttile;
	[S_bP_T_min, S_bP_T_max] = bounds(S_bP_T(:,:,1:2),'all');
	nlevs = 11;
    levs = [S_bP_T_min:0.1:S_bP_T_max];
	contourf(grd.xt,grd.yt,mean(S_bP_T(:,:,1:2),3,'omitnan'),levs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_bP_T_min, S_bP_T_max));
	colorbar;
    pbaspect(aspectR);
	xlabel('longitude (deg E)');
	ylabel('latitude (deg)');
	title('dDIP/dbPT');
end

if isfield(xhat,'bP')
	nexttile;
	[S_bP_min, S_bP_max] = bounds(S_bP(:,:,1:2),'all');
	nlevs = 11;
    levs = [S_bP_min:0.1:S_bP_max];
	contourf(grd.xt,grd.yt,mean(S_bP(:,:,1:2),3,'omitnan'),levs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_bP_min, S_bP_max));
	colorbar;
    pbaspect(aspectR);
	xlabel('longitude (deg E)');
	ylabel('latitude (deg)');
	title('dDIP/dlbP');
end

if isfield(xhat,'alpha')
	nexttile;
    [S_alpha_min, S_alpha_max] = bounds(S_alpha(:,:,1:2),'all');
	nlevs = 11;
    levs = [S_alpha_min:0.1:S_alpha_max];
	contourf(grd.xt,grd.yt,mean(S_alpha(:,:,1:2),3,'omitnan'),levs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_alpha_min, S_alpha_max));
	colorbar;
    pbaspect(aspectR);
    xlabel('longitude (deg E)');
	ylabel('latitude (deg)');
	title('dDIP/dlalpha');
end

if isfield(xhat,'beta')
	nexttile;
    [S_beta_min, S_beta_max] = bounds(S_beta(:,:,1:2),'all');
	nlevs = 11;
    levs = [S_beta_min:0.1:S_beta_max];
	contourf(grd.xt,grd.yt,mean(S_beta(:,:,1:2),3,'omitnan'),levs)
	set(gca,'color','black');
	colormap(gca,darkb2r(S_beta_min, S_beta_max));
	colorbar;
    pbaspect(aspectR);
    xlabel('longitude (deg E)');
	ylabel('latitude (deg)');
	title('dDIP/dlbeta');
end

set(gcf, 'InvertHardcopy', 'off')
figTitle = 'DIPsurface_all';
print(gcf,[figDir 'SensiTest_' figTitle '.png'],'-dpng')

%  ---------------------------------------------------
%% ---------- C2P sensitivity ------------------------
%
if par.Cellmodel ==on
%--- Load C2P Model Output ----
	C2Px = model.C2Px;

%--------define S patterns for C2P ----------
	clear S_*
	% --- Pmodel parameters -------------
	if isfield(xhat,'sigma')
		%derivative of C2P field wrt log(sigma)
		S_sigma = M3d+nan;
		S_sigma(iwet) = C2Px(:,pindx.lsigma);
	end
	if isfield(xhat,'kP_T')
		S_kP_T = M3d+nan;
		S_kP_T(iwet) = C2Px(:,pindx.kP_T);
	end
	if isfield(xhat,'kdP')
		S_kdP = M3d+nan;
		S_kdP(iwet) = C2Px(:,pindx.lkdP);
	end
	if isfield(xhat,'bP_T')
		S_bP_T = M3d+nan;
		S_bP_T(iwet) = C2Px(:,pindx.bP_T);
	end
	if isfield(xhat,'bP')
		S_bP = M3d+nan;
		S_bP(iwet) = C2Px(:,pindx.lbP);
	end
	if isfield(xhat,'alpha')
		S_alpha = M3d+nan;
		S_alpha(iwet) = C2Px(:,pindx.lalpha);
	end
	if isfield(xhat,'beta')
		S_beta = M3d+nan;
		S_beta(iwet) = C2Px(:,pindx.lbeta);
	end
	% ------- cell model parameters ----------------
	if isfield(xhat,'Q10Photo')
	   S_Q10 = M3d+nan;
	   S_Q10(iwet) = C2Px(:,pindx.lQ10Photo);
   end
   if isfield(xhat,'fStorage')
	   S_fStor = M3d+nan;
	   S_fStor(iwet) = C2Px(:,pindx.lfStorage);
   end
   if isfield(xhat,'fRibE')
	   S_fRibE = M3d+nan;
	   S_fRibE(iwet) = C2Px(:,pindx.tfRibE);
   end
   if isfield(xhat,'kST0')
	   S_kST0 = M3d+nan;
	   S_kST0(iwet) = C2Px(:,pindx.lkST0);
   end
   if isfield(xhat,'PLip_PCutoff')
	   S_PCutoff = M3d+nan;
	   S_PCutoff(iwet) = C2Px(:,pindx.lPLip_PCutoff);
   end
   if isfield(xhat,'PLip_scale')
	   S_PLipscale = M3d+nan;
	   S_PLipscale(iwet) = C2Px(:,pindx.lPLip_scale);
   end
   if isfield(xhat,'PStor_rCutoff')
	   S_rCutoff = M3d+nan;
	   S_rCutoff(iwet) = C2Px(:,pindx.lPStor_rCutoff);
   end
   if isfield(xhat,'PStor_scale')
	   S_PStorscale = M3d+nan;
	   S_PStorscale(iwet) = C2Px(:,pindx.lPStor_scale);
   end
   if isfield(xhat,'alphaS')
	   S_alphaS = M3d+nan;
	   S_alphaS(iwet) = C2Px(:,pindx.lalphaS);
   end
   if isfield(xhat,'gammaDNA')
	   S_gammaDNA = M3d+nan;
	   S_gammaDNA(iwet) = C2Px(:,pindx.tgammaDNA);
   end

   %% --- plot surface C2P sensitivity -----------
   nfig = nfig + 1 ;
   figure(nfig);
   aspectR = [1 0.58 1];  %define aspect ratio for plots
   t_surfaceC2P = tiledlayout('flow','TileSpacing','compact','Padding','compact'); 	% create tiled layout
   set(gcf,'Units','normalized','Position',[0 0.2 0.6 0.6])

   cstep = 25;
   if isfield(xhat,'sigma')
	   nexttile;
	   [S_sigma_min, S_sigma_max] = bounds(S_sigma(:,:,1:2),'all');
	   nlevs = 11;
	   levs = [S_sigma_min:cstep:S_sigma_max];
	   contourf(grd.xt,grd.yt,mean(S_sigma(:,:,1:2),3,'omitnan'),nlevs)
	   set(gca,'color','black')
	   colormap(darkb2r(S_sigma_min, S_sigma_max)), colorbar
	   xlabel('longitude (deg E)');
	   ylabel('latitude (deg)');
	   pbaspect(aspectR);
	   title('dC2P/dlsigma');
   end
   if isfield(xhat,'kP_T')
	   nexttile;
	   [S_kP_T_min, S_kP_T_max] = bounds(S_kP_T(:,:,1:2),'all');
	   nlevs = 11;
	   levs = [S_kP_T_min:cstep:S_kP_T_max];
	   contourf(grd.xt,grd.yt,mean(S_kP_T(:,:,1:2),3,'omitnan'),nlevs)
	   set(gca,'color','black');
	   colormap(gca,darkb2r(S_kP_T_min, S_kP_T_max));
	   colorbar;
	   xlabel('longitude (deg E)');
	   ylabel('latitude (deg)');
	   pbaspect(aspectR);
	   title('dC2P/dkPT');
   end
   if isfield(xhat,'kdP')
	   nexttile;
	   %[S_kdP_min, S_kdP_max] = bounds(S_kdP(:,:,1:2),'all');
	   S_kdP_min = prctile(S_kdP(:,:,1:2),1,'all');
	   S_kdP_max = prctile(S_kdP(:,:,1:2),99,'all');
	   nlevs = 11;
	   levs = [S_kdP_min:cstep:S_kdP_max];
	   contourf(grd.xt,grd.yt,mean(S_kdP(:,:,1:2),3,'omitnan'),levs)
	   set(gca,'color','black');
	   colormap(gca,darkb2r(S_kdP_min, S_kdP_max));
	   colorbar;
	   xlabel('longitude (deg E)');
	   ylabel('latitude (deg)');
	   pbaspect(aspectR);
	   title('dC2P/dlkdP');
   end
   if isfield(xhat,'bP_T')
	   nexttile;
	   %[S_bP_T_min, S_bP_T_max] = bounds(S_bP_T(:,:,1:2),'all');
	   S_bP_T_min = prctile(S_bP_T(:,:,1:2),1,'all');
	   %S_bP_T_max = prctile(S_bP_T(:,:,1:2),1,'all');
	   S_bP_T_max = 1;

	   nlevs = 11;
	   levs = [S_bP_T_min:cstep:S_bP_T_max];
	   contourf(grd.xt,grd.yt,mean(S_bP_T(:,:,1:2),3,'omitnan'),levs)
	   set(gca,'color','black');
	   colormap(gca,darkb2r(S_bP_T_min, S_bP_T_max));
	   colorbar;
	   xlabel('longitude (deg E)');
	   ylabel('latitude (deg)');
	   pbaspect(aspectR);
	   title('dC2P/dbPT');
   end
   if isfield(xhat,'bP')
	   nexttile;
	   %[S_bP_min, S_bP_max] = bounds(S_bP(:,:,1:2),'all');
	   S_bP_min = prctile(S_bP(:,:,1:2),1,'all');
	   %S_bP_max = prctile(S_bP(:,:,1:2),99,'all');
	   S_bP_max = 1;
	   nlevs = 11;
	   levs = [S_bP_min:cstep:S_bP_max];
	   contourf(grd.xt,grd.yt,mean(S_bP(:,:,1:2),3,'omitnan'),levs)
	   set(gca,'color','black');
	   colormap(gca,darkb2r(S_bP_min, S_bP_max));
	   colorbar;
	   xlabel('longitude (deg E)');
	   ylabel('latitude (deg)');
	   pbaspect(aspectR);
	   title('dC2P/dlbP');
   end
   if isfield(xhat,'alpha')
	   nexttile;
	   %[S_alpha_min, S_alpha_max] = bounds(S_alpha(:,:,1:2),'all');
	   S_alpha_max = prctile(S_alpha(:,:,1:2),99,'all');
	   S_alpha_min = min(prctile(S_alpha(:,:,1:2),1,'all'),0);
	   nlevs = 11;
	   levs = [S_alpha_min:cstep:S_alpha_max];
	   %levs = prctile(S_alpha(:,:,1:2),[5:10:95],'all');
	   %levs= 0:25:200;
	   contourf(grd.xt,grd.yt,mean(S_alpha(:,:,1:2),3,'omitnan'),levs)
	   set(gca,'color','black');
	   colormap(gca,darkb2r(S_alpha_min, S_alpha_max));
	   colorbar;
	   xlabel('longitude (deg E)');
	   ylabel('latitude (deg)');
	   pbaspect(aspectR);
	   title('dC2P/dlalpha');
   end
   if isfield(xhat,'beta')
	   nexttile;
	   %[S_beta_min, S_beta_max] = bounds(S_beta(:,:,1:2),'all');
	   S_beta_min = prctile(S_beta(:,:,1:2),1,'all');
	   S_beta_max = 1;
	   nlevs = 11;
	   %levs = [S_beta_min:cstep:S_beta_max];
	   levs = -150:25:0;
	   contourf(grd.xt,grd.yt,mean(S_beta(:,:,1:2),3,'omitnan'),levs)
	   set(gca,'color','black');
	   colormap(gca,darkb2r(S_beta_min, S_beta_max));
	   colorbar;
	   xlabel('longitude (deg E)');
	   ylabel('latitude (deg)');
	   pbaspect(aspectR);
	   title('dC2P/dlbeta');
   end
   % Cell model parameters
    if isfield(xhat,'Q10Photo')
        nexttile;
        [S_Q10_min, S_Q10_max] = bounds(S_Q10(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_Q10_min:cstep:S_Q10_max];
	    contourf(grd.xt,grd.yt,mean(S_Q10(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_Q10_min, S_Q10_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dlQ10Photo');
    end
    if isfield(xhat,'fStorage')
        nexttile;
        [S_fStor_min, S_fStor_max] = bounds(S_fStor(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_fStor_min:cstep:S_fStor_max];
	    contourf(grd.xt,grd.yt,mean(S_fStor(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_fStor_min, S_fStor_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dlfStorage');
    end
    if isfield(xhat,'fRibE')
        nexttile;
        [S_fRibE_min, S_fRibE_max] = bounds(S_fRibE(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_fRibE_min:cstep:S_fRibE_max];
	    contourf(grd.xt,grd.yt,mean(S_fRibE(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_fRibE_min, S_fRibE_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dtanh(fRibE)');
    end
    if isfield(xhat,'kST0')
        nexttile;
        [S_kST0_min, S_kST0_max] = bounds(S_kST0(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_kST0_min:cstep:S_kST0_max];
	    contourf(grd.xt,grd.yt,mean(S_kST0(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_kST0_min, S_kST0_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dlkST0');
    end
    if isfield(xhat,'PLip_PCutoff')
        nexttile;
        [S_PCutoff_min, S_PCutoff_max] = bounds(S_PCutoff(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_PCutoff_min:0.5:S_PCutoff_max];
	    contourf(grd.xt,grd.yt,mean(S_PCutoff(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_PCutoff_min, S_PCutoff_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dlPLip-PCutoff');
    end
    if isfield(xhat,'PLip_scale')
        nexttile;
        [S_PLipscale_min, S_PLipscale_max] = bounds(S_PLipscale(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_PLipscale_min:1:S_PLipscale_max];
	    contourf(grd.xt,grd.yt,mean(S_PLipscale(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_PLipscale_min, S_PLipscale_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dlPLipscale');
    end
    if isfield(xhat,'PStor_rCutoff')
        nexttile;
        [S_rCutoff_min, S_rCutoff_max] = bounds(S_rCutoff(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_rCutoff_min:1:S_rCutoff_max];
	    contourf(grd.xt,grd.yt,mean(S_rCutoff(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_rCutoff_min, S_rCutoff_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dlPStor-rCutoff');
    end
    if isfield(xhat,'PStor_scale')
        nexttile;
        [S_PStorscale_min, S_PStorscale_max] = bounds(S_PStorscale(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_PStorscale_min:cstep:S_PStorscale_max];
	    contourf(grd.xt,grd.yt,mean(S_PStorscale(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_PStorscale_min, S_PStorscale_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dlPStorscale');
    end
    if isfield(xhat,'alphaS')
        nexttile;
        [S_alphaS_min, S_alphaS_max] = bounds(S_alphaS(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_alphaS_min:cstep:S_alphaS_max];
	    contourf(grd.xt,grd.yt,mean(S_alphaS(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_alphaS_min, S_alphaS_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dlalphaS');
    end
	if isfield(xhat,'gammaDNA')
        nexttile;
        [S_gammaDNA_min, S_gammaDNA_max] = bounds(S_gammaDNA(:,:,1:2),'all');
	    nlevs = 11;
        levs = [S_gammaDNA_min:cstep:S_gammaDNA_max];
	    contourf(grd.xt,grd.yt,mean(S_gammaDNA(:,:,1:2),3,'omitnan'),nlevs)
	    set(gca,'color','black');
	    colormap(gca,darkb2r(S_gammaDNA_min, S_gammaDNA_max));
	    colorbar;
        xlabel('longitude (deg E)');
	    ylabel('latitude (deg)');
        pbaspect(aspectR);
	    title('dC2P/dtanh(gammaDNA)');
    end

    title(t_surfaceC2P,'Sensitivity of Surface C2P to model parameters');	% Add layout title
    set(gcf, 'InvertHardcopy', 'off')
    figTitle = 'C2Psurface_all';
    print(gcf,[figDir 'SensiTest_' figTitle '.png'],'-dpng')
	exportfig(gcf,[figDir 'SensiTest_fig_' figTitle],'fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')

end

% --------- make sigma figures ------------------
%{
S_sigma_min = -0.1;
S_sigma_max = 0.2;

% make a zonal cross section
crosssect_lon = 170;
nfig = nfig + 1 ;
figure(nfig);
contourf(grd.yt,-grd.zt,squeeze(S_sigma(:,crosssect_lon,:))',[S_sigma_min:0.01:S_sigma_max])
set(gca,'color','black')
colormap(darkb2r(S_sigma_min,S_sigma_max));
colorbar
xlabel('latitude (deg)');
ylabel('depth (m)')
t = sprintf('DIP sensitivity to sigma (dDIP/dsigma): lon = %4.1f deg', crosssect_lon);
title(t);
set(gcf, 'InvertHardcopy', 'off')
figTitle = ['DIPsigma' num2str(crosssect_lon)];
print(gcf,[figDir 'SensiTest_' figTitle '.png'],'-dpng')

% make a zonal average of age for the Pacific basin
nfig = nfig + 1 ;
figure(nfig);
PAC = MSKS.PAC;
PZA = squeeze(nansum(PAC.*S_sigma.*dVt,2)./sum(PAC.*dVt,2))';
% surface
subplot(2,1,1) ;
contourf(grd.yt,-grd.zt(1:9),PZA(1:9,:),[S_sigma_min:0.01:S_sigma_max]);
set(gca,'color','black')
colormap(darkb2r(S_sigma_min,S_sigma_max)), colorbar
ylabel('depth (m)');
title('Pacific zonal average DIP sensitivity to sigma (dDIP/dsigma)')
% deep
subplot(2,1,2) ;
contourf(grd.yt,-grd.zt(10:end),PZA(10:end,:),[S_sigma_min:0.01:S_sigma_max]);
set(gca,'color','black')
colormap(darkb2r(S_sigma_min,S_sigma_max)), colorbar
caxis([S_sigma_min, S_sigma_max])
colorbar
xlabel('latitude (deg)')
ylabel('depth (m)')
set(gcf, 'InvertHardcopy', 'off')
figTitle = 'DIP_sigma_PacificZonal';
print(gcf,[figDir 'SensiTest_' figTitle '.png'],'-dpng')
% exportfig(gcf,'Figs91/pza_dip','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')


% make a zonal average of age for the Atlantic basin
nfig = nfig + 1 ;
figure(nfig)

ATL = MSKS.ATL;
AZA = squeeze(nansum(ATL.*S_sigma.*dVt,2)./sum(ATL.*dVt,2))';
% surface
subplot(2,1,1) ;
contourf(grd.yt,-grd.zt(1:9),AZA(1:9,:),[S_sigma_min:0.01:S_sigma_max]);
set(gca,'color','black')
colormap(darkb2r(S_sigma_min, S_sigma_max)), colorbar
ylabel('depth (m)')
title('Atlantic zonal average DIP sigma sensitivity (dDIP/dsigma)')
% deep
subplot(2,1,2) ;
contourf(grd.yt,-grd.zt(10:end),AZA(10:end,:),[S_sigma_min:0.01:S_sigma_max]);
set(gca,'color','black')
colormap(darkb2r(S_sigma_min, S_sigma_max)), colorbar
xlabel('latitude (deg)')
ylabel('depth (m)')
set(gcf, 'InvertHardcopy', 'off')
figTitle = 'DIP_sigma_AtlanticZonal';
print(gcf,[figDir 'SensiTest_' figTitle '.png'],'-dpng')
% exportfig(gcf,'Figs91/aza_dip','fontmode','fixed','fontsize',12, ...
%           'color','rgb','renderer','painters')

% make surface sensitivity map
nfig = nfig + 1 ;
figure(nfig)
pcolor(nanmean(S_sigma(:,:,1:3),3));colorbar;shading flat
set(gca,'color','black')
colormap(darkb2r(S_sigma_min, S_sigma_max)), colorbar
set(gcf, 'InvertHardcopy', 'off')
figTitle = 'DIP_sigma_surf';
print(gcf,[figDir 'SensiTest_' figTitle '.png'],'-dpng')
% exportfig(gcf,'Figs91/surface_dip_anomaly','fontmode','fixed','fontsize',12, ...
%           'color','rgb','renderer','painters')

%}



%%% ------------------ CARBON MODEL -----------------------------------
%{
if par.Cmodel == on
	%--- Load Carbon Model Output ----
	DIC = model.DIC ;
	POC = model.POC ;
	DOC = model.DOC ;
	PIC = model.PIC ;
	ALK = model.ALK ;

	%--- Split Px into DIP, DOP, and POP
	DICx = model.Cx(0*nwet+1:1*nwet,:) ;
    POCx = model.Cx(1*nwet+1:2*nwet,:) ;
    DOCx = model.Cx(2*nwet+1:3*nwet,:) ;
    PICx = model.Cx(3*nwet+1:4*nwet,:) ;
    ALKx = model.Cx(4*nwet+1:5*nwet,:) ;

end
%}

% -------------------- old make plots --------------------
% nfig = 0 ;
% nfig = nfig + 1 ;
% figure(nfig)
% dDIP = pDIP - DIP ;
% % make a zonal cross section of the age
% contourf(grd.yt,-grd.zt,squeeze(dDIP(:,170,:))')
% set(gca,'color','black')
% colormap(darkb2r(-0.015,0.025)), colorbar
% xlabel('latitude (deg)');
% ylabel('depth (m)')
% t = sprintf('DIP anomally x = %4.1f deg', 170);
% title(t);
%
% nfig = nfig + 1 ;
% figure(nfig)
% % make a zonal average of age for the Pacific basin
% PAC = MSKS.PAC;
% PZA = squeeze(nansum(PAC.*dDIP.*dVt,2)./sum(PAC.*dVt,2))';
% subplot(2,1,1) ;
% contourf(grd.yt,-grd.zt(1:9),PZA(1:9,:),[-0.015:0.005:0.025]);
% set(gca,'color','black')
% colormap(darkb2r(-0.015,0.025)), colorbar
% ylabel('depth (m)');
% title('Pacific zonal average DIP anomaly')
% %
% subplot(2,1,2) ;
% contourf(grd.yt,-grd.zt(10:end),PZA(10:end,:),[-0.015:0.005:0.025]);
% set(gca,'color','black')
% colormap(darkb2r(-0.015,0.025)), colorbar
% caxis([-0.015 0.025])
% colorbar
% xlabel('latitutde (deg)')
% ylabel('depth (m)')
% set(gcf, 'InvertHardcopy', 'off')
% exportfig(gcf,'Figs91/pza_dip','fontmode','fixed','fontsize',12, ...
%           'color','rgb','renderer','painters')
% %

% --- old Cmodel plots ------
%{
if Cmodel == on
    if isfile(pfname)
        pDIC = pdata.DIC ;
        pDOC = pdata.DOC ;
        pPOC = pdata.POC ;
        pPIC = pdata.PIC ;
    else
        DIC = data.DIC ;
        POC = data.POC ;
        DOC = data.DOC ;
        PIC = data.PIC ;

        [par, C ] = eqCcycle(p0, par) ;
        pDIC = M3d+nan ;  pDIC(iwet) = C(0*nwet+1:1*nwet) ;
        pPOC = M3d+nan ;  pPOC(iwet) = C(1*nwet+1:2*nwet) ;
        pDOC = M3d+nan ;  pDOC(iwet) = C(2*nwet+1:3*nwet) ;
        pPIC = M3d+nan ;  pPIC(iwet) = C(3*nwet+1:4*nwet) ;
        pALK = M3d+nan ;  pALK(iwet) = C(4*nwet+1:5*nwet) ;
        par.DIC  = pDIC(iwet) ;
        par.DOC  = pDOC(iwet) ;
    end

    dDIC = pDIC + par.dicant - data.DIC ;
    nfig = nfig + 1 ;
    figure(nfig)
    % make a zonal cross section of the age
    contourf(grd.yt,-grd.zt,squeeze(dDIC(:,170,:))')
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    xlabel('latitude (deg)');
    ylabel('depth (m)')
    t = sprintf('DIC anomaly x = %4.1f deg', 170);
    title(t);

    nfig = nfig + 1 ;
    figure(nfig)
    % make a zonal average of age for the Pacific basin
    PAC = MSKS.PAC;
    PZA = squeeze(nansum(PAC.*dDIC.*dVt,2)./sum(PAC.*dVt,2))';
    subplot(2,1,1) ;
    contourf(grd.yt,-grd.zt(1:9),PZA(1:9,:),[-15:2:15]);
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    ylabel('depth (m)');
    title('Pacific zonal average DIP anomaly')
    %
    subplot(2,1,2) ;
    contourf(grd.yt,-grd.zt(10:end),PZA(10:end,:),[-15:2:15]);
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    xlabel('latitutde (deg)');
    ylabel('depth (m)');
    set(gcf, 'InvertHardcopy', 'off')
    exportfig(gcf,'Figs91/pza_dic','fontmode','fixed','fontsize',12, ...
              'color','rgb','renderer','painters')

    nfig = nfig + 1 ;
    figure(nfig)
    % make a zonal average of age for the Pacific basin
    ATL = MSKS.ATL;
    AZA = squeeze(nansum(ATL.*dDIC.*dVt,2)./sum(ATL.*dVt,2))';
    subplot(2,1,1) ;
    contourf(grd.yt,-grd.zt(1:9),AZA(1:9,:),[-15:2:15]);
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    ylabel('depth (m)');
    title('Atlantic zonal average DIC anomaly')
    %
    subplot(2,1,2) ;
    contourf(grd.yt,-grd.zt(10:end),AZA(10:end,:),[-15:2:15]);
    set(gca,'color','black')
    colormap(darkb2r(-15, 15)), colorbar
    xlabel('latitutde (deg)');
    ylabel('depth (m)');
    set(gcf, 'InvertHardcopy', 'off')
    exportfig(gcf,'Figs91/aza_dic','fontmode','fixed','fontsize',12, ...
              'color','rgb','renderer','painters')

    nfig = nfig + 1 ;
    figure(nfig)
    pcolor(nanmean(dDIC(:,:,1:3),3));colorbar;shading flat
    set(gca,'color','black')
    colormap(darkb2r(-5, 20)), colorbar
    set(gcf, 'InvertHardcopy', 'off')
    exportfig(gcf,'Figs91/surface_dic_anomaly','fontmode','fixed','fontsize',12, ...
              'color','rgb','renderer','painters')

end

if Omodel == on
    O2 = data.O2 ;
    GO = real(O2(iwet)) + 1e-9*randn(par.nwet,1);
    [par, O2] = eqOcycle(p0, par) ;
end
%}
