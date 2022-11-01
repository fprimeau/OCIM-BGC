clc; clear all; close all
% global iter
% iter = 0 ;
on   = true  ;
off  = false ;
format long

%ver = datestr(now,'mmmdd');
%RunVer = 'optGM15_CTL_He_PC_DOC0.25_DOP0'
%RunVer = 'Tv4_PC_DOC0.25_DOP0v8'
%RunVer = 'testCellinit/q0f2k1r1g1_CTL_He_PCCell_DOC0.25_DOP0'
%RunVer = 'optC_GM15_CTL_He_PC_DOC0.25_DOP0'
RunVer = 'optC_Cellv2_CTL_He_PCCell_DOC0.25_DOP0'

GridVer  = 91  ;
operator = 'A' ;

%model output directory
%outputDir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/', GridVer);
outputDir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/')
%figDir = strcat(outputDir,'FIGS_PCv9_DOC0.25_DOP0/');
%figDir = strcat(outputDir,'testCellinit/FIGS_testCellinit/q0f2k1r1g1_');
%figDir = strcat(outputDir,'FIGS_optGM15/');
figDir = strcat(outputDir,'FIGS_optC_Cell/');
outPath = figDir;

% load model output fields
fname = strcat(outputDir, RunVer, '.mat');
par.fname = fname;
load(fname);
model = data;

% load optimal parameter values
fxhat = strcat(outputDir, RunVer,'_xhat.mat');
par.fxhat = fxhat;
par.fxhatload = fxhat; % to make sure all non-optimized parameters are same as during the run.
load(fxhat);

par.Cmodel  = on ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function
%par.cscale  = 0.0 ;
par.LoadOpt = on ; % if load optimial par.
par.dynamicP = off ;


%-------------load data and set up parameters---------------------
SetUp ;
format short
xhat

%{
% P model parameters
par.opt_sigma = on ;
par.opt_kP_T  = on ;
par.opt_kdP   = on ;
par.opt_bP_T  = on ;
par.opt_bP    = on ;
par.opt_alpha = on ;
par.opt_beta  = on ;
% C model parameters
par.opt_bC_T  = on ;
par.opt_bC    = on ;
par.opt_d     = on ;
par.opt_kC_T  = on ;
par.opt_kdC   = on ;
par.opt_R_Si  = on ;
par.opt_rR    = on ;
par.opt_cc    = off ;
par.opt_dd    = off ;
% O model parameters
par.opt_O2C_T = off ;
par.opt_rO2C  = on ;
par.opt_O2P_T = off ;
par.opt_rO2P  = on ;
% Si model parameters
par.opt_dsi   = on  ;
par.opt_at    = off ;
par.opt_bt    = on  ;
par.opt_aa    = on  ;
par.opt_bb    = on  ;
%Trait Model parameters
par.opt_Q10Photo     = on ;
par.opt_fStorage     = on;
par.opt_PLip_PCutoff = on;
par.opt_PLip_scale   = off;
par.opt_PStor_rCutoff = on;
par.opt_PStor_scale  = off;
par.opt_alphaS       = on;
par.opt_fRibE 	     = on;
par.opt_kST0 	     = off;
%
%-------------load data and set up parameters---------------------
SetUp ;

% save results
% ATTENTION: Change this direcrtory to where you wanna save your output files
if ismac
    output_dir = sprintf('~/Documents/CP-model/MSK%2d/',GridVer);
elseif isunix
	output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/', GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/Cexp/']);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/TempSensi/' ...
    %                    'MSK%2d/PME4DICALK/'],GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/' ...
                        % 'TempSensi/MSK91/Zscore/'], GridVer);
    % output_dir = sprintf(['/DFS-L/DATA/primeau/weilewang/COP4WWF/' ...
                        % 'MSK%2d/'],GridVer);
	fig_dir = strcat(output_dir,'FIGS_PCCellv3b_DOC0.25_DOP0/');
	output_dir = sprintf('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK%2d/v3b_duplicate/', GridVer); %temporary
end
VER = strcat(output_dir,TRdivVer);
catDOC = sprintf('_DOC%0.2g_DOP%0.2g',par.cscale,par.pscale); % used to add scale factors to file names
% Creat output file names based on which model(s) is(are) optimized
if (par.Cmodel == off & par.Omodel == off & par.Simodel == off & par.Cellmodel == off)
	fname = strcat(VER,'_P');
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off & par.Cellmodel == off)
	base_name = strcat(VER,'_PC');
	fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == off)
	base_name = strcat(VER,'_PCO');
	fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on & par.Cellmodel == off)
	base_name = strcat(VER,'_PCSi');
	fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == off)
	base_name = strcat(VER,'_PCOSi');
	fname = strcat(base_name,catDOC);
elseif (par.Cmodel == off & par.Omodel == off & par.Simodel == off & par.Cellmodel == on) % cell model does nothing if C model is not on, so this case =Ponly
	base_name = strcat(VER,'_PCell');
	fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off & par.Cellmodel == on)
	base_name = strcat(VER,'_PCCellv3b');
	fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off & par.Cellmodel == on)
	base_name = strcat(VER,'_PCOCell');
	fname = strcat(base_name,catDOC);
elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on & par.Cellmodel == on)
	base_name = strcat(VER,'_PCOSiCell');
	fname = strcat(base_name,catDOC);
end

% if (par.Cmodel == off & par.Omodel == off & par.Simodel == off)
%     fname = strcat(VER,'_P');
% elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == off)
%     base_name = strcat(VER,'_PCv1');
%     catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
%     fname = strcat(base_name,catDOC);
% elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == off)
%     base_name = strcat(VER,'_PCOv6');
%     catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
%     fname = strcat(base_name,catDOC);
% elseif (par.Cmodel == on & par.Omodel == off & par.Simodel == on)
%     base_name = strcat(VER,'_PCSi');
%     catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
%     fname = strcat(base_name,catDOC);
% elseif (par.Cmodel == on & par.Omodel == on & par.Simodel == on)
%     base_name = strcat(VER,'_PCOSi');
%     catDOC = sprintf('_DOC%2.0e_DOP%2.0e',par.cscale,par.pscale);
%     fname = strcat(base_name,catDOC);
% end
par.fname = strcat(fname,'.mat') ;
% load optimal parameters if they exist
fxhat     = strcat(fname,'_xhat.mat');
par.fxhat = fxhat ;
load(par.fname) ;
load(par.fxhat) ;
%}
%--------------------- prepare parameters ------------------
% load optimal parameters from a file or set them to default values
par = SetPar(par) ;

%	% pack parameters into an array, assign them corresponding indices.
%	par = PackPar(par) ;

%------------------ extract parameters ---------------------------
% POP disolution constant [s^-1];
sigma = par.sigma ;
% linear parameter of npp to DIP assimilation function.
% par.alpha ;
% exponential parameter of npp to DIN assimilation function.
% par.beta ;

%------------------ extract data ---------------------------------
DIP  = model.DIP(iwet) ;
PO4  = par.po4obs(iwet)   ;
TRdiv= par.TRdiv      ;
I    = par.I          ;

if par.Cmodel == on
	DIC  = model.DIC(iwet) ;
	POC  = model.POC(iwet) ;
	DOC  = model.DOC(iwet) ;
end

% ----------------------------------------------
nn = par.nzo ; %number fo verticle boxes in euphotic zone / export depth

% -------------- C:P uptake ratio --------------------------------
W = d0(dVt(iwet)) ;
dVtwet = M3d*nan;
dVtwet(iwet) = dVt(iwet);
Wiprod = dVtwet(:,:,1:nn)/nansum(dVtwet(:,:,1:nn),'all');

C2P3D = M3d + nan ;
if par.Cellmodel==on
	C2P3D(iwet) = model.CellOut.C2P(iwet); %zero beneath the surface  layers
elseif par.Cmodel ==on
	C2P3D(iwet) = 1./(par.cc*PO4 + par.dd) ;  % DIP or PO4?
	%C2P3D(iwet) = 1./(7e-4*PO4 + 5e-3); % WL
	%C2P3D(iwet) = 1000./(6.6*PO4 + 5.3); %Qian
else
	C2P3D = M3d +nan;
	C2P3D(iwet) = 106; 		% redfield C:P
end
C2Pavg = nansum(C2P3D(:,:,1:nn).*Wiprod,'all');
fprintf('Average C:P uptake in Euphotic Layers is %4.2f \n',C2Pavg);


%------------------ prepare NPP for the model --------------------
%use cell model C2P instead of p2c function used in SetUp to define NPP
% not sure if this is a good idea or should be removed.
% if par.Cellmodel==on
% 	par.npp1  = (0.5*par.npp./grd.dzt(1)).*(1./model.CellOut.C2P(:,:,1)) ;
% 	par.npp2  = (0.5*par.npp./grd.dzt(2)).*(1./model.CellOut.C2P(:,:,2)) ;
% end

% Recall from SetUp, we defined:
% par.p2c = 0.006+0.0069*po4obs ;
% par.npp   = npp/(12*spd) ;		% units: mmol C/m^2/s
% par.npp1  = (0.5*par.npp./grd.dzt(1)).*par.p2c(:,:,1) ; % units: mmol P/m^2/s
% par.npp2  = (0.5*par.npp./grd.dzt(2)).*par.p2c(:,:,2) ;
% par.Lambda = M3d*0 ;
% par.Lambda(:,:,1) = 1./(1e-6+po4obs(:,:,1)) ;  %units: m^3/mmolP
% par.Lambda(:,:,2) = 1./(1e-6+po4obs(:,:,2)) ;

% DIP assimilation
LAM        = 0*M3d;
LAM(:,:,1) = (par.npp1.^par.beta).*par.Lambda(:,:,1);
LAM(:,:,2) = (par.npp2.^par.beta).*par.Lambda(:,:,2);
L          = d0(LAM(iwet));  % PO4 assimilation rate [s^-1];

%--------------- calculate primary production --------------------
G        = M3d*0        ;
G(iwet)  = par.alpha*L*DIP  ; % primary production [unit: mmol P/m^3/s]

% inegG = find(G<0); % should negative production values be removed?
% G(inegG)=nan;

Int_PNPP = G(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31;
Int_CNPP = G(:,:,1:nn).*grd.DZT3d(:,:,1:nn).*C2P3D(:,:,1:nn)*12;

PNPP = Int_PNPP*spa*1e-3 ; % convert to g P/m^2/yr
CNPP = Int_CNPP*spa*1e-3 ; % convert production from mg C/m^3/s to gC/m^2/year;
tem_PNPP = PNPP.*dAt(:,:,1:nn)*1e-15 ;
tem_CNPP = CNPP.*dAt(:,:,1:nn)*1e-15 ;
Sum_CNPP = nansum(tem_CNPP(:))    ;
fprintf('Model NPP (P) is %3.3e Pg P/yr \n\n',nansum(tem_PNPP(:))) ; %Pg/yr
fprintf('Model NPP is %3.3e Pg C/yr \n\n',Sum_CNPP) ; %Pg/yr

clear tem_PNPP tem_CNPP


% prod weighted C:P
WeightNPP = CNPP.*dAt(:,:,1:nn)/nansum(CNPP.*dAt(:,:,1:nn),'all'); %gC/yr
C2Pavg = nansum(C2P3D(:,:,1:nn).*WeightNPP,'all');
fprintf('Average C:P uptake in Euphotic Layers (NPP weighted) is %4.2f \n',C2Pavg);


%----------plot and save CNPP-----------------------

lat = grd.yt; lon = grd.xt;


GPtmp = sum(Int_PNPP(:,:,1:2)/31*spa*1000,3); % [mol/m^2/yr]
GPtmp(GPtmp==0)=nan;
GPtmp_latavg = mean(GPtmp,[2],'omitnan');
GPtmp_latstd = std(GPtmp,0,[2],'omitnan');
ind1 = ~isnan(GPtmp_latavg);

figure; hold on
fill([lat(ind1),fliplr(lat(ind1))],[(GPtmp_latavg(ind1)-GPtmp_latstd(ind1))', fliplr((GPtmp_latavg(ind1)+GPtmp_latstd(ind1))')],'r','LineStyle','none'); alpha(0.1); hold on
plot(lat,GPtmp_latavg,'-ro'); hold on
%plot(lat,C2P_latmedian1,'-.b','linewidth',2);
xlabel('Latitude')
ylabel('G [mol P/m^2/yr]')
title('Model zonal mean productivity')
axis tight; grid on;
figTitle = 'G_lat_avg';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')
clear h;


% figure;
% contourf(grd.xt,grd.yt,PNPP/31); c = colorbar;
% colormap(flipud(summer))
% title('Model DIP uptake rate','Fontsize',18);
% xlabel('Longitude');
% ylabel('Latitude');
% ylabel(c,'PNPP [molP/m^2/yr]');
% figTitle = 'Int_PNPP';
% print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
%
% figure;
% contourf(grd.xt,grd.yt,CNPP/12); c = colorbar;
% colormap(flipud(summer))
% caxis([0 40])
% title('Model NPP','Fontsize',18);
% xlabel('Longitude');
% ylabel('Latitude');
% ylabel(c,'NPP [molC/m^2/yr]');
% figTitle = 'Int_CNPP';
% print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
% %save([figDir 'Int_CNPP.mat'],'CNPP_surface','CNPP_Z2','DIPsurf')

%Calculate Cell model net primary production
%
% if par.Cellmodel ==on
% 	%uptake rate of P using the cell model predicted growth rate
% 	Mu = model.CellOut.mu/3600; 	% converts units to s^-1
% 	%Gcell = Mu.*model.DIP(:,:,1:nn);
% 	Gcell = Mu.*model.POP(:,:,1:nn);
%
% 	Int_PNPP = Gcell(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31;
% 	Int_CNPP = Gcell(:,:,1:nn).*grd.DZT3d(:,:,1:nn).*C2P3D(:,:,1:nn)*12;
%
% 	PNPPcell = Int_PNPPcell*spa*1e-3 ; % units: g P/m^2/yr
% 	CNPPcell = Int_CNPPcell*spa*1e-3 ; % convert production from mg C/m^3/s to g
% 	                           % C/m^2/year;
%
% 	%plot cell CNPP
% 	figure;
% 	contourf(grd.xt,grd.yt,PNPPcell/31); c = colorbar;
% 	colormap(flipud(summer))
% 	title('POP growth from cell model growth rate','Fontsize',18);
% 	xlabel('Longitude');
% 	ylabel('Latitude');
% 	ylabel(c,'PNPP [molP/m^2/yr]');
% 	figTitle = 'Int_PNPPcell';
% 	print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
%
% 	figure;
% 	CNPPcellmol = CNPPcell/12;
% 	imAlpha = ones(size(CNPPcellmol));
% 	imAlpha(isnan(CNPPcellmol)) =0;
% 	imagesc(grd.xt,grd.yt,CNPPcellmol,'AlphaData',imAlpha); hold on
% 	%contourf(grd.xt,grd.yt,CNPPcell/12);
% 	c = colorbar;
% 	colormap(flipud(summer))
% 	[CC,hh] = contour(grd.xt,grd.yt,CNPPcellmol,[-10:10:60],'k');
% 	clabel(CC,hh,'FontName','Times','Fontsize',8);
% 	axis xy
% 	title('Model NPP from POP*cell model growth rate','Fontsize',18);
% 	xlabel('Longitude');
% 	ylabel('Latitude');
% 	ylabel(c,'NPP [molC/m^2/yr]');
% 	figTitle = 'Int_CNPPcell';
% 	print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
%
% 	%tem_CNPP = CNPPcell.*dAt(:,:,1)*1e-15 ;
% 	%Sum_CNPP = nansum(tem_CNPP(:))    ;
% 	%fprintf('Cell model NPP (growth*POP*C2P) is %3.3e Pg C/yr \n',Sum_CNPP) ; %Pg/yr
%
% end

%---------------- calculate phosphorus export --------------------
PFD = buildPFD(par, 'POP') ;

F_diag_p = inv(W)*PFD'*W   ;
T_diag   = inv(W)*TRdiv'*W ;

junk = M3d ;
junk(:,:,1:nn) = 0 ;
Omega = junk(iwet) ;
Prod  = G(iwet)    ;
% adjoint method.
kP    = d0(par.kP_T * par.Tz + par.kdP) ;
Jex_P = d0(kP*Prod)*(sigma*I+par.kappa_p*(1-sigma) * ...
                     inv(F_diag_p+par.kappa_p*I))*((T_diag+kP)\Omega);

P3d = M3d+nan;
P3d(iwet) = Jex_P;

%fprintf('Integrated over the Euphotic Zone (n=%i): \n', nn)
% convert export from mmol P/m^3/s to mg P/m^2/day;
TOPexp3d = P3d(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*31*spd;
tem_Pexp = TOPexp3d.*dAt(:,:,1:nn);					% [mg P/day]
Sum_Pexp = nansum(tem_Pexp(:))*365*1e-18;			% [Pg P/yr]
fprintf('Model TOP export is %3.3e Pg P /yr   (Integrated to %4.1f m) \n',Sum_Pexp, sum(grd.dzt(1:nn)));
TOPexp = sum(TOPexp3d,3,'omitnan');

%convert export from mmol P/m^3/s to mol P/m^2/yr;
EXPORT.TOPexp3d = P3d(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
EXPORT.TOPexp = sum(EXPORT.TOPexp3d,3); %sum(TOPexp3d,3,'omitnan');

% POP export
[~,Gout] = buildPFD(par,'POP') ;
w = -Gout.w ;
% POPexp   = model.POP(:,:,nn).*w(:,:,nn)*31*spd ;
% tem_POPexp = POPexp.*dAt(:,:,nn);
% Sum_POPexp = nansum(tem_POPexp(:))*365*1e-18;
% fprintf('old Model POP export is %3.3e Pg P /yr   (at %4.1f m) \n\n',Sum_POPexp, grd.zt(nn));

%POP export: integrte POP beneath the euphotic zone
POPexp = nansum(par.kappa_p*model.POP(:,:,3:end).*dVt(:,:,3:end),3)*31*spd;
Sum_POPexp =nansum(POPexp(:))*365*1e-18;
fprintf('Model POP export is %3.3e Pg P /yr  (beneath EZ) \n\n',Sum_POPexp);

%integrated DOP remineralization below the euphotic zone. (should equal the TOP export calculated by the adjoint method)
DOPexpint = par.kdP*model.DOP(:,:,3:end).*grd.DZT3d(:,:,3:end)*31*spd;
tem_DOPexpint = nansum(DOPexpint.*dAt(:,:,3:end),3);
Sum_DOPexpint =nansum(tem_DOPexpint(:))*365*1e-18;
fprintf('Model integrated DOP below the Euphotic zone is %3.3e Pg P /yr \n\n',Sum_DOPexpint);
DOPexpint = sum(DOPexpint,3,'omitnan');

%---------------- calculate carbon export -------------------------
if par.Cmodel ==on
	PFD = buildPFD(par, 'POC') ;

	F_diag_p = inv(W)*PFD'*W   ;
	T_diag   = inv(W)*TRdiv'*W ;

	junk = M3d ;
	junk(:,:,1:nn) = 0 ;
	Omega = junk(iwet) ;
	Prod  = G(iwet).*C2P3D(iwet) ;
	% adjoint method.
	kC    = d0(par.kC_T * par.Tz + par.kdC) ;
	Jex_C = d0(kC*Prod)*(sigma*I+par.kappa_p*(1-sigma) * ...
	                     inv(F_diag_p+par.kappa_p*I))*((T_diag+kC)\Omega);

	C3d = M3d+nan ;
	C3d(iwet) = Jex_C ;

	% kC*DOC * volume weight ->integrated over z3:z24
	% does this equal DOC export? (TOCexp - DOCexp)

	%integrate over the Euphotic Zone
	% convert export from mmol C/m^3/s to mg C/m^2/day;
	TOCexp3d = C3d(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*12*spd;
	tem_Cexp = TOCexp3d.*dAt(:,:,1:nn);
	Sum_Cexp = nansum(tem_Cexp(:))*365*1e-18;
	fprintf('Model TOC export is %3.3e Pg C /yr   (Integrated to %4.1f m) \n',Sum_Cexp, sum(grd.dzt(1:nn)));
	TOCexp = sum(TOCexp3d,3,'omitnan');

	EXPORT.TOCexp3d = C3d(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*spa/1000; %[mol P/m^2/yr]
	EXPORT.TOCexp = sum(EXPORT.TOCexp3d,3); %sum(TOPexp3d,3,'omitnan');

	% inegTOCexp = find(TOCexp<0);
	% TOCexp(inegTOCexp) = nan;

	% POC export
	[~,Gout] = buildPFD(par,'POC') ;
	w = -Gout.w ;
	% POCexp   = model.POC(:,:,nn).*w(:,:,nn)*12*spd ;  %[mmol/m^3]*[m/s]*[mg/mmol]*[s/day] = [mg/m^/day]
	% tem_POCexp = POCexp.*dAt(:,:,nn);
	% Sum_POCexp = nansum(tem_POCexp(:))*365*1e-18;
	% fprintf('old Model POC export is %3.3e Pg C /yr   (at %4.1f m) \n\n',Sum_POCexp, grd.zt(nn));

	% POC export: integrate POC beneath the euphotic zone
	%POCexp = nansum(par.kappa_p*model.POC(:,:,3:end).*dVt(:,:,3:end),3)*12*spd;  % [1/s]*[mmol/m^3]*[m^3]*[mg/mmol]*[s/day] =  [mg/day]
	POCexp = par.kappa_p*model.POC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_POCexp = nansum(POCexp.*dAt(:,:,3:end),3);
	Sum_POCexp =nansum(tem_POCexp(:))*365*1e-18;
	fprintf('Model POC export is %3.3e Pg C /yr  (beneath EZ) \n\n',Sum_POCexp);
	POCexp = sum(POCexp,3,'omitnan');

	%integrated DOC remineralization below the euphotic zone. (should equal the TOC export calculated by the adjoint method)
	DOCexpint = par.kdC*model.DOC(:,:,3:end).*grd.DZT3d(:,:,3:end)*12*spd;
	tem_DOCexpint = nansum(DOCexpint.*dAt(:,:,3:end),3);
	Sum_DOCexpint =nansum(tem_DOCexpint(:))*365*1e-18;
	fprintf('Model integrated DOC below the Euphotic zone is %3.3e Pg C /yr \n\n',Sum_DOCexpint);
	DOCexpint = sum(DOCexpint,3,'omitnan');


% Save export
%% save
save([outPath 'EXPORTopt.mat'],'EXPORT');


%{
	[~,Gout] = buildPFD(par,'POC') ;
	w = -Gout.w ;
	POCexp   = data.POC(:,:,nn).*w(:,:,nn+1)*12*spd ;
	tem_POCexp = POCexp.*dAt(:,:,nn);
	Sum_POCexp = nansum(tem_POCexp(:))*365*1e-18;
	fprintf('Model POC export is %3.3e Pg C /yr   (at %4.1f m) \n\n',Sum_POCexp, sum(grd.dzt(1:nn)) );
%}

	% ----- calculate POC flux
	% fPOC = w(:,:,2:25).*POC*spd*12 ; % POC flux (mg/m^2/day)

	% ----------- C:P export ratio -------------------------
	%volume weight
	Wiexp = dVtwet(:,:,2)/nansum(dVtwet(:,:,2),'all');
	% TOC export weight
	WeightEXP =  (TOCexp/12)./(nansum(TOCexp/12,'all')); %CNPP.*dAt(:,:,1:nn)/nansum(CNPP.*dAt(:,:,1:nn),'all');

	% weight by TOC export
	% Wcexp = TOCexp./Sum_Cexp;
	% Wcexp = DOCexpint./Sum_DOCexpint;
	% mmolC/m^2/s / mmolP/m^2/s
	C2Pexp = (nansum(C3d(:,:,1:nn).*grd.DZT3d(:,:,1:nn),3))./(nansum(P3d(:,:,1:nn).*grd.DZT3d(:,:,1:nn),3));
	fprintf('Model average C:Pexport (TOCexport weighted adj.) is %4.2f \n',nansum(C2Pexp.*WeightEXP,'all'));


	% ratio of DOC remineralized below EZ to DOP remineralized below EZ
	% TOCexport/TOPexport
	C2Pexpint = (DOCexpint/12/1000*365)./(DOPexpint/31/1000*365);
	C2Pexpint_avg = nansum(C2Pexpint.*Wiexp,'all');
	fprintf('Model average C:Pexport (area weighted DOM remin below EZ) is %4.2f \n',C2Pexpint_avg);
	fprintf('Model average C:Pexport (TOCexport weighted DOM remin below EZ) is %4.2f \n',nansum(C2Pexpint.*WeightEXP,'all'));

	C2Pexp = M3d+nan;
	C2Pexp = C3d./P3d;
	%C2Pexp(C2Pexp<0) = nan;

	% weighted export. should be about 106
	C2Pexp_avg = nansum(C2Pexp(:,:,2).*Wiexp,'all');
	fprintf('Model average C:Pexport (volume weighted  adj.) is %4.2f \n\n',C2Pexp_avg);

	% Plot C2Pexpint instead
	figure;
	%pcolor(grd.xt,grd.yt,C2Pexp(:,:,2));
	imAlpha = ones(size(C2Pexp(:,:,nn)));
	imAlpha(isnan(C2Pexp(:,:,nn))) =0;
	imagesc(grd.xt,grd.yt,C2Pexp(:,:,nn),'AlphaData',imAlpha); hold on
	cb=colorbar;
	colormap(flipud(summer))
	[CC,hh] = contour(grd.xt,grd.yt,C2Pexp(:,:,nn),[106 106],'k');
	clabel(CC,hh,'FontName','Times');
	axis xy
	title(['C:P export at Z=' num2str(grd.zt(nn))])
	%title('C:P export')
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(cb,'C2P [molC/molP]');
	figTitle = 'C2Pexp_z2';
	print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')


	% average across Latitude or basin?
	% plot all vertical profiles of C2Pexp
	% plot profiles within a basin

	% figure; hold on
	% for ii = 1:180
	% 	plot(squeeze(C2Pexp(30,ii,:)),grd.zt); hold on
	% end
	% title(['C:P export profiles at lat =' num2str(grd.yt(30))])
	% xlabel('C:P total organic matter export')
	% ylabel('depth')
	% figTitle = 'C2Pexport_profiles30';
	% print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')

	%C2Pexp = POCexp./POPexp;

	%keyboard;
	%%% plot export
	% this idoes not work yet
	%{
	figure;
	contourf(grd.xt,grd.yt,C2Pexp); c = colorbar;
	title('Model C:P export','Fontsize',18);
	xlabel('Longitude');
	ylabel('Latitude');
	ylabel(c,'C:P [gC/gP]');
	grid off

	figTitle = 'C2Pexport';
	print(gcf,[fig_dir 'FIG_' figTitle '.png'],'-dpng')
	%}

else
	%------------- for SONT -----------------------------------
	fprintf('Carbon model is off. Converting Pexp using constant C:P = 106 \n')
	% convert export from mmol P/m^3/s to mg C/m^2/day;
	TOCexp3d = P3d(:,:,1:nn).*C2P3D(:,:,1:nn).*grd.DZT3d(:,:,1:nn)*12*spd;
	tem_Cexp = TOCexp3d.*dAt(:,:,1:nn);
	Sum_Cexp = nansum(tem_Cexp(:))*365*1e-18;
	fprintf('Model TOC export is %3.3e Pg C /yr   (Integrated to %4.1f m) \n',Sum_Cexp, sum(grd.dzt(1:nn)));
	TOCexp = sum(TOCexp3d,3,'omitnan');

	POCexp   = data.POP(:,:,nn).*C2P3D(:,:,nn).*w(:,:,nn)*12*spd ;
	tem_POCexp = POCexp.*dAt(:,:,nn);
	Sum_POCexp = nansum(tem_POCexp(:))*365*1e-18;
	fprintf('Model POC export is %3.3e Pg C /yr   (at %4.1f m) \n\n',Sum_POCexp, grd.zt(nn));
end

%%------- Smooth export fields -------------------
DOCexp = TOCexp - POCexp; DOCexp(DOCexp(:)<0) = 0 ;

 % TOCexp = smoothit(grd,M3d,TOCexp,3,1e5);
 % POCexp = smoothit(grd,M3d,POCexp,3,1e5);
 % DOCexp = smoothit(grd,M3d,DOCexp,3,1e5);
 % str_smoothparams = 'n3r1e5';
 str_smoothparams = 'NoSmoothing';

%%------- Plot TOC export -------------------------

iwetsurf = find(M3d(:,:,1));
idrysurf = find(~M3d(:,:,1));
TOCexp(idrysurf) = nan;

TOCexp_std = std(TOCexp(iwetsurf)./12/1000*365,[],'all');
TOCexp_mean = nanmean(TOCexp(iwetsurf)./12/1000*365,'all');
fprintf('Mean TOCexp is %4.3f ; One standard deviation of TOCexp is %4.3f \n', TOCexp_mean, TOCexp_std);
[TOCmin, TOCmax] = bounds(TOCexp./12/1000*365,'all');
TOCmin = -4;
TOCmax = 8;
%TOCmin = TOCexp_mean-2*TOCexp_std ;
%TOCmax = TOCexp_mean+2*TOCexp_std;

figure
imAlpha = ones(size(TOCexp));
imAlpha(isnan(TOCexp)) =0;
imagesc(grd.xt,grd.yt,TOCexp/12/1000*365,'AlphaData',imAlpha); hold on
cb = colorbar; colormap(gca,darkb2r(TOCmin, TOCmax));
%colormap(gca,darkb2r(-1.4, 7));
%[CC,hh] = contour(grd.xt,grd.yt,TOCexp/12/1000*365,[0:2:8],'k');
%clabel(CC,hh,'FontName','Times','FontSize',6);
axis xy;
title('Model Total Organic Carbon Export')
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'TOC export [mol C/m^2/yr]');
figTitle = ['TOCexp_' str_smoothparams];
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')


%--------------- compare to ANCP -----------------------------

Lat_HOTS = 22+45/60; Lon_HOTS = mod(-158,360);
Lat_BATS = 31+40/60; Lon_BATS = mod((-64-10/60),360);
Lat_OSP  = 50+1/60;  Lon_OSP  = mod((-144-9/60),360);

indx_hots = length(find(grd.xt<Lon_HOTS));
indy_hots = length(find(grd.yt<Lat_HOTS));

indx_bats = length(find(grd.xt<Lon_BATS));
indy_bats = length(find(grd.yt<Lat_BATS));

indx_osp = length(find(grd.xt<Lon_OSP));
indy_osp = length(find(grd.yt<Lat_OSP));

% find ANCP at specific location and convert unit from
% mg/m2/day to mol/m2/year;
TOCexp_HOTS = TOCexp(indy_hots,indx_hots)/12/1000*365;
TOCexp_BATS = TOCexp(indy_bats,indx_bats)/12/1000*365;
TOCexp_OSP  = TOCexp(indy_osp,indx_osp)/12/1000*365;
fprintf('TOC export at HOT is %2.2f mol/m2/year\n', TOCexp_HOTS)
fprintf('TOC export at BATS is %2.2f mol/m2/year \n',TOCexp_BATS)
fprintf('TOC export at OSP is %2.2f mol/m2/year \n\n', TOCexp_OSP)

POCexp_HOTS = POCexp(indy_hots,indx_hots)/12/1000*365;
POCexp_BATS = POCexp(indy_bats,indx_bats)/12/1000*365;
POCexp_OSP  = POCexp(indy_osp,indx_osp)/12/1000*365;
fprintf('POC export at HOT is %2.2f mol/m2/year\n', POCexp_HOTS)
fprintf('POC export at BATS is %2.2f mol/m2/year \n',POCexp_BATS)
fprintf('POC export at OSP is %2.2f mol/m2/year \n\n', POCexp_OSP)

DOCexp_HOTS = DOCexp(indy_hots,indx_hots)/12/1000*365;
DOCexp_BATS = DOCexp(indy_bats,indx_bats)/12/1000*365;
DOCexp_OSP  = DOCexp(indy_osp,indx_osp)/12/1000*365;
fprintf('DOC export at HOT is %2.2f mol/m2/year\n', DOCexp_HOTS)
fprintf('DOC export at BATS is %2.2f mol/m2/year \n',DOCexp_BATS)
fprintf('DOC export at OSP is %2.2f mol/m2/year \n\n', DOCexp_OSP)

%% ----latitude regions --------------------------------

msk_tropical = M3d(:,:,1)*0;
msk_tropical(length(find(grd.yt<-15)):length(find(grd.yt<15)),:) = 1;

junk1 = M3d(:,:,1)*0;
junk1(length(find(grd.yt<-30)):length(find(grd.yt<30)),:) = 1;
msk_subtro = junk1-msk_tropical;

junk2  = M3d(:,:,1)*0;
junk2(length(find(grd.yt<-45)):length(find(grd.yt<45)),:) = 1;
msk_subtro_subpo = junk2-junk1;

junk3 =  M3d(:,:,1)*0;
junk3(length(find(grd.yt<-60)):length(find(grd.yt<60)),:) = 1;
msk_subpolar = junk3-junk2;

junk4 =  M3d(:,:,1)*0;
junk4(length(find(grd.yt<-88)):length(find(grd.yt<90)),:) = 1;
msk_polar = junk4-junk3;

% units mg/m2/day;
TOCexp_tropical = msk_tropical.*TOCexp;
TOCexp_subtro = msk_subtro.*TOCexp;
TOCexp_subtro_subpo = msk_subtro_subpo.*TOCexp;
TOCexp_subpolar = msk_subpolar.*TOCexp;
TOCexp_polar = msk_polar.*TOCexp;

TOCexp_tropical(TOCexp_tropical(:)==0) = nan;
TOCexp_subtro(TOCexp_subtro(:)==0) = nan;
TOCexp_subtro_subpo(TOCexp_subtro_subpo(:)==0) = nan;
TOCexp_subpolar(TOCexp_subpolar(:)==0) = nan;
TOCexp_polar(TOCexp_polar(:)==0) = nan;

mean_TOC_tropical = nanmean(TOCexp_tropical(:))/12/1000*365;
mean_TOC_subtro = nanmean(TOCexp_subtro(:))/12/1000*365;
mean_TOC_subtro_subpo = nanmean(TOCexp_subtro_subpo(:))/12/1000*365;
mean_TOC_subpolar = nanmean(TOCexp_subpolar(:))/12/1000*365;
mean_TOC_polar = nanmean(TOCexp_polar(:))/12/1000*365;
fprintf('TOC export at tropical is %2.2f mol/m2/yr\n', mean_TOC_tropical)
fprintf('TOC export at subtropical is %2.2f mol/m2/yr \n',mean_TOC_subtro)
fprintf('TOC export at subtropical-subpolar is %2.2f mol/m2/yr \n', mean_TOC_subtro_subpo)
fprintf('TOC export at subpolar is %2.2f mol/m2/year \n', mean_TOC_subpolar)
fprintf('TOC export at polar is %2.2f mol/m2/year \n\n', mean_TOC_polar)

TOCexp_tropical_tmp = TOCexp_tropical.*dAt(:,:,1);
Sum_Cexp_tropical = nansum(TOCexp_tropical_tmp(:))*365*1e-18;

TOCexp_subtro_tmp = TOCexp_subtro.*dAt(:,:,1);
Sum_Cexp_subtro = nansum(TOCexp_subtro_tmp(:))*365*1e-18;

TOCexp_subtro_subpo_tmp = TOCexp_subtro_subpo.*dAt(:,:,1);
Sum_Cexp_subtro_subpo = nansum(TOCexp_subtro_subpo_tmp(:))*365*1e-18;

TOCexp_subpolar_tmp = TOCexp_subpolar.*dAt(:,:,1);
Sum_Cexp_subpolar = nansum(TOCexp_subpolar_tmp(:))*365*1e-18;

TOCexp_polar_tmp = TOCexp_polar.*dAt(:,:,1);
Sum_Cexp_polar = nansum(TOCexp_polar_tmp(:))*365*1e-18;

fprintf('tropical TOC export is %4.1f%% of total TOC export \n',Sum_Cexp_tropical/Sum_Cexp*100)
fprintf('subtropical TOC export is %4.1f%% of total TOC export \n',Sum_Cexp_subtro/Sum_Cexp*100)
fprintf('subtropical-subpolar TOC export is %4.1f%% of total TOC export \n',Sum_Cexp_subtro_subpo/Sum_Cexp*100)
fprintf('subpolar TOC export is %4.1f%% of total TOC export \n',Sum_Cexp_subpolar/Sum_Cexp*100)
fprintf('polar TOC export is %4.1f%% of total TOC export \n \n',Sum_Cexp_polar/Sum_Cexp*100)


% DOC export units mol/m2/yr;
DOCexp_tropical = msk_tropical.*DOCexp;
DOCexp_subtro = msk_subtro.*DOCexp;
DOCexp_subtro_subpo = msk_subtro_subpo.*DOCexp;
DOCexp_subpolar = msk_subpolar.*DOCexp;
DOCexp_polar = msk_polar.*DOCexp;

DOCexp_tropical(DOCexp_tropical(:)==0) = nan;
DOCexp_subtro(DOCexp_subtro(:)==0) = nan;
DOCexp_subtro_subpo(DOCexp_subtro_subpo(:)==0) = nan;
DOCexp_subpolar(DOCexp_subpolar(:)==0) = nan;
DOCexp_polar(DOCexp_polar(:)==0) = nan;

mean_DOC_tropical = nanmean(DOCexp_tropical(:))/12/1000*365;
mean_DOC_subtro = nanmean(DOCexp_subtro(:))/12/1000*365;
mean_DOC_subtro_subpo = nanmean(DOCexp_subtro_subpo(:))/12/1000*365;
mean_DOC_subpolar = nanmean(DOCexp_subpolar(:))/12/1000*365;
mean_DOC_polar = nanmean(DOCexp_polar(:))/12/1000*365;
fprintf('DOC export at tropical is %2.2f mol/m2/yr\n', mean_DOC_tropical)
fprintf('DOC export at subtropical is %2.2f mol/m2/yr \n',mean_DOC_subtro)
fprintf('DOC export at subtropical-subpolar is %2.2f mol/m2/yr \n', mean_DOC_subtro_subpo)
fprintf('DOC export at subpolar is %2.2f mol/m2/year \n', mean_DOC_subpolar)
fprintf('DOC export at polar is %2.2f mol/m2/year \n\n', mean_DOC_polar)


% calculate DOC to TOC export ratio for the four biomes.
D2T = DOCexp./TOCexp;
D2T(isinf(D2T)) = nan;
D2T_tropical = (msk_tropical.*D2T);
D2T_subtro = (msk_subtro.*D2T);
D2T_subtro_subpo = (msk_subtro_subpo.*D2T);
D2T_subpolar = (msk_subpolar.*D2T);
D2T_polar = (msk_polar.*D2T);

D2T_tropical(D2T_tropical(:)==0) = nan;
D2T_subtro(D2T_subtro(:)==0) = nan;
D2T_subtro_subpo(D2T_subtro_subpo(:)==0) = nan;
D2T_subpolar(D2T_subpolar(:)==0) = nan;
D2T_polar(D2T_polar(:)==0) = nan;

mean_D2T_tropical = nanmean(D2T_tropical(:));
mean_D2T_subtro = nanmean(D2T_subtro(:));
mean_D2T_subtro_subpo = nanmean(D2T_subtro_subpo(:));
mean_D2T_subpolar = nanmean(D2T_subpolar(:));
mean_D2T_polar = nanmean(D2T_polar(:));

fprintf('tropical zonal mean DOC to TOC export ratio is %2.2f percent\n', ...
        mean_D2T_tropical*100)
fprintf('subtropical zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subtro*100)
fprintf('subtropical subpolar zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subtro_subpo*100)
fprintf('subpolar zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subpolar*100)
fprintf('polar zonal mean DOC to TOC export ratio is %2.2f percent\n\n', ...
        mean_D2T_polar*100)
