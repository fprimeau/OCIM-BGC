% Cell model driver script to make a map of total respiration quotients, 
% to be read in the OCIM biogeochemical model. 
% This script loads optimal cell model parameter values from the output 
% published in Sulivan et al., 2024 
% RQO2C is the total respiration quotient of reminerlization of organic matter, including nitrification
% Run CellCHNOP
clear all; close all
spd  = 24*60^2 ;
spa  = 365*spd ;
on   = true    ;
off  = false   ;


% load in optimal parameter values
% greenplanet path
% output_dir = '/DFS-L/DATA/primeau/meganrs/C2P-Model-published/output/'
% local path
output_dir = '../../../../C2P_Paper/Sullivan-C2P-Model-2023GB007986/output/';
fname_opt = 'optC_Cell_CTL_He_PCCell_DOC0.25_DOP0';
	par.fxhatload = strcat(output_dir,fname_opt,'_xhat.mat');
	par.fmodelload = strcat(output_dir,fname_opt,'.mat') ;


%% load SetUp data
% load T,N,P,and Light inputs for cell model 
addpath('../../DATA/BGC_24layer/')
load('OCIM2_CTL_He.mat','output') ;
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
    clear output 
load TS_WOA_91x180x24.mat tempobs % WOA temperature & salinity
load O2_Nut_WOA_91x180x24.mat DIN_obs DIP_obs % WOA O2 Si DIN DIP observations
load PARobs_processed_91x180x24.mat PARobs

par.no3obs = DIN_obs;
par.po4obs = DIP_obs;
par.PARobs = PARobs;
par.Temp = tempobs;


%% Run SetUp
addpath('../src/')
par.optim   = off ;
par.Cmodel  = off ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.C2P_Tzmodel = off;
par.C2P_PO4model = off; 
par.C2P_constant = off;
par.LoadOpt = on ; % if load optimial par.
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigP  = off ; 
par.opt_Q10P  = off ;
par.opt_kdP   = off ;
par.opt_bP_T  = off ; 
par.opt_bP    = off ;
par.opt_alpha = off ;
par.opt_beta  = off ;
% C model parameter
par.opt_sigC  = off ; 
par.opt_kru   = off ;
par.opt_krd   = off ;
par.opt_etau  = off ;
par.opt_etad  = off ;
par.opt_bC_T  = off ;
par.opt_bC    = off ; 
par.opt_d     = off ;
par.opt_Q10C  = off ;
par.opt_kdC   = off ; 
par.opt_R_Si  = off ; 
par.opt_rR    = off ; 
% --- C:P function parameters -----
% phosphate-dependent function parameters
par.opt_cc    = off ;
par.opt_dd    = off ; 
% temperature-dependent function parameters
par.opt_ccT   = off; 
par.opt_ddT   = off;
%Cell model parameters
par.opt_Q10Photo = off;
par.opt_fRibE 	 = off;
par.opt_fStorage = off;
par.opt_kST0 	 = off;
par.opt_PStor_rCutoff = off;
par.opt_PStor_scale = off;
par.opt_PLip_PCutoff = off;
par.opt_PLip_scale = off;
par.opt_alphaS = off;
par.opt_gammaDNA = off;

%% SetUp
%SetUpCell

% run set parameter values
par = SetPar(par)

% PackPar
%par = PackPar(par)
%x0 = par.p0;
x0 = [];

%% run cell model
x=x0;
% ii = 10
% P0 = inputs.P(:); N0 = inputs.N(:);
% T0 = inputs.T(:); Irr0 = inputs.PAR(:);

iprod = find(M3d(:,:,1:2)); %production in top two layers
N0 = par.no3obs(iprod)./10^6;   % *par.permil ; convert [umol/kg --> mmol/m^3 --> mol/L]
T0 = par.Temp(iprod);
Irr0 = par.PARobs(iprod);
P0 = par.po4obs(iprod)./10^6;	% *par.permil ; convert [umol/kg --> mmol/m^3 --> mol/L]


[CellOut, parBIO] = CellCHNOP(par,x, P0,N0,T0,Irr0);
fprintf('done  \n')

% regrid output
    CellO2C = M3d*NaN;
    CellO2C(iprod) = CellOut.RQtotalO2toC;

%% average resporation quotient in the top 2 layers in the euphotic zone; givee all veritcal layers this map;   

%average 
O2Csurf = mean(CellO2C(:,:,1:2),3,'omitnan');

% extend verically 
O2C = M3d*NaN; 
for ii = 1:length(grd.dzt)
    O2C(:,:,ii) = O2Csurf; 
end


%% save file to DATA
CellO2C_parameters = parBIO;
save('../../DATA/BGC_24layer/CellO2C_91x180x24.mat','O2C','CellO2C_parameters');


%% quick check output visually
%%% plot surface RQ 
%%% surface layer only
% lon = grd.xt;
% lat = grd.yt;
% RQ1 = O2C(:,:,1);
% [RQmin RQmax] = bounds(RQ1,'all')
% clevs = [1.35:0.05:1.50];
% 
% figure(1);
% imAlpha = ~isnan(RQ1);
% imagesc(lon,lat,RQ1,'AlphaData',imAlpha)
% axis xy; hold on
% cb = colorbar;
% xlabel('Longitude');
% ylabel('Latitude');
% ylabel(cb,'Respiration Quotient [mol O_2/mol C]','Rotation',-90,'VerticalAlignment','bottom');
% title('Cell Model surface total RQ (-O2:C)','Fontsize',18);
% % add some contours labels
% [CC,hh] = contour(lon,lat,RQ1,clevs,'k');
% clabel(CC,hh,'FontName','Times');
