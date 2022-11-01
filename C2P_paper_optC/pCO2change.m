% pCO2change
% impact on atmospheric CO2

% Integrate total DIC in ocean in each model. See change between optimized
% and future scenarios in each of the models

%%
clc; clear all; close all
on   = true  ;
off  = false ;
spd  = 24*60^2 ;
spa  = 365*spd ;

GridVer  = 91  ;


%% ---------- load model output --------------------------
%% load model output
outputDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/C2P_paper_optC/';
figDir = strcat(outputDir,'FIGS_CellGM15Tz/');

% load cell model output fields
fname = strcat(outputDir,'optC_Cell_CTL_He_PCCell_DOC0.25_DOP0.mat');
load(fname);
model_cell = data;
% load optimal parameter values
fxhat = strcat(outputDir,'optC_Cell_CTL_He_PCCell_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_cell = xhat.allparams;
clear data xhat

% load Tup 5C cell model output fields
fname_future = strcat(outputDir, 'futureproj/Tup5Cv2_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_cell_Tup5C = d.model;
clear d

% load P04obs down 20 percent cell model output fields
fname_future = strcat(outputDir, 'futureproj/Pdn20pct_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_cell_Pdn20 = d.model;
clear d

% load CMIP5 projection to 2100 cell model output fields
fname_future = strcat(outputDir, 'futureproj/CMIP2100_Cell_CTL_He_PCCell.mat');
d = load(fname_future);
model_cell_CMIP = d.model;
clear d

% load GM15 model output
fname = strcat(outputDir,'optC_GM15_CTL_He_PC_DOC0.25_DOP0.mat');
load(fname);
model_GM15 = data;
fxhat = strcat(outputDir,'optC_GM15_CTL_He_PC_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_GM15 = xhat.allparams;
clear data xhat

% load Temp increase by 5C GM15 output fields
fname_future = strcat(outputDir, 'futureproj/Tup5Cv2_GM15_CTL_He_PC.mat');
d = load(fname_future);
model_GM15_Tup5C = d.model;
clear d

% load PO4obs decrease by 20% GM15 output fields
fname_future = strcat(outputDir, 'futureproj/Pdn20pct_GM15_CTL_He_PC.mat');
d = load(fname_future);
model_GM15_Pdn20 = d.model;
clear d

% load CMIP5 projection to 2100 GM15 output fields
fname_future = strcat(outputDir, 'futureproj/CMIP2100_GM15_CTL_He_PC.mat');
d = load(fname_future);
model_GM15_CMIP = d.model;
clear d

% load Tonly model output
% load Tz model output
fname = strcat(outputDir,'optC_Tz_CTL_He_PC_DOC0.25_DOP0.mat');
load(fname);
model_Tz = data;
fxhat = strcat(outputDir,'optC_Tz_CTL_He_PC_DOC0.25_DOP0_xhat.mat');
load(fxhat);
xhat_Tz = xhat.allparams;
clear data xhat

% load Temp increase by 5C GM15 output fields
fname_future = strcat(outputDir, 'futureproj/Tup5Cv2_Tz_CTL_He_PC.mat');
d = load(fname_future);
model_Tz_Tup5C = d.model;
clear d

% load PO4obs decrease by 20% GM15 output fields
fname_future = strcat(outputDir, 'futureproj/Pdn20pct_Tz_CTL_He_PC.mat');
d = load(fname_future);
model_Tz_Pdn20 = d.model;
clear d

% load CMIP2100 Tz output fields
fname_future = strcat(outputDir, 'futureproj/CMIP2100_Tz_CTL_He_PC.mat');
d = load(fname_future);
model_Tz_CMIP = d.model;
clear d

%% -------------load grid data ---------------------

% load grid
addpath('/DFS-L/DATA/primeau/weilewang/DATA/')
addpath('/DFS-L/DATA/primeau/weilewang/my_func/'  )

load M3d91x180x24.mat MSKS
OperName = sprintf('/DFS-L/DATA/primeau/weilewang/DATA/OCIM2/OCIM2_CTL_He.mat');
    load(OperName,'output') ;
    M3d = output.M3d;
    grd = output.grid;
    TR  = output.TR/spa;
    %
iwet = find(M3d(:)) ;
dAt  = grd.DXT3d.*grd.DYT3d;
dVt  = dAt.*grd.DZT3d;
I = speye(length(iwet)) ;

%volume weight diagonal matrix
W = d0(dVt(iwet)) ;

% par.M3d = M3d;
% par.grd = grd;
% par.TRdiv = -TR     ;
% TRdiv= par.TRdiv      ;

% number of euphotic zone layers
nn = 2;
% par.nzo = nn;

%%
% Given total molecules in the atmosphere, we can recompute ppm of co2. molecules co2 to molecules atmosphere
% mass of the atmosphere = 5.148 x 10^21 g, and that air's average molar mass = 28.97 g/mol.
% moles of air in atmosphere = about 1.777e20 mol air in atmosphere
molAir = 5.148e21./28.97;

%% compute total DIC in ocean
% convert DIC units from mmolC/m^3 to mol C

% total DIC in optCell
    DICtmp = model_cell.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_Cellopt = sum(DICtmp(iwet),'all');
    fprintf('Cell model (optimal) : integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_Cellopt,tDIC_Cellopt*12*1e-15)
    % Cell model (optimal) : integrated total DIC in ocean = 3.0464e+18 mol C    (3.656e+04 Pg C)

        % all forms of carbon
        DICtmp = (model_cell.DIC+model_cell.DOC + model_cell.POC + model_cell.ALK +model_cell.PIC ).*dVt.*1e-3;
        tC_Cellopt = sum(DICtmp(iwet),'all');
    fprintf('Cell model (optimal) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n\n',tC_Cellopt,tC_Cellopt*12*1e-15)

% total DIC in Tup5C_cell
    DICtmp = model_cell_Tup5C.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_CellTup5C = sum(DICtmp(iwet),'all');
    fprintf('Cell model (Temp +5C) : integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_CellTup5C,tDIC_CellTup5C*12*1e-15)
    %Cell model (Temp +5C) : integrated total DIC in ocean = 3.0515e+18 mol C    (3.662e+04 Pg C)
    dDIC_CellTup5C = tDIC_CellTup5C - tDIC_Cellopt;
    fprintf('Cell model (Temp +5C) : change in DIC from optimal =   %+9.3e mol C   (%+6.3f Pg C) \n',dDIC_CellTup5C,dDIC_CellTup5C*12*1e-15)
    % Cell model (Temp +5C) : change in DIC from optimal =  + 5.064e+15 mol C   (+ 60.769 Pg C)

    % all forms of carbon
        DICtmp = (model_cell_Tup5C.DIC+model_cell_Tup5C.DOC + model_cell_Tup5C.POC + model_cell_Tup5C.ALK +model_cell_Tup5C.PIC ).*dVt.*1e-3;
        tC_CellTup5C = sum(DICtmp(iwet),'all');
        fprintf('Cell model (Temp +5C) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_CellTup5C,tC_CellTup5C*12*1e-15)
        dC_CellTup5C = tC_CellTup5C - tC_Cellopt;
        fprintf('Cell model (Temp +5C) : change in total C from optimal = %+9.3e mol C  (%+6.3f Pg C) \n',dC_CellTup5C,dC_CellTup5C*12*1e-15)

        dCO2ppm_CellT = -dC_CellTup5C/molAir*1e6;
        fprintf('Cell model (Temp +5C) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_CellT)
    % this would lower atmospheric CO2 by 28.93 ppm

% total DIC in Pdn20_Cell
    DICtmp = model_cell_Pdn20.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_CellPdn20 = sum(DICtmp(iwet),'all');
    fprintf('Cell model (PO4 -20%%) : integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_CellPdn20,tDIC_CellPdn20*12*1e-15)
    dDIC_CellPdn20 = tDIC_CellPdn20 - tDIC_Cellopt;
    fprintf('Cell model (PO4 -20%%) : change in DIC from optimal =    %+9.3e mol C  (%+6.3f Pg C) \n',dDIC_CellPdn20,dDIC_CellPdn20*12*1e-15)

    % all forms of carbon
        DICtmp = (model_cell_Pdn20.DIC+model_cell_Pdn20.DOC + model_cell_Pdn20.POC + model_cell_Pdn20.ALK +model_cell_Pdn20.PIC ).*dVt.*1e-3;
        tC_CellPdn20 = sum(DICtmp(iwet),'all');
        fprintf('Cell model (PO4 -20%%) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_CellPdn20,tC_CellPdn20*12*1e-15)
        dC_CellPdn20 = tC_CellPdn20 - tC_Cellopt;
        fprintf('Cell model (PO4 -20%%) : change in total C from optimal = %+9.3e mol C  (%+6.3f Pg C) \n',dC_CellPdn20,dC_CellPdn20*12*1e-15)

		dCO2ppm_CellP = -dC_CellPdn20/molAir*1e6;
        fprintf('Cell model (PO4 -20%%) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_CellP)

% total DIC in CMIP2100
    DICtmp = model_cell_CMIP.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_CellCMIP = sum(DICtmp(iwet),'all');
    fprintf('Cell model (CMIP2100) : integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_CellCMIP,tDIC_CellCMIP*12*1e-15)
    dDIC_CellCMIP = tDIC_CellCMIP - tDIC_Cellopt;
    fprintf('Cell model (CMIP2100) : change in DIC from optimal =    %+9.3e mol C  (%+6.3f Pg C) \n',dDIC_CellCMIP,dDIC_CellPdn20*12*1e-15)

    % all forms of carbon
        DICtmp = (model_cell_CMIP.DIC+model_cell_CMIP.DOC + model_cell_CMIP.POC + model_cell_CMIP.ALK +model_cell_CMIP.PIC ).*dVt.*1e-3;
        tC_CellCMIP = sum(DICtmp(iwet),'all');
        fprintf('Cell model (CMIP2100) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_CellCMIP,tC_CellCMIP*12*1e-15)
        dC_CellCMIP = tC_CellCMIP - tC_Cellopt;
        fprintf('Cell model (CMIP2100) : change in total C from optimal = %+9.3e mol C  (%+6.3f Pg C) \n',dC_CellCMIP,dC_CellCMIP*12*1e-15)

		dCO2ppm_CellC = -dC_CellCMIP/molAir*1e6;
        fprintf('Cell model (CMIP2100) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_CellC)

%-------------- GM15 model ---------------------------
% total DIC in optGM15
    DICtmp = model_GM15.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_GM15opt = sum(DICtmp(iwet),'all');
    fprintf('P-only model (optimal): integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_GM15opt,tDIC_GM15opt*12*1e-15)
    %P-only model (optimal): integrated total DIC in ocean = 3.0466e+18 mol C    (3.656e+04 Pg C)
    % all forms of carbon
        DICtmp = (model_GM15.DIC+model_GM15.DOC + model_GM15.POC + model_GM15.ALK +model_GM15.PIC ).*dVt.*1e-3;
        tC_GM15opt = sum(DICtmp(iwet),'all');
        fprintf('P-only model (optimal) : integrated total C in ocean = %10.4e mol C   (%6.3e Pg C) \n\n',tC_GM15opt,tC_GM15opt*12*1e-15)

% total DIC in Tup5C_GM15
% total DIC in Pdn20_GM15
    DICtmp = model_GM15_Pdn20.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_GM15Pdn20 = sum(DICtmp(iwet),'all');
    fprintf('P-only model (PO4 -20%%): integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_GM15Pdn20,tDIC_GM15Pdn20*12*1e-15)
    dDIC_GM15Pdn20 = tDIC_GM15Pdn20 - tDIC_GM15opt;
    fprintf('P-only model (PO4 -20%%): change in DIC from optimal =  %+10.4e mol C   (%+6.3f Pg C) \n',dDIC_GM15Pdn20,dDIC_GM15Pdn20*12*1e-15)

    % all forms of carbon
        DICtmp = (model_GM15_Pdn20.DIC+model_GM15_Pdn20.DOC + model_GM15_Pdn20.POC + model_GM15_Pdn20.ALK +model_GM15_Pdn20.PIC ).*dVt.*1e-3;
        tC_GM15Pdn20 = sum(DICtmp(iwet),'all');
        fprintf('P-only model (PO4 -20%%) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_GM15Pdn20,tC_GM15Pdn20*12*1e-15)
        dC_GM15Pdn20 = tC_GM15Pdn20 - tC_GM15opt;
        fprintf('P-only model (PO4 -20%%) : change in total C from optimal = %+9.4e mol C  (%+6.3f Pg C) \n',dC_GM15Pdn20,dC_GM15Pdn20*12*1e-15)

		dCO2ppm_GM15P = -dC_GM15Pdn20/molAir*1e6;
        fprintf('P-only model (PO4 -20%%) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_GM15P)

%%
% total DIC in Tup5C_GM15
    DICtmp = model_GM15_Tup5C.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_GM15Tup5C = sum(DICtmp(iwet),'all');
    fprintf('P-only model (Temp +5C) : integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_GM15Tup5C,tDIC_GM15Tup5C*12*1e-15)
    %GM15 model (Temp +5C) : integrated total DIC in ocean = 3.054e+18 mol C    (3.665e+04 Pg C)
    dDIC_GM15Tup5C = tDIC_GM15Tup5C - tDIC_GM15opt;
    fprintf('P-only model (Temp +5C) : change in DIC from optimal =   %+9.3e mol C   (%+6.3f Pg C) \n',dDIC_GM15Tup5C,dDIC_GM15Tup5C*12*1e-15)
    % GM15 model (Temp +5C) : change in DIC from optimal =  + 9.301e+15 mol C   (+ 111 Pg C)

    % all forms of carbon
        DICtmp = (model_GM15_Tup5C.DIC+model_GM15_Tup5C.DOC + model_GM15_Tup5C.POC + model_GM15_Tup5C.ALK +model_GM15_Tup5C.PIC ).*dVt.*1e-3;
        tC_GM15Tup5C = sum(DICtmp(iwet),'all');
        fprintf('P-only model (Temp +5C) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_GM15Tup5C,tC_GM15Tup5C*12*1e-15)
        dC_GM15Tup5C = tC_GM15Tup5C - tC_GM15opt;
        fprintf('P-only model (Temp +5C) : change in total C from optimal = %+9.3e mol C  (%+6.3f Pg C) \n',dC_GM15Tup5C,dC_GM15Tup5C*12*1e-15)

        dCO2ppm_GM15T = -dC_GM15Tup5C/molAir*1e6;
        fprintf('P-only model (Temp +5C) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_GM15T)

% Total DIC in CMIP2100 GM15
    DICtmp = model_GM15_CMIP.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_GM15CMIP = sum(DICtmp(iwet),'all');
    fprintf('P-only model (CMIP2100): integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_GM15CMIP,tDIC_GM15CMIP*12*1e-15)
    dDIC_GM15CMIP = tDIC_GM15CMIP - tDIC_GM15opt;
    fprintf('P-only model (CMIP2100): change in DIC from optimal =  %+10.4e mol C   (%+6.3f Pg C) \n',dDIC_GM15CMIP,dDIC_GM15CMIP*12*1e-15)

    % all forms of carbon
        DICtmp = (model_GM15_CMIP.DIC+model_GM15_CMIP.DOC + model_GM15_CMIP.POC + model_GM15_CMIP.ALK +model_GM15_CMIP.PIC ).*dVt.*1e-3;
        tC_GM15CMIP = sum(DICtmp(iwet),'all');
        fprintf('P-only model (CMIP2100) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_GM15CMIP,tC_GM15CMIP*12*1e-15)
        dC_GM15CMIP = tC_GM15CMIP - tC_GM15opt;
        fprintf('P-only model (CMIP2100) : change in total C from optimal = %+9.4e mol C  (%+6.3f Pg C) \n',dC_GM15CMIP,dC_GM15CMIP*12*1e-15)

		dCO2ppm_GM15C = -dC_GM15CMIP/molAir*1e6;
        fprintf('P-only model (CMIP2100) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_GM15C)

% ----------------- Tz model ---------------
% total DIC in optTz
    DICtmp = model_Tz.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_Tzopt = sum(DICtmp(iwet),'all');
    fprintf('T-only model (optimal): integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_Tzopt,tDIC_Tzopt*12*1e-15)
    %T-only model (optimal): integrated total DIC in ocean =  mol C    ( Pg C)
    % all forms of carbon
        DICtmp = (model_Tz.DIC+model_Tz.DOC + model_Tz.POC + model_Tz.ALK +model_Tz.PIC ).*dVt.*1e-3;
        tC_Tzopt = sum(DICtmp(iwet),'all');
        fprintf('T-only model (optimal) : integrated total C in ocean = %10.4e mol C   (%6.3e Pg C) \n\n',tC_Tzopt,tC_Tzopt*12*1e-15)

% total DIC in Tup5C_Tz
    DICtmp = model_Tz_Tup5C.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_TzTup5C = sum(DICtmp(iwet),'all');
    fprintf('Tz model (Temp +5C) : integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_TzTup5C,tDIC_TzTup5C*12*1e-15)
    %Tz model (Temp +5C) : integrated total DIC in ocean = 3.054e+18 mol C    (3.665e+04 Pg C)
    dDIC_TzTup5C = tDIC_TzTup5C - tDIC_Tzopt;
    fprintf('Tz model (Temp +5C) : change in DIC from optimal =   %+9.3e mol C   (%+6.3f Pg C) \n',dDIC_TzTup5C,dDIC_TzTup5C*12*1e-15)
    % Tz model (Temp +5C) : change in DIC from optimal =  + 9.301e+15 mol C   (+ 111 Pg C)

    % all forms of carbon
        DICtmp = (model_Tz_Tup5C.DIC+model_Tz_Tup5C.DOC + model_Tz_Tup5C.POC + model_Tz_Tup5C.ALK +model_Tz_Tup5C.PIC ).*dVt.*1e-3;
        tC_TzTup5C = sum(DICtmp(iwet),'all');
        fprintf('Tz model (Temp +5C) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_TzTup5C,tC_TzTup5C*12*1e-15)
        dC_TzTup5C = tC_TzTup5C - tC_Tzopt;
        fprintf('Tz model (Temp +5C) : change in total C from optimal = %+9.3e mol C  (%+6.3f Pg C) \n',dC_TzTup5C,dC_TzTup5C*12*1e-15)

        dCO2ppm_TzT = -dC_TzTup5C/molAir*1e6;
        fprintf('Tz model (Temp +5C) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_TzT)
    % this would lower atmospheric CO2 by 53.03 ppm
%%
% total DIC in Pdn20_Tz
    DICtmp = model_Tz_Pdn20.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_TzPdn20 = sum(DICtmp(iwet),'all');
    fprintf('T-only model (PO4 -20%%): integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_TzPdn20,tDIC_TzPdn20*12*1e-15)
    dDIC_TzPdn20 = tDIC_TzPdn20 - tDIC_Tzopt;
    fprintf('T-only model (PO4 -20%%): change in DIC from optimal =  %+10.4e mol C   (%+6.3f Pg C) \n',dDIC_TzPdn20,dDIC_TzPdn20*12*1e-15)

    % all forms of carbon
        DICtmp = (model_Tz_Pdn20.DIC+model_Tz_Pdn20.DOC + model_Tz_Pdn20.POC + model_Tz_Pdn20.ALK +model_Tz_Pdn20.PIC ).*dVt.*1e-3;
        tC_TzPdn20 = sum(DICtmp(iwet),'all');
        fprintf('T-only model (PO4 -20%%) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_TzPdn20,tC_TzPdn20*12*1e-15)
        dC_TzPdn20 = tC_TzPdn20 - tC_Tzopt;
        fprintf('T-only model (PO4 -20%%) : change in total C from optimal = %+9.4e mol C  (%+6.3f Pg C) \n',dC_TzPdn20,dC_TzPdn20*12*1e-15)

		dCO2ppm_TzP = -dC_TzPdn20/molAir*1e6;
        fprintf('T-only model (PO4 -20%%) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_TzP)

%%
% total DIC in CMIP_Tz
    DICtmp = model_Tz_CMIP.DIC.*dVt.*1e-3;  %units = mol C
    tDIC_TzCMIP = sum(DICtmp(iwet),'all');
    fprintf('T-only model (CMIP2100): integrated total DIC in ocean = %10.4e mol C   (%6.3e Pg C) \n',tDIC_TzCMIP,tDIC_TzCMIP*12*1e-15)
    dDIC_TzCMIP = tDIC_TzCMIP - tDIC_Tzopt;
    fprintf('T-only model (CMIP2100): change in DIC from optimal =  %+10.4e mol C   (%+6.3f Pg C) \n',dDIC_TzCMIP,dDIC_TzCMIP*12*1e-15)

    % all forms of carbon
        DICtmp = (model_Tz_CMIP.DIC+model_Tz_CMIP.DOC + model_Tz_CMIP.POC + model_Tz_CMIP.ALK +model_Tz_CMIP.PIC ).*dVt.*1e-3;
        tC_TzCMIP = sum(DICtmp(iwet),'all');
        fprintf('T-only model (CMIP2100) : integrated total C in ocean =   %10.4e mol C   (%6.3e Pg C) \n',tC_TzCMIP,tC_TzCMIP*12*1e-15)
        dC_TzCMIP = tC_TzCMIP - tC_Tzopt;
        fprintf('T-only model (CMIP2100) : change in total C from optimal = %+9.4e mol C  (%+6.3f Pg C) \n',dC_TzCMIP,dC_TzCMIP*12*1e-15)

		dCO2ppm_TzC = -dC_TzCMIP/molAir*1e6;
        fprintf('T-only model (CMIP2100) : atmospheric CO2 concentration changed by %+4.2f ppm \n\n',dCO2ppm_TzC)

%% convert to change in atmospher pCO2 (ppm)
% increase in ocean DIC = decrease in atm CO2
% Can convert moles of carbon added to atmosphere to molecules.
%dDIC_CellTup5C*6.022e23
% Then, given total molecules in the atmosphere, we can recompute ppm of co2. molecules co2 to molecules atmosphere

% mass of the atmosphere = 5.148 x 10^21 g, and that air's average molar mass = 28.97 g/mol.
% moles of air in atmosphere = about 1.777e20 mol air in atmosphere
%molAir = 5.148e21./28.97;

%dCO2ppm_CellT = -dC_CellTup5C/molAir*1e6
% this would lower atmospheric CO2 by 28.93 ppm
