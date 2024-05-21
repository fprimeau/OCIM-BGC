%% Process_Results.m 
%Author: Shayma Al Ali 3/14/2024. Chaned in 04/30/2024 by Hojong.

%Use script to process the results from the OCIM model, compare with
%GLODAP observations, calculate the Root Mean Square Error, 
% and optimize the model

%Step 1: Read in GLODAP Observations

%Step 2: Read in the results and calculate D14C

%Step 3: Calculate the RMSE 

%Step 4: Optimize the model
clc; close all; clear all;
%% Step 1: Read in GLODAP Observations
t0=1850;
tf=2021; % up to Dec.2021

% load observational data of C13-DIC, C14-DIC, and observation from GLODAP
data = getglodap_Transient(t0,tf);
dic13=data.c13;
H1_dic13=data.dic13h1;
H2_dic13=data.dic13h2;
dic14=data.c14;
H1_dic14=data.dic14h1;
H2_dic14=data.dic14h2;
%o2 =data.o2;
%H1_o2 = data.o2h1;
%H2_o2 = data.o2h2;
clear data

% load observational data of C13-DOC and C14-DOC from data merged files (Hansel et al. (2021) and Verwega et al. (2021))
data = getOC_Transient(t0,tf);
doc13=data.doc13;
H1_doc13=data.doc13h1;
H2_doc13=data.doc13h2;
doc14=data.doc14;
H1_doc14=data.doc14h1;
H2_doc14=data.doc14h2;
poc13 =data.poc13;
H1_poc13 = data.poc13h1;
H2_poc13 = data.poc13h2;
clear data

%find the observations that are in the OCIM grid
Obs.dic13 = H2_dic13*dic13(:); %permil
Obs.doc13 = H2_doc13*doc13(:); %permil
Obs.poc13 = H2_poc13*poc13(:); %permil
Obs.dic14 = H2_dic14*dic14(:); %permil
Obs.doc14 = H2_doc14*doc14(:); %permil
%Obs.o2    = H2_o2 *o2(:) ;



%-------------------------------------------------------
%% Step 2: Read in the results and calculate del13C and Delta14C
%Model results will be stored in the output directory

fprintf('load the transient model results and OCIM grid')
addpath('../Results/Transient_Cisotope')
addpath('../../DATA/BGC_2023Nature/')
load Transient_0.7777xkw_fras=1.00_frpho=0.50_fc14=2.00.mat
load OCIM2_CTL_He.mat output % call the M3d and output
M3d = output.M3d;
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
 
% Set the calculating delta13c and Deltah14c
R14oxa=1.176*1e-12;
R13pdb=1.12372*1e-2;
R13 = @(c13,c12) c13./(c13+c12);
R14 = @(c14,c13,c12) c14./( c13 + c12 );
delta13c = @(c13,c12) ( R13(c13,c12) ./ ( (1 - R13(c13,c12) ) * R13pdb ) - 1 ) * 1e3;
Delta14c = @(c14,c13,c12) (( R14(c14,c13,c12) ./ R14oxa ) .* ( 0.975 ./ ( 1 + delta13c(c13,c12) / 1e3 ) ).^2 -1 ) * 1e3;
    
%
DI12C=Xout.C12(0*nwet+1:1*nwet,:);
PO12C=Xout.C12(1*nwet+1:2*nwet,:);
DO12C=Xout.C12(2*nwet+1:3*nwet,:);
%
DI13C=Xout.C13(0*nwet+1:1*nwet,:);
PO13C=Xout.C13(1*nwet+1:2*nwet,:);
DO13C=Xout.C13(2*nwet+1:3*nwet,:);
%
DI14C=Xout.C14(0*nwet+1:1*nwet,:);
PO14C=Xout.C14(1*nwet+1:2*nwet,:);
DO14C=Xout.C14(2*nwet+1:3*nwet,:);
%
%O2   =Xout.O2(1:nwet,:) ;
    
%-------Calculate the delta values for C isotopes-----------
deltaDIC13=delta13c(DI13C,DI12C) ;
deltaPOC13=delta13c(PO13C,PO12C) ;
deltaDOC13=delta13c(DO13C,DO12C) ;
%
DeltaDIC14=Delta14c(DI14C,DI13C,DI12C) ;
DeltaPOC14=Delta14c(PO14C,PO13C,PO12C) ;
DeltaDOC14=Delta14c(DO14C,DO13C,DO12C) ;

% final column means Jan.2022. BUt GLODAP only has data up to Dec.2021
deltaDIC13=deltaDIC13(:,1:end-1)  ; 
deltaPOC13=deltaPOC13(:,1:end-1)  ; 
deltaDOC13=deltaDOC13(:,1:end-1)  ; 
DeltaDIC14=DeltaDIC14(:,1:end-1)  ;
DeltaPOC14=DeltaPOC14(:,1:end-1)  ;
DeltaDOC14=DeltaDOC14(:,1:end-1)  ;
%O2=O2(:,1:end-1) ;

Model.dic13 = H2_dic13*H1_dic13*deltaDIC13(:) ;
Model.poc13 = H2_poc13*H1_poc13*deltaPOC13(:) ;
Model.doc13 = H2_doc13*H1_doc13*deltaDOC13(:) ;
Model.dic14 = H2_dic14*H1_dic14*DeltaDIC14(:) ;
% Model.poc14 = H2_poc14*H1_poc14*DeltaPOC14(:) ;   % If we can compile the poc14 data, it would work
Model.doc14 = H2_doc14*H1_doc14*DeltaDOC14(:) ;
%Model.o2    = H2_o2*H1_o2*O2(:) ;
    
%save the model and obs data
fileName  = 'ModelvsObs_0.7777xkw_fras=1.00_frpho=0.50_fc=2.00.mat'
directory = '../Results/Transient_Cisotope'
filePath  = fullfile(directory, fileName) ;
save(filePath,'Obs','Model','-v7.3');

%% Step 3: Calculate the RMSE
%sprintf('Starting RMSE calc')
%RMSE_c13=rmse(Model_c13,Obs_c13);
%RMSE_c14=rmse(Model_c14,Obs_c14);
%RMSE_o2 =rmse(Model_o2,Obs_o2);

%% Optimize the model
%x=[0.20;0.251;0.30;0.32];
%y=[RMSE_020;RMSE_025;RMSE_030;RMSE_032];
%sprintf('starting optimization')
%[alpha,min_x]=solveparabola(x,y,0);
%x2=linspace(0.200,0.4,1000);
%y2=polyval(alpha,x2);

%find min point of parabola
%sprintf('Finding min point on parbola')
%f=@(x) (alpha(1)*x^2)+(alpha(2)*x)+alpha(3)
%initial=0;
%min_x=fminsearch(f,initial)
%min_y=f(min_x)
    
%sprintf('Drawing fig')
%figure(1)
%plot(x,y,'or','MarkerSize',10); hold on
%plot(min_x,min_y,'ob','MarkerSize',10); hold on
%plot(x2,y2,'-.k'); hold off
%xlim([0.200 0.4]);
%ylim([50 250])
%grid on
%xlabel('Dimensional Fitting Coefficients')
%ylabel('RMSE')
%legend('Root Mean Square Error','Minimum Point','Parabola')
%title('Model Optimization')
%saveas(gcf,'Wan_opt_04222024.fig') 

%sprintf('All done! :)')