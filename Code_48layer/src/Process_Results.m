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
dic14=data.c14;
H1_dic14=data.dic14h1;
H2_dic14=data.dic14h2;
Obs.dic14 = H2_dic14*dic14(:); %permil
%o2 =data.o2;
%H1_o2 = data.o2h1;
%H2_o2 = data.o2h2;
clear data


%-------------------------------------------------------
%% Step 2: Read in the results and calculate del13C and Delta14C
% Model results will be stored in the output directory

%
Path_067='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_48layer/Results/Transient_Cisotope/0.6700xkw'
Path_070='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_48layer/Results/Transient_Cisotope/0.7070xkw'
Path_080='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_48layer/Results/Transient_Cisotope/0.8000xkw'
Path_090='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_48layer/Results/Transient_Cisotope/0.9000xkw'
Path_100='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_48layer/Results/Transient_Cisotope/1.0000xkw'
%
sprintf('Doing 0.6700xkw')
[D14c_067]=load_results(Path_067,0.6700);
Model.dic14_067  = H2_dic14*H1_dic14*D14c_067(:) ;
%
sprintf('Doing 0.7070xkw')
[D14c_070]=load_results(Path_070,0.7070);
Model.dic14_070  = H2_dic14*H1_dic14*D14c_070(:) ;
%
sprintf('Doing 0.8000xkw')
[D14c_080]=load_results(Path_080,0.8000);
Model.dic14_080  = H2_dic14*H1_dic14*D14c_080(:) ;
%
sprintf('Doing 0.9000xkw')
[D14c_090]=load_results(Path_090,0.9000);
Model.dic14_090  = H2_dic14*H1_dic14*D14c_090(:) ;
%
sprintf('Doing 1.0000xkw')
[D14c_100]=load_results(Path_100,1.0000);
Model.dic14_100  = H2_dic14*H1_dic14*D14c_100(:) ;
%


%-------------------Calculate the rmse------------------
fprintf('load the transient model results and OCIM grid... \n')
addpath('../Results/Transient_Cisotope')
addpath('../../DATA/BGC_48layer/')
load OCIM2_CTL_He_48layer.mat output % call the M3d and output
M3d = output.M3d;
grd = output.grid;
dAt  = output.grid.dAt;
dVt  = output.grid.dVt;
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
dVt_ocn  = dVt(iwet) ;
dVt_docn = repmat(dVt_ocn, 1, 6*172) ; % 2month time steppinf for 172 years
Model.dVt    = H2_dic14*H1_dic14*dVt_docn(:) ;

%
vwRMSE_058 = sqrt(sum(((Obs.dic14 - Model.dic14_058).^2).*Model.dVt / sum(Model.dVt(:)))) ;
vwRMSE_067 = sqrt(sum(((Obs.dic14 - Model.dic14_067).^2).*Model.dVt / sum(Model.dVt(:)))) ;
vwRMSE_070 = sqrt(sum(((Obs.dic14 - Model.dic14_070).^2).*Model.dVt / sum(Model.dVt(:)))) ;
vwRMSE_080 = sqrt(sum(((Obs.dic14 - Model.dic14_080).^2).*Model.dVt / sum(Model.dVt(:)))) ;
vwRMSE_090 = sqrt(sum(((Obs.dic14 - Model.dic14_090).^2).*Model.dVt / sum(Model.dVt(:)))) ;
vwRMSE_100 = sqrt(sum(((Obs.dic14 - Model.dic14_100).^2).*Model.dVt / sum(Model.dVt(:)))) ;
vwRMSE_120 = sqrt(sum(((Obs.dic14 - Model.dic14_120).^2).*Model.dVt / sum(Model.dVt(:)))) ;

%% Optimize the model
x=[0.6700; 0.7070; 0.8000; 0.9000; 1.0000];
y=[vwRMSE_067; vwRMSE_070; vwRMSE_080; vwRMSE_090; vwRMSE_100];

sprintf('starting optimization')
lx=length(x);
[alpha,min_x]=solveparabola(x,y,0,lx);
x2=linspace(0.300,1.3,1000);
y2=polyval(alpha,x2);


%find min point of parabola
sprintf('Finding min point on parbola')
f=@(x) (alpha(1)*x^2)+(alpha(2)*x)+alpha(3)
initial=0;
min_x=fminsearch(f,initial)
min_y=f(min_x)
%
sprintf('Drawing fig')
figure(1)
set(gcf,'Color','w')
plot(x,y,'or','MarkerSize',10); hold on
plot(min_x,min_y,'ob','MarkerSize',10); hold on
plot(x2,y2,'-.k'); hold off
axis square
xlim([0.300 1.3]);
ylim([0 100])
grid on
xlabel('Dimensional Fitting Coefficients')
ylabel('RMSE')
legend('Root Mean Square Error','Minimum Point','Parabola')
title('Model Optimization')