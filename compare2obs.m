clc; clear all; close all
addpath('/DFS-L/DATA/primeau/weilewang/DATA/');
addpath('/DFS-L/DATA/primeau/weilewang/my_func');
addpath('/DFS-L/DATA/primeau/weilewang/GREG/Couple_CP/')
load transport_v4.mat
load po4obs_90x180x24.mat
load GLODAP_grid_dic
load human_co2
load theta_from_qqq.mat theta

load o2obs_90x180x24 o2obs 
load splco2_mod_monthly % monthly CO2 data
load sio4obs_90x180x24.mat sio4obs
load Sobs_90x180x24.mat
load tempobs_90x180x24.mat
load SOxhalf_P.mat
load SOxhalf_C.mat
load SOxhalf_O2.mat

grd = grid;
iwet = find(M3d(:));
nwet = length(iwet);

%%%%%%%%%%%%%%%%% compare DIP  %%%%%%%%%%%%%%
figure(1)
DIPobs = po4obs; 
ipo4 = find(DIP(iwet)>0);
W = (dVt(iwet(ipo4))./sum(dVt(iwet(ipo4))));

cr = 5:5:95;
data = [DIPobs(iwet(ipo4)),DIP(iwet(ipo4))];
O = DIPobs(iwet(ipo4));
M = DIP(iwet(ipo4));
rsquare(O,M)
[bandwidth,density,X,Y] = mykde2d(data,200,[0 0],[4 4],W);

dx = X(3,5)-X(3,4); 
dy = Y(4,2)-Y(3,2);
[q,ii] = sort(density(:)*dx*dy,'descend');
D = density;
D(ii) = cumsum(q);
subplot('position',[0.2 0.2 0.6 0.6])
contourf(X,Y,100*(1-D),cr); hold on
contour(X,Y,100*(1-D),cr);

caxis([5 95])
%set(gca,'FontSize',16);
grid on
axis square
xlabel('Observed DIP (mmol/m^3)');
ylabel('Model DIP (mmol/m^3)');
% title('model V.S. observation')
plot([0 4],[0 4],'r--','linewidth',2);

subplot('position',[0.82 0.2 0.05 0.6]);
contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
contour([1 2],cr,[cr(:),cr(:)],cr);
hold off
%set(gca,'FontSize',14);
set(gca,'XTickLabel',[]);
set(gca,'YAxisLocation','right');
set(gca,'TickLength',[0 0])
ylabel('(percentile)')

exportfig(gcf,'DIP_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
%% ---------------------------------------------------
figure(2)
o2obs = o2obs;
% convert unit form [ml/l] to [umol/l].
o2obs = o2obs.*44.661;  
% o2 correction based on Bianchi et al.(2012) [umol/l] .
o2obs = o2obs.*1.009-2.523;
io2 = find(O2(iwet)>0);
W = (dVt(iwet(io2))./sum(dVt(iwet(io2))));

cr = 5:5:95;
M = O2(iwet(io2));
O = o2obs(iwet(io2));
data = [O, M];

rsquare(O,M)
[bandwidth,density,X,Y] = mykde2d(data,200,[0 0],[300 300],W);

dx = X(3,5)-X(3,4); 
dy = Y(4,2)-Y(3,2);
[q,ii] = sort(density(:)*dx*dy,'descend');
D = density;
D(ii) = cumsum(q);
subplot('position',[0.2 0.2 0.6 0.6])
contourf(X,Y,100*(1-D),cr); hold on
contour(X,Y,100*(1-D),cr);

caxis([5 95])
%set(gca,'FontSize',16);
grid on
axis square
xlabel('Observed O2 (mmol/m^3)');
ylabel('Model O2 (mmol/m^3)');
% title('model V.S. observation')
plot([0 300],[0 300],'r--','linewidth',2);

subplot('position',[0.82 0.2 0.05 0.6]);
contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
contour([1 2],cr,[cr(:),cr(:)],cr);
hold off
%set(gca,'FontSize',14);
set(gca,'XTickLabel',[]);
set(gca,'YAxisLocation','right');
set(gca,'TickLength',[0 0])
ylabel('(percentile)')
exportfig(gcf,'O2_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
%%-----------------------------------------------------
figure(3)
rho = 1024.5     ; % seawater density;
permil = rho*1e-3; % from umol/kg to mmol/m3;

DICobs = dic*permil; % GLODAP dic obs [mmol/m3];
human_co2 = human_co2*permil;

DICobs = DICobs;
iDIC = find(DIC(iwet)>0);
W = (dVt(iwet(iDIC))./sum(dVt(iwet(iDIC))));

cr = 5:5:95;
M = DIC(iwet(iDIC)) + human_co2(iwet(iDIC));
O = DICobs(iwet(iDIC));
data = [O, M];
rsquare(O,M)
[bandwidth,density,X,Y] = mykde2d(data,200,[2000 2000],[2500 2500],W);

dx = X(3,5)-X(3,4); 
dy = Y(4,2)-Y(3,2);
[q,ii] = sort(density(:)*dx*dy,'descend');
D = density;
D(ii) = cumsum(q);
subplot('position',[0.2 0.2 0.6 0.6])
contourf(X,Y,100*(1-D),cr); hold on
contour(X,Y,100*(1-D),cr);

caxis([5 95])
%set(gca,'FontSize',16);
grid on
axis square
xlabel('Observed DIC (mmol/m^3)');
ylabel('Model DIC (mmol/m^3)');
% title('model V.S. observation')
plot([2000 2500],[2000 2500],'r--','linewidth',2);

subplot('position',[0.82 0.2 0.05 0.6]);
contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
contour([1 2],cr,[cr(:),cr(:)],cr);
hold off
%set(gca,'FontSize',14);
set(gca,'XTickLabel',[]);
set(gca,'YAxisLocation','right');
set(gca,'TickLength',[0 0])
ylabel('(percentile)')

exportfig(gcf,'DIC_MvsO','fontmode','fixed','fontsize',12,'color','rgb','renderer','painters')
