clc; close all; clear all;
%% Step 1: Read in GLODAP Observations
t0=1850;
tf=2021; % up to Dec.2021

% load observational data of C13-DIC, C14-DIC, and observation from GLODAP
data = getglodap_Transient(t0,tf);
dic=data.dic;
H1_dic=data.dich1;
H2_dic=data.dich2;
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
Obs.dic   = H2_dic*dic(:); %permil
Obs.dic13 = H2_dic13*dic13(:); %permil
Obs.doc13 = H2_doc13*doc13(:); %permil
Obs.poc13 = H2_poc13*poc13(:); %permil
Obs.dic14 = H2_dic14*dic14(:); %permil
Obs.doc14 = H2_doc14*doc14(:); %permil
%Obs.o2    = H2_o2 *o2(:) ;

%-------------------------------------------------------
%% Step 2: Read in the results and calculate del13C and Delta14C
% Model results will be stored in the output directory

fprintf('load the transient model results and OCIM grid... \n')
addpath('../Results/Transient_Cisotope')
addpath('../../DATA/BGC_2023Nature/')
load OCIM2_CTL_He.mat output % call the M3d and output
M3d = output.M3d;
grd = output.grid;
iwet = find(M3d(:)) ;
nwet = length(iwet) ;
%
Path_058='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_Merge/Results/Transient_Cisotope/0.5800xkw'
Path_067='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_Merge/Results/Transient_Cisotope/0.6700xkw'
Path_070='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_Merge/Results/Transient_Cisotope/0.7070xkw'
Path_080='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_Merge/Results/Transient_Cisotope/0.8000xkw'
Path_090='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_Merge/Results/Transient_Cisotope/0.9000xkw'
Path_100='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_Merge/Results/Transient_Cisotope/1.0000xkw'
Path_120='/DFS-L/DATA/primeau/hojons1/OCIM-BGC/Code_Merge/Results/Transient_Cisotope/1.2000xkw'
%
%sprintf('Doing 0.5800xkw')
%[~,~,DIC12_058,DIC13_058,DIC14_058]=load_results(Path_058,0.5800);
%sprintf('Doing 0.6700xkw')
%[~,~,DIC12_067,DIC13_067,DIC14_067]=load_results(Path_067,0.6700);
%sprintf('Doing 0.7070xkw')
%[~,~,DIC12_070,DIC13_070,DIC14_070]=load_results(Path_070,0.7070);
sprintf('Doing 0.8000xkw')
[D14c,d13c,DI12C,DI13C,~,ddo13c,dpo13c,Ddo14c,Dpo14c]=load_results(Path_080,0.8000);
%sprintf('Doing 0.9000xkw')
%[~,~,DIC12_090,DIC13_090,DIC14_090]=load_results(Path_090,0.9000);
%sprintf('Doing 1.0000xkw')
%[~,~,DIC12_100,DIC13_100,DIC14_100]=load_results(Path_100,1.0000);
%sprintf('Doing 1.2000xkw')
%[~,~,DIC12_120,DIC13_120,DIC14_120]=load_results(Path_120,1.2000);
DICt = DI12C + DI13C ;
%
Model.dic    = H2_dic*H1_dic*DICt(:) ;
Model.dic13  = H2_dic13*H1_dic13*d13c(:) ;
Model.poc13  = H2_poc13*H1_poc13*dpo13c(:) ;
Model.doc13t = H2_doc13*H1_doc13*ddo13c(:) ;
%
Model.dic14  = H2_dic14*H1_dic14*D14c(:) ;
%Model.poc14 = H2_poc14*H1_poc14*Dpo14c(:) ;
Model.doc14t = H2_doc14*H1_doc14*Ddo14c(:) ;
%
%Model.o2    = H2_o2*H1_o2*O2(:) ;


%-------------------Calculate the rmse------------------
    
[S_dic, I_dic, rl_dic, rsqu_dic, rmse_dic] = regression(Obs.dic, Model.dic) ;
[S_dic13, I_dic13, rl_dic13, rsqu_dic13, rmse_dic13] = regression(Obs.dic13, Model.dic13) ;
[S_dic14, I_dic14, rl_dic14, rsqu_dic14, rmse_dic14] = regression(Obs.dic14, Model.dic14) ;
[S_doc13t, I_doc13t, rl_doc13t, rsqu_doc13t, rmse_doc13t] = regression(Obs.doc13, Model.doc13t) ;
[S_doc14t, I_doc14t, rl_doc14t, rsqu_doc14t, rmse_doc14t] = regression(Obs.doc14, Model.doc14t) ;
[S_poc13, I_poc13, rl_poc13, rsqu_poc13, rmse_poc13] = regression(Obs.poc13, Model.poc13) ;


%-------------Draw figure model vs obsercation--------
figure (1)
set(gcf, 'Color', 'w');
subplot(2,3,1)
scatter(Obs.dic13(:), Model.dic13(:), 'MarkerEdgeColor', 'b', 'LineWidth', 0.2 ) 
title ('\delta DI^{13}C (permil)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Observed','FontSize', 12,'FontName', 'Arial');
ylabel('Model'   ,'FontSize', 12,'FontName', 'Arial')
axis square;
xlim ([-3, 5]);
ylim ([-3, 5]);
xticks(-3:2:5);
yticks(-3:2:5);
hold on;
plot(Obs.dic13, rl_dic13, 'r-') ;
text(1, -1.7, sprintf('y = %.2fx+%.4f',S_dic13, I_dic13), 'FontSize', 12) ;
text(1, -2.3, sprintf('r^2 = %.2f, rmse = %.4f', rsqu_dic13, rmse_dic13), 'FontSize', 12) ;
%
subplot(2,3,2)
scatter(Obs.doc13(:), Model.doc13t(:), 'MarkerEdgeColor', 'b', 'LineWidth', 0.2 ) 
title ('\delta DO^{13}C (permil)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Observed','FontSize', 12,'FontName', 'Arial');
ylabel('Model'   ,'FontSize', 12,'FontName', 'Arial')
axis square;
xlim ([-30, -18]);
ylim ([-30, -18]);
xticks(-30:3:-18);
yticks(-30:3:-18);
hold on;
plot(Obs.doc13, rl_doc13t, 'r-') ;
text(-27, -26, sprintf('y = %.2fx+%.4f',S_doc13t, I_doc13t), 'FontSize', 12) ;
text(-27, -27, sprintf('r^2 = %.2f, rmse = %.4f', rsqu_doc13t, rmse_doc13t), 'FontSize', 12) ;
%
subplot(2,3,3)
scatter(Obs.poc13(:), Model.poc13(:), 'MarkerEdgeColor', 'b', 'LineWidth', 0.2 ) 
title ('\delta PO^{13}C (permil)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Observed','FontSize', 12,'FontName', 'Arial');
ylabel('Model'   ,'FontSize', 12,'FontName', 'Arial')
axis square;
xlim ([-40, -15]);
ylim ([-40, -15]);
xticks(-40:5:-15);
yticks(-40:5:-15);
hold on;
plot(Obs.poc13, rl_poc13, 'r-') ;
text(-24, -35, sprintf('y = %.2fx%.4f',S_poc13, I_poc13), 'FontSize', 12) ;
text(-24, -37, sprintf('r^2 = %.2f, rmse = %.4f', rsqu_poc13, rmse_poc13), 'FontSize', 12) ;
%
subplot(2,3,4)
scatter(Obs.dic14(:), Model.dic14(:), 'MarkerEdgeColor', 'b', 'LineWidth', 0.2 ) 
title ('\Delta DI^{14}C (permil)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Observed','FontSize', 12,'FontName', 'Arial');
ylabel('Model'   ,'FontSize', 12,'FontName', 'Arial')
axis square;
xlim ([-250, 200]);
ylim ([-250, 200]);
xticks(-250:50:200);
yticks(-250:50:200);
hold on;
plot(Obs.dic14, rl_dic14, 'r-') ;
text(-50, -180, sprintf('y = %.2fx%.4f',S_dic14, I_dic14), 'FontSize', 12) ;
text(-50, -220, sprintf('r^2 = %.2f, rmse = %.4f', rsqu_dic14, rmse_dic14), 'FontSize', 12) ;
%
subplot(2,3,5)
scatter(Obs.doc14(:), Model.doc14t(:), 'MarkerEdgeColor', 'b', 'LineWidth', 0.2 ) 
title ('\Delta DO^{14}C (permil)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Observed','FontSize', 12,'FontName', 'Arial');
ylabel('Model'   ,'FontSize', 12,'FontName', 'Arial')
axis square;
xlim ([-600, -200]);
ylim ([-600, -200]);
xticks(-600:50:-200);
yticks(-600:50:-200);
hold on;
plot(Obs.doc14, rl_doc14t, 'r-') ;
text(-480, -530, sprintf('y = %.2fx%.4f',S_doc14t, I_doc14t), 'FontSize', 12) ;
text(-480, -560, sprintf('r^2 = %.2f, rmse = %.4f', rsqu_doc14t, rmse_doc14t), 'FontSize', 12) ;
%
subplot(2,3,6)
scatter(Obs.dic(:), Model.dic(:), 'MarkerEdgeColor', 'b', 'LineWidth', 0.2 ) 
title ('DIC (\mu mol kg^{-1})', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Observed','FontSize', 12,'FontName', 'Arial');
ylabel('Model'   ,'FontSize', 12,'FontName', 'Arial')
axis square;
xlim ([1500, 2500]);
ylim ([1500, 2500]);
xticks(1500:500:2500);
yticks(1500:500:2500);
hold on;
plot(Obs.dic, rl_dic, 'r-') ;
text(2000, 1650, sprintf('y = %.2fx+%.4f',S_dic, I_dic), 'FontSize', 10) ;
text(2000, 1600, sprintf('r^2 = %.2f, rmse = %.4f', rsqu_dic, rmse_dic), 'FontSize', 10) ;

%--------------Sectional distribution--------------
delDIC13  = M3d + nan;  delDIC13(iwet)  = deltaDIC13(:,end) ;
DelDIC14  = M3d + nan;  DelDIC14(iwet)  = DeltaDIC14(:,end) ;
delDOC13t = M3d + nan;  delDOC13t(iwet) = deltaDOC13t(:,end) ;
DelDOC14t  = M3d + nan; DelDOC14t(iwet) = DeltaDOC14t(:,end) ;


figure (2)
subplot (2,2,1);
contourf(grd.yt,-grd.zt,squeeze(delDIC13(:,165,:))')
title ('\deltaDI^1^3C (Atlantic: 165^oE)','FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Latitude(^oN)','FontSize', 14,'FontName', 'Arial');
ylabel('Depth (m)','FontSize', 14,'FontName', 'Arial')
colorbar
caxis([-2, 4])


subplot (2,2,2);
contourf(grd.yt,-grd.zt,squeeze(DelDIC14(:,165,:))')
title ('\DeltaDI^1^4C (Atlantic: 165^oE)','FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Latitude(^oN)','FontSize', 14,'FontName', 'Arial');
ylabel('Depth (m)','FontSize', 14,'FontName', 'Arial')
colorbar
caxis([-150, 50])

subplot (2,2,3);
%contourf(grd.yt,-grd.zt,squeeze(delDOC13t(:,165,:))')
%title ('\deltaDO^1^3C (Atlantic: 165^oE)')
%xlabel('Latitude(^oN)');
%ylabel('Depth (m)')
%colorbar
pcolor(M3d(:,:,1))

subplot (2,2,4);
contourf(grd.yt,-grd.zt,squeeze(DelDOC14t(:,165,:))')
title ('\DeltaDO^1^4C (Atlantic: 165^oE)','FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial')
xlabel('Latitude(^oN)','FontSize', 14,'FontName', 'Arial');
ylabel('Depth (m)')
colorbar
caxis([-400, -200])

%--------------Vertical Structure------------------

mean_delDIC13 = squeeze(nanmean(delDIC13, [1, 2])) ;
std_delDIC13  = squeeze(nanstd(delDIC13, 0, [1, 2])) ;

mean_DelDIC14 = squeeze(nanmean(DelDIC14, [1, 2])) ;
std_DelDIC14  = squeeze(nanstd(DelDIC14, 0, [1, 2])) ;

mean_delDOC13t = squeeze(nanmean(delDOC13t, [1, 2])) ;
std_delDOC13t  = squeeze(nanstd(delDOC13t, 0, [1, 2])) ;

mean_DelDOC14t = squeeze(nanmean(DelDOC14t, [1, 2])) ;
std_DelDOC14t  = squeeze(nanstd(DelDOC14t, 0, [1, 2])) ;


%
figure (3) ;
subplot (1,2,1) ;
scatter(mean_delDIC13(:), -grd.zt(:), 100, 's', 'filled') ;
axis square;
hold on;
%errorbar(mean_delDIC13, -grd.zt, std_delDIC13, 'horizontal', 'k', 'LineStyle','none','LineWidth',0.5) ;
xlabel('\deltaDI^1^3C (Model global mean)','FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial') ;
ylabel('Depth (m)','FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial') ;
ax = gca ;
ax.XAxisLocation = 'top' ;
xlim ([0.5, 1.5]);
xticks([0.5:0.25:1.5]);
hold off ;

subplot (1,2,2) ;
scatter(mean_DelDIC14(:), -grd.zt(:), 100, 's', 'filled') ;
axis square;
hold on;
errorbar(mean_DelDIC14, -grd.zt, std_DelDIC14, 'horizontal', 'k', 'LineStyle','none','LineWidth',0.5) ;
xlabel('\Delta^{14}C (Model global mean)', 'FontWeight', 'bold') ;
ylabel('Depth (m)','FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial') ;
ax = gca ;
ax.XAxisLocation = 'top' ;
xlim ([-500, 100]);
hold on ;
scatter(mean_DelDOC14t(:), -grd.zt(:), 100, 's', 'filled') ;
axis square;
hold on;
errorbar(mean_DelDOC14t, -grd.zt, std_DelDOC14t, 'horizontal', 'k', 'LineStyle','none','LineWidth',0.5) ;
xlabel('\Delta^{14}C (Model global mean)','FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial') ;
ylabel('Depth (m)','FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial') 
ax = gca ;
ax.XAxisLocation = 'top' ;
hold off ;

%-----------Save the figures------------
fname = 'Results_1.000xkw.tiff'
%
figure(1);
frame = getframe(gcf);
imwrite(frame.cdata, fname, 'WriteMode', 'overwrite', 'Compression', 'none');

% Figure 2 저장
figure(2);
frame = getframe(gcf);
imwrite(frame.cdata, fname, 'WriteMode', 'append', 'Compression', 'none');

% Figure 3 저장
figure(3);
frame = getframe(gcf);
imwrite(frame.cdata, fname, 'WriteMode', 'append', 'Compression', 'none');

disp(['Figures saved as multi-page TIFF in ', fname]);
