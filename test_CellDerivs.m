% CellCNP derivative test
%% load inputs
load('/Users/megansullivan/Documents/UC Irvine/DATASETS/TraitModel_TestData/inputs_surf_CellCNP.mat')
%load('/Users/megansullivan/Documents/UC Irvine/DATASETS/TraitModel_TestData/parBIO.mat')
%par.BIO  =  BIO;
%par.BIO.r0Cutoff = 2.25;
%par = param;
par.BIO = parBIO;
%par.BIO.lPCutoff = -7.252;

on = true; off = false;
par.opt_Q10Photo = off;
par.opt_fRibE 	 = off;
par.opt_fStorage = off;
par.opt_kST0 	 = off;
par.opt_PStor_rCutoff = off;
par.opt_PStor_scale = off;
par.opt_PLip_PCutoff = off;
par.opt_PLip_scale = off;
par.opt_alphaS = off;


par.pindx = pindx;
par.x0=x0;

P0(P0<0)= real(min(P0(P0>=0)));

%% testing with model P output
load('/Users/megansullivan/Documents/UC Irvine/GitHub/TraitModel_output/FIGS_02May21/Int_CNPP.mat')

j=1:length(N0);
DIP = DIPsurf(iprod);
[CellOut, parBIO] = CellCNP(par,x0, DIP(j),N0(j),T0(j),Irr0(j))

mu_model = M3dsurf*0;
mu_model(iprod) = CellOut.mu;

%save('mu_model.mat','mu_model');
%% set figure properties
%%%% only run once per matlab session

set(groot,'defaultAxesFontName','Times',...
    'defaultAxesFontSize',14,...
    'defaultAxesTickLabelInterpreter','latex',...   
    'defaultAxesXMinorTick','on',...
    'defaultAxesYMinorTick','on');

% TEXT PROPERTIES
set(groot,'defaultTextFontName','Times')%,...
    %'defaultTextInterpreter','tex');

%%
j=1:length(P0);
%j=ibad(1:30);
[CellOut, parBIO] = CellCNP(par,x0, P0(j),N0(j),T0(j),Irr0(j))


%%
		%par.BIO = parBIO;
		%clear parBIO;
		par.CellOut.C2P = M3dsurf*0;
		par.CellOut.N2P = M3dsurf*0;
		par.CellOut.C2N = M3dsurf*0;
		par.CellOut.LimType = M3dsurf*NaN;
		par.CellOut.r = M3dsurf*NaN;

		par.CellOut.C2P(iprod) = CellOut.CP;
        par.CellOut.N2P(iprod) = CellOut.NP;
		par.CellOut.C2N(iprod) = CellOut.CN;
        
        par.CellOut.LimType(iprod) = CellOut.LimType;
		par.CellOut.r(iprod) = CellOut.r;

		%par.CellOut.C2P(isnan(par.CellOut.C2P)) = 0; %remove NaNs



%% kST0
param_name = 'kST0';
lparam_name = 'lkST0';
dparam_name = 'dC2P_dkST0';
opt_name = 'opt_kST0'

par.opt_kST0 	 = on;
par.BIO.kST0 = 0.185;

x0=[];
        logparam = log(par.BIO.(param_name));
		strt = length(x0) + 1;
		x0 = [x0; logparam];
		par.pindx.(lparam_name) = strt : length(x0);
        xim = zeros(size(x0));
        xim(par.pindx.(lparam_name)) = sqrt(-1)*eps^3;
        
j=1:length(P0);
[CellOut, parBIO] = CellCNP(par,x0+xim, P0(j),N0(j),T0(j),Irr0(j))
dC2P_dlkST0 =real(CellOut.dC2P_dkST0)*par.BIO.kST0;
dC2P_dlKST0_CSD=imag(CellOut.CP)./eps.^3;

%diff = (dC2P_dkST0 - dC2P_dKST0_CSD)./(dC2P_dKST0_CSD);
ekST0 = (dC2P_dlkST0 - dC2P_dlKST0_CSD);

ibad = find(abs(ekST0)>1);
maxk(abs(ekST0),20)

%% fRibE
param_name = 'fRibE';
lparam_name = 'tfRibE';
dparam_name = 'dC2P_dfRibE';
opt_name = 'opt_fRibE'

par.opt_fRibE 	 = on;
par.BIO.fRibE = .618;
x0=[];
        tfRibE = atanh(2*par.BIO.fRibE-1);
		strt = length(x0) + 1;
		x0 = [x0; tfRibE];
		par.pindx.(lparam_name) = strt : length(x0);
        xim = zeros(size(x0));
        xim(par.pindx.(lparam_name)) = sqrt(-1)*eps^3;
j=1:length(P0);
[CellOut, parBIO] = CellCNP(par,x0+xim, P0(j),N0(j),T0(j),Irr0(j))
dC2P_dlfRibE =real(CellOut.dC2P_dfRibE)*0.5*sech(tfRibE)^2; 	% dfRibE/dtfRibE = 0.5*sech(tfRibE)^2
dC2P_dlfRibE_CSD=imag(CellOut.CP)./eps.^3;

e = dC2P_dlfRibE - dC2P_dlfRibE_CSD;
%%
mink(abs(e),10)
ibad2 = find(abs(e)>0.1);  %find where derivative is wrong

unique(CellOut.LimType(ibad2)); %fRibE is wrong in co-limited case only
        
%% check C2P derivatives w.r.t BIO parameters
ver = datestr(now,'mmmdd');
outPath='/Users/megansullivan/Documents/UC Irvine/GitHub/OCIM-BGC-Cell/FIGS/'
%turn everything off
on = true; off = false;
par.opt_Q10Photo = off;
par.opt_fRibE 	 = off;
par.opt_fStorage = off;
par.opt_kST0 	 = off;
par.opt_PStor_rCutoff = off;
par.opt_PStor_scale = off;
par.opt_PLip_PCutoff = off;
par.opt_PLip_scale = off;
par.opt_alphaS = off;



%% choose parameter
param_name = 'alphaS';
lparam_name = 'lalphaS';
dparam_name = 'dC2P_dalphaS';
opt_name = 'opt_alphaS'
step = 0.01;
param_range = [0.01:step:4];

% choose a point
j = 100;

%% dC2P_dkST0
param_name = 'kST0';
lparam_name = 'lkST0';
dparam_name = 'dC2P_dkST0';
opt_name = 'opt_kST0'
step = 0.01;
param_range = [0.01:step:1];

% choose a point
j = 100;  %looks good
%j = 5100; %looks good (limType 2,0)
j = 308;  % limType 1: wrong
j = 2534;

%% dC2P_dfRibE %%%% need to modify to match tanh version of fRibE
param_name = 'fRibE';
lparam_name = 'tfRibE';
dparam_name = 'dC2P_dfRibE';
opt_name = 'opt_fRibE'
step = 0.01;
param_range = [0.01:step:1]; %0.618 

% choose a point
j = 5100;

%% turn on chosen parameter opt
par.(opt_name) = on;

XX = NaN([length(param_range),1]); 
C2P  = XX;
dC2P = XX;
dC2P_CSD = XX;
ddC2P_CSD = XX;
LimType = XX;

for i=1:length(param_range)
    par.BIO.(param_name) = param_range(i);
    x0=[];
        logparam = log(par.BIO.(param_name));
		strt = length(x0) + 1;
		x0 = [x0; logparam];
		par.pindx.(lparam_name) = strt : length(x0);
        xim = zeros(size(x0));
        xim(par.pindx.(lparam_name)) = sqrt(-1)*eps^3;
        
    [out, ~] = CellCNP(par,x0+xim, P0(j),N0(j),T0(j),Irr0(j));
    C2P(i)=real(out.CP(1));
    dC2P(i) =real(out.(dparam_name));
    dC2P_CSD(i)=imag(out.CP(1))/eps^3;
    ddC2P_CSD(i) = imag(out.(dparam_name))/eps^3;
    LimType(i)=out.LimType;
end

[~,limorder]=unique(LimType,'first');
L = LimType(sort(limorder))';

% finite difference
dC2P_diff=diff(C2P)./step;
%diffjump = maxk(dC2P_diff,length(L));
%dC2P_diff(dC2P_diff>300)=NaN;
%dC2P_diff(dC2P_diff<-300)=NaN;

% deriv wrt log(param)
dC2P_dlogparam = dC2P.*param_range';


% picking where to display annotations
len = length(param_range);

%% plot
figure; hold on
yyaxis left
plot(param_range,C2P,'-b')
    tind = floor(len/6);
    text(param_range(tind),C2P(tind),'\leftarrow C:P+stor','Color','b')
ylabel('C:P')

yyaxis right
plot(param_range,dC2P,'-r','LineWidth',1); hold on
    tind = floor(len/4);
    text(param_range(tind),dC2P(tind),'\leftarrow explicit derivative','Color','red')  
plot(param_range(2:end)-0.5*step,dC2P_diff,'--m','LineWidth',1.5)
    tind = floor(len/3);
    text(param_range(tind),dC2P_diff(tind),'\leftarrow finite diff','Color','m')
plot(param_range,dC2P_CSD,'-g')
    tind = floor(len/2);
    text(param_range(tind),dC2P_CSD(tind),['\leftarrow Complex step deriv of ', lparam_name],'Color',[0 0.5 0])       
plot(param_range,dC2P_dlogparam,'-.r','LineWidth',1.5)
    tind = floor(len*2/3);
    text(param_range(tind),dC2P_dlogparam(tind),['\leftarrow explicit deriv of ',lparam_name],'Color','red')       
ylabel(['dC:P / d', param_name])
xlabel(param_name)
grid on

title(['dC2P/d' param_name ' Test: indx = ' num2str(j) '; LimType = ' num2str(L)])   

figTitle = [dparam_name '_test_ind' num2str(j)];
%print(gcf,[outPath 'FIG_' figTitle '_' ver '.png'],'-dpng')
