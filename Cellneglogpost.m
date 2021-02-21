%% Cell model only neglogpost optim

% save more fields from par
% out.iwet = par.iwet;
% out.nwet = par.nwet;
% out.dVt  = par.dVt   ;
% out.M3d  = par.M3d   ;
% out.po4raw = par.po4raw;
% out.no3raw = par.no3raw;
% out.no3obs = par.no3obs;
% out.Tobs = par.Tobs;
% out.PARobs = par.PARobs;
% out.BIO = par.BIO;


on = true; off = false;
%% load inputs
addpath('/Users/megansullivan/Documents/UC Irvine/DATASETS/TraitModel_TestData/')
load('/Users/megansullivan/Documents/UC Irvine/DATASETS/TraitModel_TestData/inputs_surf_CellCNP.mat')
%load('/Users/megansullivan/Documents/UC Irvine/DATASETS/TraitModel_TestData/parBIO.mat')

on = true; off = false;

par.BIO = parBIO;
par.pindx = pindx;
par.x0=x0;
par.optim = on;

%%% setup
par.P0 = P0; par.N0=N0; par.T0=T0; par.Irr0=Irr0; 
par.M3d = M3dsurf; 
par.dVt = dVtsurf;
par.iprod = iprod;
par.iwet = iprod;
par.nwet = length(iprod);


% if loading in Greenplanet model
%    M3dsurf=par.M3dsurf;
%    iwet = par.iwet;
%    nwet = par.nwet;
%    
%    iprod = find(M3d(:,:,1:2)); %production in top two layers
%    par.iprod = iprod;
%    
%	P0 = par.DIPgrd(iprod)./10^6; 			% convert ug/m^3 to g/m^3
%	N0 = par.no3obs(iprod)./10^6;   % convert ug/m^3 to g/m^3
%	T0 = par.Tobs(iprod);
%	Irr0 = par.PARobs(iprod);
dVt  = par.dVt   ;

%%
par.opt_Q10Photo = off;
par.opt_fStorage = off;
par.opt_fRibE 	 = on;
par.opt_kST0 	 = off;
par.opt_PStor_rCutoff = off;
par.opt_PStor_scale = off;
par.opt_PLip_PCutoff = off;
par.opt_PLip_scale = off;
par.opt_alphaS = off;
%% create a map from Teng2014 C2P

Teng = load('teng_regions.mat')
%load('teng_region_91x180.mat')
Reg = Teng.R(:,:,1:2);
CPTeng = NaN(size(Reg));
i1 = find(Reg==1);
i2 = find(Reg==2);
i3 = find(Reg==3);
i4 = find(Reg==4);
i5 = find(Reg==5);
i6 = find(Reg==6);
i7 = find(Reg==7);
i8 = find(Reg==8);
i9 = find(Reg==9);
i10 = find(Reg==10);
i11 = find(Reg==11);
i12 = find(Reg==12);

CPTeng(i1) = 355; %NAtl gyre
CPTeng(i2) = 81; %EqAtl
CPTeng(i3) = 163; %SAtl
CPTeng(i4) = 91; %southern ocean
CPTeng(i5) = 115; % madagascar-aus
CPTeng(i6) = 103; %indian
CPTeng(i7) = 138; %aus-pacific
CPTeng(i8) = 83; %eq. pacific
CPTeng(i9) = 176; %N-central pacific
CPTeng(i10) = 86; %N pacific
CPTeng(i11) = 80;    %acrctic  %no value here. arbitrarily added
CPTeng(i12) = 63; % N atl

par.CPTeng = CPTeng;


% load M3d91x180x24.mat MSKS
% ARC  = MSKS.ARC     ;
% MED  = MSKS.MED     ;
% PAC  = MSKS.PAC     ;
% ATL  = MSKS.ATL     ;
% IND  = MSKS.IND     ;
% iarc = find(ARC(:)) ;
% imed = find(MED(:)) ;
% ipac = find(PAC(:)) ;
% iatl = find(ATL(:)) ;
% iind = find(IND(:)) ;


%% weighting function
iteng = find(CPTeng(iwet)>0) ;
    Wb   = d0(dVt(iwet(iteng))/sum(dVt(iwet(iteng)))) ;
    mu   = sum(Wb*CPTeng(iwet(iteng)))/sum(diag(Wb)) ;
    var  = sum(Wb*(CPTeng(iwet(iteng))-mu).^2)/sum(diag(Wb));
    Wb   = Wb/var ;
    
    
    Wb   = d0(dVt(iprod)/sum(dVt(iprod))) ;


%%
%[par2, C2P] = C2Puptake(x0,par);

x0=[];
kk=0;
if par.opt_Q10Photo
    kk = kk + 1;
    par.pindx.lQ10Photo = kk;
    x0(kk) = log(par.BIO.Q10Photo);
end
if par.opt_fStorage
    kk = kk + 1;
    par.pindx.lfStorage = kk;
    x0(kk) = log(par.BIO.fStorage);
end
if par.opt_fRibE
    kk = kk + 1;
    par.pindx.lfRibE = kk;
    x0(kk) = log(par.BIO.fRibE);
end
if par.opt_PLip_PCutoff
    kk = kk + 1;
    par.pindx.lPLip_PCutoff = kk;
    x0(kk) = log(par.BIO.PLip_PCutoff);
end
if par.opt_PLip_scale
    kk = kk + 1;
    par.pindx.lPLip_scale = kk;
    x0(kk) = log(par.BIO.PLip_scale);
end
if par.opt_PStor_rCutoff
    kk = kk + 1;
    par.pindx.lPStor_rCutoff = kk;
    x0(kk) = log(par.BIO.PStor_rCutoff);
end
if par.opt_PStor_scale
    kk = kk + 1;
    par.pindx.lPStor_scale = kk;
    x0(kk) = log(par.BIO.PStor_scale);
end
if par.opt_alphaS
    kk = kk + 1;
    par.pindx.lalphaS = kk;
    x0(kk) = log(par.BIO.alphaS);
end

nbx = kk;
par.nbx = nbx;

%% Optimize cell model parameters
global iter
iter = 0 ;
par_old = par;      
par.ftest='testx0.mat';
myfun = @(x) C2Pneglogpost(x, par,iter);
%
options = optimoptions(@fminunc                  , ...   
                       'Display','iter'          , ...
                       'DerivativeCheck','off'   , ... 
                       'Algorithm','trust-region', ...
                       'GradObj','on'            , ...
                       'Hessian','on'            , ...
                       'MaxFunEvals',2000        , ...
                       'MaxIter',2000            , ...
                       'TolX',5e-7               , ...
                       'TolFun',5e-7             , ...
                       'PrecondBandWidth',Inf);
                   
                   % 'FinDiffType','central'   , ...
                    

[xhat,fval,exitflag] = fminunc(myfun,x0,options);
%[f,fx,fxx] = neglogpost(xhat,par,par_old);

% xhat: logQ10 = -1.8577    logfStorage 5.2378
%          Q10   0.1560        fStorage 188.2589

%lQ10photo =0.0133, lfstorage = 0.0298, PCutoff=2.3196e-07, rCutoff =
%121.8822, alphaS =3.5513e-04
% fit to teng
% xhat =
%    0.6668   -0.3440  -14.4332    0.8034   -1.6546

%%
[par, C2P] = C2Puptake(xhat,par);
C2P(C2P==0) = NaN;

GridVer = 90;
if GridVer == 90;
    lat = [-89:2:89];
    lon = [1:2:360];
end
%% C2P surface plot
%Zlevs = [100:20:300];
figure;
contourf(lon,lat,C2P(:,:,1)); hold on
cb=colorbar;
colormap(flipud(summer))
title('Cell Model C:P Uptake - Fit to Teng: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'C:P [gC/gP]');
axis tight; grid off

figTitle = 'C2Psurface';
%print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%figure; hold on; contourf(1:180,-89:2:89,C2P(:,:,1)); colorbar; title('surface C2P - fit to Teng 2014')
%%
function [f,fx,fxx] = C2Pneglogpost(x, par,iter)
on = true; off = false;
global iter
    M3d = par.M3d;
    dVt = par.dVt;
    iprod = par.iprod;
    
    if iter == 0
        printpara(x, par);
    end
    
    % reset para
    if iter < 4
        load(par.ftest,'x0')
        if (par.opt_fRibE == on)
	        ifRibE = par.pindx.lfRibE      ;
	        xnew   = exp(x(ifRibE))    ;
	        xold   = exp(x0(ifRibE))   ;
	        if (xnew > 1 | xnew <= 0)
	            % x(isigma) = x0(isigma) ;
	            x(ifRibE) = log(0.3+rand*0.2) ;
            end
        end
    end
    
    printpara(x, par);
    
    iter = iter + 1;    
    
    [par, C2P, C2Px, C2Pxx] = C2Puptake(x,par);
    %[par, C2P] = C2Puptake(x,par);
    
    % quick test - difference from CP=120
     Wb   = d0(dVt(iprod)/sum(dVt(iprod))) ;
     %eb = C2P(iprod) - 120.*ones(size(iprod));
     eb = C2P(iprod) - par.CPTeng(iprod);
     
     f = 0.5*(eb.'*Wb*eb);
     %f  = f + 0.5*(eb.'*Wb*eb) ;
        
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% calculate gradient
 if (nargout > 1)
        fx = zeros(length(x), 1)   ;

        nbx = par.nbx              ;
        % ---------------------------------
        for ji = 1 : nbx
            %c2px = C2Px(1:nwet,:);
            %fx(ji) = eb.'*Wb*C2Px(iteng,ji);
            fx(ji) = eb.'*Wb*C2Px(:,ji);
        end    
 end
    
% calculate hessian
 if nargout>2
        fxx = sparse(nbx, nbx)  ;

        % ----------------------------------------------------------------
        kk = 0;
        for ju = 1:nbx
            for jo = ju:nbx
                kk = kk + 1 ;
                %c2pxx = C2Pxx(1:nwet,:);
                
                fxx(ju,jo) = C2Px(:,ju).'*Wb*C2Px(:,jo) + eb.'*Wb*C2Pxx(:,kk);

                % make Hessian symetric;
                fxx(jo, ju) = fxx(ju, jo);
            end
        end
        %kpp = kk;
 end
    

end



function [par, C2P, C2Px, C2Pxx] = C2Puptake(x,par)
on = true; off = false;
    M3d =par.M3d;
    %pindx = par.pindx;
    iprod = par.iprod;
    iwet = par.iwet;
    nwet = par.nwet;
    
    P0=par.P0;
    N0=par.N0;
    T0=par.T0;
    Irr0=par.Irr0;

		[CellOut, parBIO] = CellCNP(par,x, P0,N0,T0,Irr0);
		par.BIO = parBIO;
		clear parBIO;
		par.CellOut.C2P = M3d*0;
		par.CellOut.N2P = M3d*0;
		par.CellOut.C2N = M3d*0;
		par.CellOut.LimType = NaN(size(M3d));
		par.CellOut.r = M3d*0;

		par.CellOut.C2P(iprod) = CellOut.CP;
		par.CellOut.N2P(iprod) = CellOut.NP;
		par.CellOut.C2N(iprod) = CellOut.CN;
		par.CellOut.LimType(iprod) = CellOut.LimType;
		par.CellOut.r(iprod) = CellOut.r;

		par.CellOut.C2P(isnan(par.CellOut.C2P)) = 0; %remove NaNs
        
        C2P = par.CellOut.C2P;
        
    if par.optim == off
        C2Px = [];
    elseif (par.optim & nargout > 2)
        % gradient of uptake operator
        nbx  = par.nbx; %ncx=par.ncx; npx =par.npx;
        C2Px  = zeros(nwet,nbx); % C2Px  = zeros(nwet,npx+ncx+nbx);

        if (par.opt_Q10Photo)
            dC2P_dQ10Photo = M3d*0;
			dC2P_dQ10Photo(iprod) = CellOut.dC2P_dQ10Photo;
            C2Px(:,par.pindx.lQ10Photo) = dC2P_dQ10Photo(iwet);
        end
        if (par.opt_fStorage)
            dC2P_dfStorage = M3d*0;
			dC2P_dfStorage(iprod) = CellOut.dC2P_dfStorage;
            C2Px(:,par.pindx.lfStorage) = dC2P_dfStorage(iwet);
        end
        if par.opt_fRibE
            dC2P_dfRibE  = M3d*0;
			dC2P_dfRibE(iprod)  = CellOut.dC2P_dfRibE;
            C2Px(:,par.pindx.lfRibE) = dC2P_dfRibE(iwet);
        end
		if (par.opt_PLip_PCutoff)
            dC2P_dPCutoff  = M3d*0;
			dC2P_dPCutoff(iprod)  = CellOut.dC2P_dPCutoff;
			C2Px(:,par.pindx.lPLip_PCutoff) = dC2P_dPCutoff(iwet);
        end
        if (par.opt_PLip_scale)
            dC2P_dPLip_scale = M3d*0;
			dC2P_dPLip_scale(iprod) = CellOut.dC2P_dPLipscale;
			C2Px(:,par.pindx.lPLip_scale) = dC2P_dPLip_scale(iwet);
		end
		if (par.opt_PStor_rCutoff)
            dC2P_drCutoff  = M3d*0;
			dC2P_drCutoff(iprod)  = CellOut.dC2P_drCutoff;
			C2Px(:,par.pindx.lPStor_rCutoff) = dC2P_drCutoff(iwet);
        end
        if par.opt_PStor_scale
            dC2P_dPStor_scale = M3d*0;
			dC2P_dPStor_scale(iprod) = CellOut.dC2P_dPStorscale;
			C2Px(:,par.pindx.lPStor_scale) = dC2P_dPStor_scale(iwet);
        end
        if par.opt_alphaS
            dC2P_dalphaS = M3d*0;
			dC2P_dalphaS(iprod) = CellOut.dC2P_dalphaS;
            C2Px(:,par.pindx.lalphaS) =  dC2P_dalphaS(iwet);
        end
    end   

        
 %%%%%% hessian v2 %%%%%%%%%
    if par.optim == off
        C2Pxx = [];
    elseif (par.optim & nargout > 3)
        kk = 0;
        
        % Q10Photo Derivatives
		if (par.opt_Q10Photo)
            kk = kk + 1;

			%second Derivatives w.r.t. Q10Photo
			xim = zeros(size(x));
			xim(par.pindx.lQ10Photo) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dQ10Photo2 = M3d*0;
			d2C2P_dQ10Photo2(iprod) = imag(CellOut.dC2P_dQ10Photo)./eps^3;
            C2Pxx(:,kk) = d2C2P_dQ10Photo2(iwet);

			if (par.opt_fStorage)
                kk = kk + 1;
				d2C2P_dfStorage_dQ10Photo = M3d*0;
				d2C2P_dfStorage_dQ10Photo(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3;
                C2Pxx(:,kk) = d2C2P_dfStorage_dQ10Photo(iwet);
            end
            
            if (par.opt_fRibE)
                kk = kk + 1;
                d2C2P_dfRibE_dQ10Photo = M3d*0;
				d2C2P_dfRibE_dQ10Photo(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
                C2Pxx(:,kk) = d2C2P_dfRibE_dQ10Photo(iwet);
            end
            
			if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dQ10Photo = M3d*0;
				d2C2P_dPCutoff_dQ10Photo(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPCutoff_dQ10Photo(iwet);
            end
            
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dQ10Photo = M3d*0;
				d2C2P_drCutoff_dQ10Photo(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dQ10Photo(iwet);
            end
            
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dQ10Photo = M3d*0;
				d2C2P_dPStorscale_dQ10Photo(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dQ10Photo(iwet);
            end
            
			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dQ10Photo = M3d*0;
				d2C2P_dPLipscale_dQ10Photo(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPLipscale_dQ10Photo(iwet);
            end
            
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dQ10Photo = M3d*0;
				d2C2P_dalphaS_dQ10Photo(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dQ10Photo(iwet);
			end
        end
        
        % fStorage
		if (par.opt_fStorage)
			kk = kk + 1;
            
			%second Derivatives w.r.t. fStorage
			xim = zeros(size(x));
			xim(par.pindx.lfStorage) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dfStorage2 = M3d*0;
			d2C2P_dfStorage2(iprod) = imag(CellOut.dC2P_dfStorage)./eps^3;
            C2Pxx(:,kk) = d2C2P_dfStorage2(iwet);

            if (par.opt_fRibE)
                kk = kk + 1;
                d2C2P_dfRibE_dfStorage = M3d*0;
				d2C2P_dfRibE_dfStorage(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
                C2Pxx(:,kk) = d2C2P_dfRibE_dfStorage(iwet);
            end
			if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dfStorage = M3d*0;
				d2C2P_dPCutoff_dfStorage(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPCutoff_dfStorage(iwet);
			end
			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dfStorage = M3d*0;
				d2C2P_dPLipscale_dfStorage(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPLipscale_dfStorage(iwet);
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dfStorage = M3d*0;
				d2C2P_drCutoff_dfStorage(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dfStorage(iwet);
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dfStorage = M3d*0;
				d2C2P_dPStorscale_dfStorage(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dfStorage(iwet);
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dfStorage = M3d*0;
				d2C2P_dalphaS_dfStorage(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dfStorage(iwet);
			end
        end
        
        %fRibE derivatives
		if (par.opt_fRibE)
            kk = kk + 1;

			%second Derivatives w.r.t. PLip_PCutoff
			xim = zeros(size(x));
			xim(par.pindx.lfRibE) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dfRibE2 = M3d*0;
			d2C2P_dfRibE2(iprod) = imag(CellOut.dC2P_dfRibE)./eps^3;
            C2Pxx(:,kk) = d2C2P_dfRibE2(iwet);

            if (par.opt_PLip_PCutoff)
                kk = kk + 1;
				d2C2P_dPCutoff_dfRibE = M3d*0;
				d2C2P_dPCutoff_dfRibE(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPCutoff_dfRibE(iwet);
			end
			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dfRibE = M3d*0;
				d2C2P_dPLipscale_dfRibE(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPLipscale_dfRibE(iwet);
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dfRibE = M3d*0;
				d2C2P_drCutoff_dfRibE(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dfRibE(iwet);
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dfRibE = M3d*0;
				d2C2P_dPStorscale_dfRibE(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dfRibE(iwet);
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dfRibE = M3d*0;
				d2C2P_dalphaS_dfRibE(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dfRibE(iwet);
			end
        end
        
        %PLip_PCutoff
        % PLip_PCutoff derivatives
        if (par.opt_PLip_PCutoff)
			kk = kk + 1;
            
			xim = zeros(size(x));
			xim(par.pindx.lPLip_PCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dPCutoff2 = M3d*0;
			d2C2P_dPCutoff2(iprod) = imag(CellOut.dC2P_dPCutoff)./eps^3;
            C2Pxx(:,kk) = d2C2P_dPCutoff2(iwet);

			if (par.opt_PLip_scale)
                kk = kk + 1;
				d2C2P_dPLipscale_dPCutoff = M3d*0;
				d2C2P_dPLipscale_dPCutoff(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPLipscale_dPCutoff(iwet);
			end
			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dPCutoff = M3d*0;
				d2C2P_drCutoff_dPCutoff(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dPCutoff(iwet);
			end
			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_dPCutoff = M3d*0;
				d2C2P_dPStorscale_dPCutoff(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dPCutoff(iwet);
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dPCutoff = M3d*0;
				d2C2P_dalphaS_dPCutoff(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dPCutoff(iwet);
            end
        end
        
        %PLip_scale
        if (par.opt_PLip_scale)
			kk = kk + 1;
            
			xim = zeros(size(x));
			xim(par.pindx.lPLip_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dPLipscale2 = M3d*0;
			d2C2P_dPLipscale2(iprod) = imag(CellOut.dC2P_dPLipscale)./eps^3;
            C2Pxx(:,kk) = d2C2P_dPLipscale2(iwet);

			if (par.opt_PStor_rCutoff)
                kk = kk + 1;
				d2C2P_drCutoff_dPLipscale = M3d*0;
				d2C2P_drCutoff_dPLipscale(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
                C2Pxx(:,kk) = d2C2P_drCutoff_dPLipscale(iwet);
			end
			if (par.opt_PStor_scale)
				kk = kk + 1;
                d2C2P_dPStorscale_dPLipscale = M3d*0;
				d2C2P_dPStorscale_dPLipscale(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_dPLipscale(iwet);
			end
			if (par.opt_alphaS)
				kk = kk + 1;
                d2C2P_dalphaS_dPLipscale = M3d*0;
				d2C2P_dalphaS_dPLipscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dPLipscale(iwet);
            end
        end
        
        %PStor_rCutoff
        if (par.opt_PStor_rCutoff)
			kk = kk + 1;
            
			xim = zeros(size(x));
			xim(par.pindx.lPStor_rCutoff) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_drCutoff2 = M3d*0;
			d2C2P_drCutoff2(iprod) = imag(CellOut.dC2P_drCutoff)./eps^3;
            C2Pxx(:,kk) = d2C2P_drCutoff2(iwet);

			if (par.opt_PStor_scale)
                kk = kk + 1;
				d2C2P_dPStorscale_drCutoff = M3d*0;
				d2C2P_dPStorscale_drCutoff(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
                C2Pxx(:,kk) = d2C2P_dPStorscale_drCutoff(iwet);
			end
			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_drCutoff = M3d*0;
				d2C2P_dalphaS_drCutoff(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_drCutoff(iwet);
            end
        end
        
        %PStor_scale
        if (par.opt_PStor_scale)
			kk = kk + 1;
            
			xim = zeros(size(x));
			xim(par.pindx.lPStor_scale) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dPStorscale2 = M3d*0;
			d2C2P_dPStorscale2(iprod) = imag(CellOut.dC2P_dPStorscale)./eps^3;
            C2Pxx(:,kk) = d2C2P_dPStorscale2(iwet);

			if (par.opt_alphaS)
                kk = kk + 1;
				d2C2P_dalphaS_dPStorscale = M3d*0;
				d2C2P_dalphaS_dPStorscale(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
                C2Pxx(:,kk) = d2C2P_dalphaS_dPStorscale(iwet);
			end
        end
        
        % alphaS
        if (par.opt_alphaS)
			kk = kk + 1;
            
			xim = zeros(size(x));
			xim(par.pindx.lalphaS) = sqrt(-1)*eps^3;
			[CellOut, ~] = CellCNP(par,x+xim, P0,N0,T0,Irr0);

			d2C2P_dalphaS2 = M3d*0;
			d2C2P_dalphaS2(iprod) = imag(CellOut.dC2P_dalphaS)./eps^3;
            C2Pxx(:,kk) = d2C2P_dalphaS2(iwet);
        end        
        
    end    
end

function printpara(x, par)
    on = true; off = false;
        if (par.opt_Q10Photo == on)
            iQ10Photo = par.pindx.lQ10Photo;
            fprintf('current Q10Photo       is  % 3.2e \n', exp(x(iQ10Photo)));
            xhat.Q10Photo = exp(x(iQ10Photo));
        end
		if (par.opt_fStorage == on)
            ifStorage = par.pindx.lfStorage;
            fprintf('current fStorage       is  % 3.2e \n', exp(x(ifStorage)));
            xhat.fStorage = exp(x(ifStorage));
        end
        if (par.opt_fRibE == on)
            ifRibE = par.pindx.lfRibE;
            fprintf('current fRibE         is  % 3.2e \n', exp(x(ifRibE)));
            xhat.fRibE = exp(x(ifRibE));
        end
		if (par.opt_PLip_PCutoff == on)
            iPLip_PCutoff = par.pindx.lPLip_PCutoff;
            fprintf('current PLip_PCutoff   is  % 3.2e \n', exp(x(iPLip_PCutoff)));
            xhat.PLip_PCutoff = exp(x(iPLip_PCutoff));
        end
		if (par.opt_PLip_scale == on)
            iPLip_scale = par.pindx.lPLip_scale;
            fprintf('current PLip_scale     is  % 3.2e \n', exp(x(iPLip_scale)));
            xhat.PLip_scale = exp(x(iPLip_scale));
        end
		if (par.opt_PStor_rCutoff == on)
            iPStor_rCutoff = par.pindx.lPStor_rCutoff;
            fprintf('current PStor_rCutoff  is  % 3.2e \n', exp(x(iPStor_rCutoff)));
            xhat.PStor_rCutoff = exp(x(iPStor_rCutoff));
        end
		if (par.opt_PStor_scale == on)
            iPStor_scale = par.pindx.lPStor_scale;
            fprintf('current PStor_scale    is  % 3.2e \n', exp(x(iPStor_scale)));
            xhat.PStor_scale = exp(x(iPStor_scale));
        end
		if (par.opt_alphaS == on)
            ialphaS = par.pindx.lalphaS;
            fprintf('current alphaS         is  % 3.2e \n', exp(x(ialphaS)));
            xhat.alphaS = exp(x(ialphaS));
        end
    
    x0 = x ;
    if (par.optim == on)
        save(par.ftest, 'x0','xhat')
    end
    
end